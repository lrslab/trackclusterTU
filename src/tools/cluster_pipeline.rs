use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use clap::{Parser, ValueEnum};
use rayon::prelude::*;

use crate::model::{Bed12Attrs, Coord, Interval, Strand, Transcript};
use crate::tu::multi::{
    parse_manifest, pooled_read_id, MultiSampleCounter, MultiSampleCounts, SampleManifestRecord,
};
use crate::tu::{ReadRecord, Tu, TuClusteringOptions};

#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
enum InputFormat {
    Auto,
    Bed6,
    Bed12,
    Tsv,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ResolvedInputFormat {
    Bed6,
    Bed12,
    Tsv,
    Fastq,
}

#[derive(Parser, Debug, Clone)]
#[command(
    name = "trackclustertu cluster",
    version,
    about = "Cluster bacterial directRNA reads into transcript units (TUs)"
)]
struct ClusterCli {
    /// Single input reads file path.
    #[arg(long = "in", conflicts_with = "manifest")]
    input: Option<PathBuf>,

    /// Multi-sample manifest TSV with columns: sample, reads, [group].
    #[arg(long, conflicts_with = "input")]
    manifest: Option<PathBuf>,

    /// Input format. Defaults to auto-detecting from file extensions.
    #[arg(long, value_enum, default_value_t = InputFormat::Auto)]
    format: InputFormat,

    /// Output directory. Missing output paths default to files inside this directory.
    #[arg(long)]
    out_dir: Option<PathBuf>,

    /// score1 threshold (overlap / union).
    #[arg(long, default_value_t = 0.95)]
    score1_threshold: f64,

    /// score2 threshold (overlap / max_len).
    #[arg(long, default_value_t = 0.6)]
    score2_threshold: f64,

    /// Allowed strand-aware 3 prime mismatch during second-pass score2 attachment (bp).
    #[arg(long, default_value_t = 12)]
    three_prime_tolerance_bp: u32,

    /// Optional maximum strand-aware 5 prime delta allowed during second-pass attachment (bp).
    ///
    /// When set, pairs within this 5 prime cap and the 3 prime tolerance may still merge even
    /// if score2 falls below the score2 threshold.
    #[arg(long = "max-5p-delta")]
    max_five_prime_delta_bp: Option<u32>,

    /// Skip the second-pass score2 merge and keep score1 seed clusters as final TUs.
    #[arg(long)]
    skip_score2_attachment: bool,

    /// Optional minimum read length filter (bp).
    #[arg(long)]
    min_read_len: Option<u32>,

    /// Output shared TU BED6.
    #[arg(long)]
    out_tu: Option<PathBuf>,

    /// Output pooled membership TSV.
    #[arg(long)]
    out_membership: Option<PathBuf>,

    /// Optional pooled BED6 used for clustering after sample-tagging read IDs.
    #[arg(long, requires = "manifest")]
    out_pooled_reads: Option<PathBuf>,

    /// Optional TU expression CSV (`tu_id,count`).
    #[arg(long)]
    out_tu_count: Option<PathBuf>,

    /// Optional per-sample TU long table (`tu_id, sample, count` as TSV).
    #[arg(long, requires = "manifest")]
    out_tu_sample_count_long: Option<PathBuf>,

    /// Optional per-sample TU matrix (`tu_id` + manifest-order sample columns as TSV).
    #[arg(long, requires = "manifest")]
    out_tu_sample_count_matrix: Option<PathBuf>,

    /// Optional per-group TU matrix (`tu_id` + group columns as TSV).
    ///
    /// This file is only written when the manifest has a non-empty `group` column.
    #[arg(long, requires = "manifest")]
    out_tu_group_count_matrix: Option<PathBuf>,

    /// Optional minimum reads per TU (filters outputs).
    #[arg(long)]
    min_tu_count: Option<u64>,

    /// Optional annotation BED6 (e.g., genes) used to anchor TU results.
    #[arg(long)]
    annotation_bed: Option<PathBuf>,

    /// Optional TU-to-gene overlap TSV (one line per TU×gene overlap).
    #[arg(long, requires = "annotation_bed")]
    out_tu_gene: Option<PathBuf>,

    /// Optional TU BED12 output, with blocks clipped to overlapping annotated genes.
    ///
    /// Adds extra columns:
    /// - name2: comma-separated subreads, with `|<read_count>` suffix (TrackCluster-style)
    /// - gene_list: comma-separated overlapping genes (or '.')
    #[arg(long, requires = "annotation_bed")]
    out_tu_bed12: Option<PathBuf>,

    /// Optional gene counts CSV (`gene_id,count`) from TU-to-gene overlaps.
    #[arg(long, requires = "annotation_bed")]
    out_gene_count: Option<PathBuf>,

    /// Optional per-sample gene count matrix (`gene_id` + manifest-order sample columns as TSV).
    #[arg(long, requires_all = ["annotation_bed", "manifest"])]
    out_gene_sample_count_matrix: Option<PathBuf>,

    /// Optional per-group gene count matrix (`gene_id` + group columns as TSV).
    #[arg(long, requires_all = ["annotation_bed", "manifest"])]
    out_gene_group_count_matrix: Option<PathBuf>,

    /// Number of worker threads to use (default: all logical CPUs).
    #[arg(long)]
    threads: Option<usize>,

    /// Print a timing breakdown to stderr.
    #[arg(long)]
    timings: bool,
}

#[derive(Parser, Debug, Clone)]
#[command(
    name = "trackclustertu recount",
    version,
    about = "Recompute TU count tables from pooled membership assignments"
)]
struct RecountCli {
    /// Multi-sample manifest TSV with columns: sample, reads, [group].
    #[arg(long)]
    manifest: PathBuf,

    /// Existing pooled membership TSV used to recompute multi-sample counts without reclustering.
    #[arg(long)]
    pooled_membership: PathBuf,

    /// Output directory. Missing output paths default to files inside this directory.
    #[arg(long)]
    out_dir: Option<PathBuf>,

    /// Optional TU expression CSV (`tu_id,count`).
    #[arg(long)]
    out_tu_count: Option<PathBuf>,

    /// Optional per-sample TU long table (`tu_id, sample, count` as TSV).
    #[arg(long)]
    out_tu_sample_count_long: Option<PathBuf>,

    /// Optional per-sample TU matrix (`tu_id` + manifest-order sample columns as TSV).
    #[arg(long)]
    out_tu_sample_count_matrix: Option<PathBuf>,

    /// Optional per-group TU matrix (`tu_id` + group columns as TSV).
    ///
    /// This file is only written when the manifest has a non-empty `group` column.
    #[arg(long)]
    out_tu_group_count_matrix: Option<PathBuf>,

    /// Optional minimum reads per TU (filters outputs).
    #[arg(long)]
    min_tu_count: Option<u64>,

    /// Number of worker threads to use (default: all logical CPUs).
    #[arg(long)]
    threads: Option<usize>,

    /// Print a timing breakdown to stderr.
    #[arg(long)]
    timings: bool,
}

#[derive(Debug, Clone)]
struct Cli {
    input: Option<PathBuf>,
    manifest: Option<PathBuf>,
    pooled_membership: Option<PathBuf>,
    format: InputFormat,
    out_dir: Option<PathBuf>,
    score1_threshold: f64,
    score2_threshold: f64,
    three_prime_tolerance_bp: u32,
    max_five_prime_delta_bp: Option<u32>,
    skip_score2_attachment: bool,
    min_read_len: Option<u32>,
    out_tu: Option<PathBuf>,
    out_membership: Option<PathBuf>,
    out_pooled_reads: Option<PathBuf>,
    out_tu_count: Option<PathBuf>,
    out_tu_sample_count_long: Option<PathBuf>,
    out_tu_sample_count_matrix: Option<PathBuf>,
    out_tu_group_count_matrix: Option<PathBuf>,
    min_tu_count: Option<u64>,
    annotation_bed: Option<PathBuf>,
    out_tu_gene: Option<PathBuf>,
    out_tu_bed12: Option<PathBuf>,
    out_gene_count: Option<PathBuf>,
    out_gene_sample_count_matrix: Option<PathBuf>,
    out_gene_group_count_matrix: Option<PathBuf>,
    threads: Option<usize>,
    timings: bool,
}

impl From<ClusterCli> for Cli {
    fn from(cli: ClusterCli) -> Self {
        Self {
            input: cli.input,
            manifest: cli.manifest,
            pooled_membership: None,
            format: cli.format,
            out_dir: cli.out_dir,
            score1_threshold: cli.score1_threshold,
            score2_threshold: cli.score2_threshold,
            three_prime_tolerance_bp: cli.three_prime_tolerance_bp,
            max_five_prime_delta_bp: cli.max_five_prime_delta_bp,
            skip_score2_attachment: cli.skip_score2_attachment,
            min_read_len: cli.min_read_len,
            out_tu: cli.out_tu,
            out_membership: cli.out_membership,
            out_pooled_reads: cli.out_pooled_reads,
            out_tu_count: cli.out_tu_count,
            out_tu_sample_count_long: cli.out_tu_sample_count_long,
            out_tu_sample_count_matrix: cli.out_tu_sample_count_matrix,
            out_tu_group_count_matrix: cli.out_tu_group_count_matrix,
            min_tu_count: cli.min_tu_count,
            annotation_bed: cli.annotation_bed,
            out_tu_gene: cli.out_tu_gene,
            out_tu_bed12: cli.out_tu_bed12,
            out_gene_count: cli.out_gene_count,
            out_gene_sample_count_matrix: cli.out_gene_sample_count_matrix,
            out_gene_group_count_matrix: cli.out_gene_group_count_matrix,
            threads: cli.threads,
            timings: cli.timings,
        }
    }
}

impl From<RecountCli> for Cli {
    fn from(cli: RecountCli) -> Self {
        Self {
            input: None,
            manifest: Some(cli.manifest),
            pooled_membership: Some(cli.pooled_membership),
            format: InputFormat::Auto,
            out_dir: cli.out_dir,
            score1_threshold: 0.95,
            score2_threshold: 0.6,
            three_prime_tolerance_bp: 12,
            max_five_prime_delta_bp: None,
            skip_score2_attachment: false,
            min_read_len: None,
            out_tu: None,
            out_membership: None,
            out_pooled_reads: None,
            out_tu_count: cli.out_tu_count,
            out_tu_sample_count_long: cli.out_tu_sample_count_long,
            out_tu_sample_count_matrix: cli.out_tu_sample_count_matrix,
            out_tu_group_count_matrix: cli.out_tu_group_count_matrix,
            min_tu_count: cli.min_tu_count,
            annotation_bed: None,
            out_tu_gene: None,
            out_tu_bed12: None,
            out_gene_count: None,
            out_gene_sample_count_matrix: None,
            out_gene_group_count_matrix: None,
            threads: cli.threads,
            timings: cli.timings,
        }
    }
}

#[derive(Clone, Debug)]
struct GeneRecord {
    contig: String,
    strand: Strand,
    interval: Interval,
    id: String,
}

#[derive(Clone, Debug)]
struct PreparedOutputs {
    tus: Vec<Tu>,
    old_to_new: Vec<Option<usize>>,
    counts_kept: Vec<u64>,
}

enum InputMode {
    Single,
    ManifestCluster {
        samples: Vec<SampleManifestRecord>,
    },
    ManifestCountOnly {
        samples: Vec<SampleManifestRecord>,
        membership: PathBuf,
    },
}

fn default_output_dir(cli: &Cli, mode: &InputMode) -> PathBuf {
    if let Some(out_dir) = cli.out_dir.as_ref() {
        return out_dir.clone();
    }

    let source_path = match mode {
        InputMode::Single => cli
            .input
            .as_ref()
            .expect("single-input mode must have --in"),
        InputMode::ManifestCluster { .. } | InputMode::ManifestCountOnly { .. } => cli
            .manifest
            .as_ref()
            .expect("manifest modes must have --manifest"),
    };

    let stem = source_path
        .file_stem()
        .and_then(|stem| stem.to_str())
        .filter(|stem| !stem.is_empty())
        .unwrap_or("trackclustertu");
    source_path
        .parent()
        .unwrap_or(Path::new("."))
        .join(format!("{stem}.trackclustertu"))
}

fn set_default_path(slot: &mut Option<PathBuf>, out_dir: &Path, filename: &str) {
    if slot.is_none() {
        *slot = Some(out_dir.join(filename));
    }
}

fn apply_default_outputs(cli: &mut Cli, mode: &InputMode) {
    let out_dir = default_output_dir(cli, mode);
    cli.out_dir = Some(out_dir.clone());

    match mode {
        InputMode::Single => {
            set_default_path(&mut cli.out_tu, &out_dir, "tus.bed");
            set_default_path(&mut cli.out_membership, &out_dir, "membership.tsv");
            set_default_path(&mut cli.out_tu_count, &out_dir, "tu_count.csv");

            if cli.annotation_bed.is_some() {
                set_default_path(&mut cli.out_tu_gene, &out_dir, "tu_gene.tsv");
                set_default_path(&mut cli.out_tu_bed12, &out_dir, "tus.anchored.bed12");
                set_default_path(&mut cli.out_gene_count, &out_dir, "gene_count.csv");
            }
        }
        InputMode::ManifestCluster { .. } => {
            set_default_path(&mut cli.out_tu, &out_dir, "tus.bed");
            set_default_path(&mut cli.out_membership, &out_dir, "membership.tsv");
            set_default_path(&mut cli.out_pooled_reads, &out_dir, "pooled.bed");
            set_default_path(&mut cli.out_tu_count, &out_dir, "tu_count.csv");
            set_default_path(
                &mut cli.out_tu_sample_count_long,
                &out_dir,
                "tu_sample_long.tsv",
            );
            set_default_path(
                &mut cli.out_tu_sample_count_matrix,
                &out_dir,
                "tu_sample_matrix.tsv",
            );
            set_default_path(
                &mut cli.out_tu_group_count_matrix,
                &out_dir,
                "tu_group_matrix.tsv",
            );

            if cli.annotation_bed.is_some() {
                set_default_path(&mut cli.out_tu_gene, &out_dir, "tu_gene.tsv");
                set_default_path(&mut cli.out_tu_bed12, &out_dir, "tus.anchored.bed12");
                set_default_path(&mut cli.out_gene_count, &out_dir, "gene_count.csv");
                set_default_path(
                    &mut cli.out_gene_sample_count_matrix,
                    &out_dir,
                    "gene_sample_matrix.tsv",
                );
                set_default_path(
                    &mut cli.out_gene_group_count_matrix,
                    &out_dir,
                    "gene_group_matrix.tsv",
                );
            }
        }
        InputMode::ManifestCountOnly { .. } => {
            set_default_path(&mut cli.out_tu_count, &out_dir, "tu_count.csv");
            set_default_path(
                &mut cli.out_tu_sample_count_long,
                &out_dir,
                "tu_sample_long.tsv",
            );
            set_default_path(
                &mut cli.out_tu_sample_count_matrix,
                &out_dir,
                "tu_sample_matrix.tsv",
            );
            set_default_path(
                &mut cli.out_tu_group_count_matrix,
                &out_dir,
                "tu_group_matrix.tsv",
            );
        }
    }
}

fn prepare_output_dirs(cli: &Cli) -> anyhow::Result<()> {
    if let Some(out_dir) = cli.out_dir.as_ref() {
        fs::create_dir_all(out_dir)?;
    }

    let file_slots = [
        cli.out_tu.as_ref(),
        cli.out_membership.as_ref(),
        cli.out_pooled_reads.as_ref(),
        cli.out_tu_count.as_ref(),
        cli.out_tu_sample_count_long.as_ref(),
        cli.out_tu_sample_count_matrix.as_ref(),
        cli.out_tu_group_count_matrix.as_ref(),
        cli.out_tu_gene.as_ref(),
        cli.out_tu_bed12.as_ref(),
        cli.out_gene_count.as_ref(),
        cli.out_gene_sample_count_matrix.as_ref(),
        cli.out_gene_group_count_matrix.as_ref(),
    ];

    for path in file_slots.into_iter().flatten() {
        fs::create_dir_all(path.parent().unwrap_or(Path::new(".")))?;
    }

    Ok(())
}

fn cmp_read(a: &ReadRecord, b: &ReadRecord) -> std::cmp::Ordering {
    a.contig
        .cmp(&b.contig)
        .then_with(|| a.strand.cmp(&b.strand))
        .then_with(|| a.interval.start.cmp(&b.interval.start))
        .then_with(|| a.interval.end.cmp(&b.interval.end))
        .then_with(|| a.id.cmp(&b.id))
}

fn read_bed6(path: &Path) -> anyhow::Result<Vec<ReadRecord>> {
    let reader = crate::io::bed::read_bed6(path)?;
    let mut reads: Vec<ReadRecord> = Vec::new();
    for result in reader {
        let record = result?;
        reads.push(ReadRecord {
            contig: record.chrom,
            strand: record.strand,
            interval: Interval::new(record.start, record.end)?,
            id: record.name,
        });
    }
    Ok(reads)
}

fn read_bed12(path: &Path) -> anyhow::Result<Vec<ReadRecord>> {
    let reader = crate::io::bed::read_bed12(path)?;
    let mut reads: Vec<ReadRecord> = Vec::new();
    for result in reader {
        let tx = result?;
        // TU clustering is intentionally span-based, so BED12 block structure is
        // collapsed to the outer transcript interval here.
        reads.push(ReadRecord {
            contig: tx.chrom,
            strand: tx.strand,
            interval: Interval::new(tx.tx_start, tx.tx_end)?,
            id: tx.name,
        });
    }
    Ok(reads)
}

fn read_tsv(path: &Path) -> anyhow::Result<Vec<ReadRecord>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut reads: Vec<ReadRecord> = Vec::new();

    for (line_number, line) in reader.lines().enumerate() {
        let line_number = line_number + 1;
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let mut fields: Vec<&str> = line.split('\t').collect();
        if fields.len() == 1 {
            fields = line.split_whitespace().collect();
        }

        if fields.len() < 5 {
            anyhow::bail!(
                "{path:?}:{line_number}: expected at least 5 columns (contig, start, end, id, strand), got {}",
                fields.len()
            );
        }

        let contig = fields[0].to_owned();
        let start: u32 = fields[1].parse().map_err(|_| {
            anyhow::anyhow!(
                "{path:?}:{line_number}: invalid integer for start: {:?}",
                fields[1]
            )
        })?;
        let end: u32 = fields[2].parse().map_err(|_| {
            anyhow::anyhow!(
                "{path:?}:{line_number}: invalid integer for end: {:?}",
                fields[2]
            )
        })?;
        let id = fields[3].to_owned();
        let strand = Strand::try_from(fields[4]).map_err(|e| {
            anyhow::anyhow!(
                "{path:?}:{line_number}: invalid strand {:?}: {e}",
                fields[4]
            )
        })?;

        reads.push(ReadRecord {
            contig,
            strand,
            interval: Interval::new(Coord::new(start), Coord::new(end))?,
            id,
        });
    }

    Ok(reads)
}

fn resolve_explicit_input_format(format: InputFormat) -> Option<ResolvedInputFormat> {
    match format {
        InputFormat::Auto => None,
        InputFormat::Bed6 => Some(ResolvedInputFormat::Bed6),
        InputFormat::Bed12 => Some(ResolvedInputFormat::Bed12),
        InputFormat::Tsv => Some(ResolvedInputFormat::Tsv),
    }
}

fn read_reads(path: &Path, format: ResolvedInputFormat) -> anyhow::Result<Vec<ReadRecord>> {
    match format {
        ResolvedInputFormat::Bed6 => read_bed6(path),
        ResolvedInputFormat::Bed12 => read_bed12(path),
        ResolvedInputFormat::Tsv => read_tsv(path),
        ResolvedInputFormat::Fastq => anyhow::bail!(
            "FASTQ inputs are not supported directly; use `trackclustertu map` or `trackclustertu run`, or convert sorted BAMs with `trackclustertu bam-to-bed` before clustering"
        ),
    }
}

fn infer_format_from_path(path: &Path) -> Option<ResolvedInputFormat> {
    let file_name = path.file_name()?.to_string_lossy().to_ascii_lowercase();
    if file_name.ends_with(".fastq")
        || file_name.ends_with(".fq")
        || file_name.ends_with(".fastq.gz")
        || file_name.ends_with(".fq.gz")
    {
        return Some(ResolvedInputFormat::Fastq);
    }
    if file_name.ends_with(".bed12") {
        return Some(ResolvedInputFormat::Bed12);
    }
    if file_name.ends_with(".bed") {
        return Some(ResolvedInputFormat::Bed6);
    }
    if file_name.ends_with(".tsv") || file_name.ends_with(".txt") {
        return Some(ResolvedInputFormat::Tsv);
    }
    None
}

fn resolve_single_input_format(cli: &Cli) -> anyhow::Result<ResolvedInputFormat> {
    if let Some(format) = resolve_explicit_input_format(cli.format) {
        return Ok(format);
    }

    let input = cli
        .input
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("single-input mode requires --in"))?;
    Ok(infer_format_from_path(input).unwrap_or(ResolvedInputFormat::Bed6))
}

fn resolve_manifest_input_format(
    cli: &Cli,
    samples: &[SampleManifestRecord],
) -> anyhow::Result<ResolvedInputFormat> {
    if let Some(format) = resolve_explicit_input_format(cli.format) {
        return Ok(format);
    }

    let mut inferred: Option<ResolvedInputFormat> = None;
    for sample in samples {
        let Some(sample_format) = infer_format_from_path(&sample.reads) else {
            continue;
        };

        if let Some(existing) = inferred {
            if existing != sample_format {
                anyhow::bail!(
                    "auto-detected mixed manifest input formats ({existing:?} and {sample_format:?}); pass --format explicitly"
                );
            }
        } else {
            inferred = Some(sample_format);
        }
    }

    Ok(inferred.unwrap_or(ResolvedInputFormat::Bed6))
}

fn read_annotation_bed6(path: &Path) -> anyhow::Result<Vec<GeneRecord>> {
    let reader = crate::io::bed::read_bed6(path)?;
    let mut genes: Vec<GeneRecord> = Vec::new();
    for result in reader {
        let record = result?;
        genes.push(GeneRecord {
            contig: record.chrom,
            strand: record.strand,
            interval: Interval::new(record.start, record.end)?,
            id: record.name,
        });
    }
    Ok(genes)
}

fn add_active(
    end: u32,
    idx: usize,
    active: &mut HashSet<usize>,
    ends: &mut BinaryHeap<Reverse<(u32, usize)>>,
) {
    active.insert(idx);
    ends.push(Reverse((end, idx)));
}

fn expire_active(
    current_start: u32,
    active: &mut HashSet<usize>,
    ends: &mut BinaryHeap<Reverse<(u32, usize)>>,
) {
    while let Some(Reverse((end, idx))) = ends.peek().copied() {
        if end <= current_start {
            ends.pop();
            active.remove(&idx);
        } else {
            break;
        }
    }
}

fn tu_gene_overlaps(tus: &[Tu], genes: &[GeneRecord]) -> Vec<Vec<(usize, u32)>> {
    let mut overlaps: Vec<Vec<(usize, u32)>> = vec![Vec::new(); tus.len()];

    let mut tus_by_key: BTreeMap<(String, Strand), Vec<usize>> = BTreeMap::new();
    for (idx, tu) in tus.iter().enumerate() {
        tus_by_key
            .entry((tu.contig.clone(), tu.strand))
            .or_default()
            .push(idx);
    }

    let mut genes_by_key: BTreeMap<(String, Strand), Vec<usize>> = BTreeMap::new();
    for (idx, gene) in genes.iter().enumerate() {
        genes_by_key
            .entry((gene.contig.clone(), gene.strand))
            .or_default()
            .push(idx);
    }

    for (key, mut tu_indices) in tus_by_key {
        let Some(gene_indices_ref) = genes_by_key.get(&key) else {
            continue;
        };
        let mut gene_indices = gene_indices_ref.clone();

        tu_indices.sort_by(|&a, &b| {
            tus[a]
                .interval
                .start
                .cmp(&tus[b].interval.start)
                .then_with(|| tus[a].interval.end.cmp(&tus[b].interval.end))
                .then_with(|| tus[a].id.cmp(&tus[b].id))
        });

        gene_indices.sort_by(|&a, &b| {
            genes[a]
                .interval
                .start
                .cmp(&genes[b].interval.start)
                .then_with(|| genes[a].interval.end.cmp(&genes[b].interval.end))
                .then_with(|| genes[a].id.cmp(&genes[b].id))
        });

        let mut ti: usize = 0;
        let mut gi: usize = 0;

        let mut active_tus: HashSet<usize> = HashSet::new();
        let mut active_genes: HashSet<usize> = HashSet::new();

        let mut tu_ends: BinaryHeap<Reverse<(u32, usize)>> = BinaryHeap::new();
        let mut gene_ends: BinaryHeap<Reverse<(u32, usize)>> = BinaryHeap::new();

        loop {
            let next_tu_start = tu_indices
                .get(ti)
                .map(|&idx| tus[idx].interval.start.get())
                .unwrap_or(u32::MAX);
            let next_gene_start = gene_indices
                .get(gi)
                .map(|&idx| genes[idx].interval.start.get())
                .unwrap_or(u32::MAX);

            if next_tu_start == u32::MAX && next_gene_start == u32::MAX {
                break;
            }

            let current_start = next_tu_start.min(next_gene_start);
            expire_active(current_start, &mut active_tus, &mut tu_ends);
            expire_active(current_start, &mut active_genes, &mut gene_ends);

            if next_tu_start <= next_gene_start {
                let tu_idx = tu_indices[ti];
                ti += 1;

                let tu = &tus[tu_idx];
                if tu.interval.is_empty() {
                    continue;
                }

                add_active(tu.interval.end.get(), tu_idx, &mut active_tus, &mut tu_ends);

                for &gene_idx in &active_genes {
                    let gene = &genes[gene_idx];
                    let overlap = tu.interval.overlap_len(gene.interval);
                    if overlap > 0 {
                        overlaps[tu_idx].push((gene_idx, overlap));
                    }
                }
            } else {
                let gene_idx = gene_indices[gi];
                gi += 1;

                let gene = &genes[gene_idx];
                if gene.interval.is_empty() {
                    continue;
                }

                add_active(
                    gene.interval.end.get(),
                    gene_idx,
                    &mut active_genes,
                    &mut gene_ends,
                );

                for &tu_idx in &active_tus {
                    let tu = &tus[tu_idx];
                    let overlap = tu.interval.overlap_len(gene.interval);
                    if overlap > 0 {
                        overlaps[tu_idx].push((gene_idx, overlap));
                    }
                }
            }
        }
    }

    for overlaps_for_tu in &mut overlaps {
        overlaps_for_tu.sort_by(|(gene_a, _), (gene_b, _)| {
            genes[*gene_a]
                .interval
                .start
                .cmp(&genes[*gene_b].interval.start)
                .then_with(|| {
                    genes[*gene_a]
                        .interval
                        .end
                        .cmp(&genes[*gene_b].interval.end)
                })
                .then_with(|| genes[*gene_a].id.cmp(&genes[*gene_b].id))
        });
        overlaps_for_tu.dedup_by(|(gene_a, _), (gene_b, _)| gene_a == gene_b);
    }

    overlaps
}

fn format_name2(read_indices: &[usize], reads: &[ReadRecord], read_count: usize) -> String {
    if read_indices.is_empty() {
        return format!(",|{read_count}");
    }

    let mut name2 = String::new();
    for (i, &read_idx) in read_indices.iter().enumerate() {
        if i > 0 {
            name2.push(',');
        }
        name2.push_str(&reads[read_idx].id);
    }
    name2.push_str(",|");
    name2.push_str(&read_count.to_string());
    name2
}

fn format_tu_id(index: usize, width: usize) -> String {
    format!("TU{:0width$}", index + 1, width = width)
}

fn default_threads() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

fn resolve_mode(cli: &Cli) -> anyhow::Result<InputMode> {
    match (&cli.input, &cli.manifest, &cli.pooled_membership) {
        (Some(_), None, None) => Ok(InputMode::Single),
        (None, Some(manifest), None) => Ok(InputMode::ManifestCluster {
            samples: parse_manifest(manifest)?,
        }),
        (None, Some(manifest), Some(membership)) => Ok(InputMode::ManifestCountOnly {
            samples: parse_manifest(manifest)?,
            membership: membership.clone(),
        }),
        (None, None, _) => anyhow::bail!("exactly one of --in or --manifest is required"),
        (Some(_), None, Some(_)) => anyhow::bail!("--pooled-membership cannot be used with --in"),
        (Some(_), Some(_), _) => anyhow::bail!("--in and --manifest are mutually exclusive"),
    }
}

fn want_multi_sample_outputs(cli: &Cli) -> bool {
    cli.out_tu_sample_count_long.is_some()
        || cli.out_tu_sample_count_matrix.is_some()
        || cli.out_tu_group_count_matrix.is_some()
}

fn want_gene_multi_sample_outputs(cli: &Cli) -> bool {
    cli.out_gene_sample_count_matrix.is_some() || cli.out_gene_group_count_matrix.is_some()
}

fn want_any_multi_sample_outputs(cli: &Cli) -> bool {
    want_multi_sample_outputs(cli) || want_gene_multi_sample_outputs(cli)
}

fn want_gene_outputs(cli: &Cli) -> bool {
    cli.out_gene_count.is_some() || want_gene_multi_sample_outputs(cli)
}

fn want_any_tu_count_outputs(cli: &Cli) -> bool {
    cli.out_tu_count.is_some() || want_multi_sample_outputs(cli)
}

fn load_manifest_reads(
    samples: &[SampleManifestRecord],
    format: ResolvedInputFormat,
) -> anyhow::Result<Vec<ReadRecord>> {
    let mut pooled_reads: Vec<ReadRecord> = Vec::new();

    for sample in samples {
        let mut reads = read_reads(&sample.reads, format)?;
        for read in &mut reads {
            read.id = pooled_read_id(&sample.sample, &read.id);
        }
        pooled_reads.extend(reads);
    }

    Ok(pooled_reads)
}

fn filter_reads(reads: &mut Vec<ReadRecord>, min_read_len: Option<u32>) {
    reads.retain(|read| !read.interval.is_empty());
    if let Some(min_len) = min_read_len {
        reads.retain(|read| read.interval.len() >= min_len);
    }
}

fn write_pooled_reads_bed6(path: &Path, reads: &[ReadRecord]) -> anyhow::Result<()> {
    let records: Vec<crate::io::bed::Bed6Record> = reads
        .iter()
        .map(|read| crate::io::bed::Bed6Record {
            chrom: read.contig.clone(),
            start: read.interval.start,
            end: read.interval.end,
            name: read.id.clone(),
            score: 0,
            strand: read.strand,
            extra_fields: Vec::new(),
        })
        .collect();
    crate::io::bed::write_bed6(path, records.iter())?;
    Ok(())
}

fn prepare_outputs(
    result: &crate::tu::TuClusteringResult,
    min_tu_count: Option<u64>,
) -> PreparedOutputs {
    let mut tu_counts: Vec<u64> = vec![0; result.tus.len()];
    for &tu_idx in &result.read_to_tu {
        tu_counts[tu_idx] += 1;
    }

    let mut kept_old_tu_indices: Vec<usize> = Vec::new();
    for (tu_idx, &count) in tu_counts.iter().enumerate() {
        if min_tu_count.is_none_or(|min| count >= min) {
            kept_old_tu_indices.push(tu_idx);
        }
    }

    let new_width = 6usize.max(kept_old_tu_indices.len().to_string().len());
    let mut old_to_new: Vec<Option<usize>> = vec![None; result.tus.len()];
    let mut tus: Vec<Tu> = Vec::with_capacity(kept_old_tu_indices.len());
    let mut counts_kept: Vec<u64> = Vec::with_capacity(kept_old_tu_indices.len());

    for (new_idx, &old_idx) in kept_old_tu_indices.iter().enumerate() {
        old_to_new[old_idx] = Some(new_idx);
        let tu = &result.tus[old_idx];
        tus.push(Tu {
            id: format_tu_id(new_idx, new_width),
            contig: tu.contig.clone(),
            strand: tu.strand,
            interval: tu.interval,
            rep_read_index: tu.rep_read_index,
        });
        counts_kept.push(tu_counts[old_idx]);
    }

    PreparedOutputs {
        tus,
        old_to_new,
        counts_kept,
    }
}

fn write_total_counts_csv(
    path: &Path,
    id_header: &str,
    ids: &[String],
    counts: &[u64],
) -> anyhow::Result<()> {
    let mut writer = std::io::BufWriter::new(File::create(path)?);
    writeln!(writer, "{id_header},count")?;
    for (id, count) in ids.iter().zip(counts.iter()) {
        writeln!(writer, "{id},{count}")?;
    }
    Ok(())
}

fn write_entity_matrix(
    path: &Path,
    id_header: &str,
    ids: &[String],
    column_names: &[String],
    counts: &[Vec<u64>],
) -> anyhow::Result<()> {
    let mut writer = std::io::BufWriter::new(File::create(path)?);
    write!(writer, "{id_header}")?;
    for column_name in column_names {
        write!(writer, "\t{column_name}")?;
    }
    writeln!(writer)?;

    for (row_idx, id) in ids.iter().enumerate() {
        write!(writer, "{id}")?;
        for count in &counts[row_idx] {
            write!(writer, "\t{count}")?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

fn write_entity_long_counts(
    path: &Path,
    id_header: &str,
    column_header: &str,
    ids: &[String],
    column_names: &[String],
    counts: &[Vec<u64>],
) -> anyhow::Result<()> {
    let mut writer = std::io::BufWriter::new(File::create(path)?);
    writeln!(writer, "{id_header}\t{column_header}\tcount")?;
    for (row_idx, id) in ids.iter().enumerate() {
        for (column_idx, column_name) in column_names.iter().enumerate() {
            let count = counts[row_idx][column_idx];
            if count > 0 {
                writeln!(writer, "{id}\t{column_name}\t{count}")?;
            }
        }
    }

    Ok(())
}

fn total_counts_from_multi(counts: &MultiSampleCounts) -> Vec<u64> {
    counts
        .counts
        .iter()
        .map(|row| row.iter().sum())
        .collect::<Vec<u64>>()
}

fn filter_multi_counts(counts: MultiSampleCounts, min_tu_count: Option<u64>) -> MultiSampleCounts {
    let Some(min_tu_count) = min_tu_count else {
        return counts;
    };

    let mut kept_tu_ids: Vec<String> = Vec::new();
    let mut kept_counts: Vec<Vec<u64>> = Vec::new();
    let mut kept_group_counts: Vec<Vec<u64>> = Vec::new();

    for tu_idx in 0..counts.tu_ids.len() {
        let total: u64 = counts.counts[tu_idx].iter().sum();
        if total >= min_tu_count {
            kept_tu_ids.push(counts.tu_ids[tu_idx].clone());
            kept_counts.push(counts.counts[tu_idx].clone());
            kept_group_counts.push(counts.group_counts[tu_idx].clone());
        }
    }

    MultiSampleCounts {
        tu_ids: kept_tu_ids,
        samples: counts.samples,
        counts: kept_counts,
        groups: counts.groups,
        group_counts: kept_group_counts,
    }
}

fn write_multi_sample_outputs(cli: &Cli, counts: &MultiSampleCounts) -> anyhow::Result<()> {
    if let Some(path) = cli.out_tu_sample_count_long.as_ref() {
        write_entity_long_counts(
            path,
            "tu_id",
            "sample",
            &counts.tu_ids,
            &counts.samples,
            &counts.counts,
        )?;
    }

    if let Some(path) = cli.out_tu_sample_count_matrix.as_ref() {
        write_entity_matrix(
            path,
            "tu_id",
            &counts.tu_ids,
            &counts.samples,
            &counts.counts,
        )?;
    }

    if let Some(path) = cli.out_tu_group_count_matrix.as_ref() {
        if counts.has_groups() {
            write_entity_matrix(
                path,
                "tu_id",
                &counts.tu_ids,
                &counts.groups,
                &counts.group_counts,
            )?;
        }
    }

    Ok(())
}

fn build_gene_total_counts(
    genes: &[GeneRecord],
    overlaps: &[Vec<(usize, u32)>],
    tu_counts: &[u64],
) -> Vec<u64> {
    let mut gene_counts = vec![0u64; genes.len()];
    for (tu_idx, gene_hits) in overlaps.iter().enumerate() {
        let mut seen = HashSet::new();
        for &(gene_idx, _) in gene_hits {
            if seen.insert(gene_idx) {
                gene_counts[gene_idx] += tu_counts[tu_idx];
            }
        }
    }
    gene_counts
}

fn build_gene_multi_counts(
    genes: &[GeneRecord],
    overlaps: &[Vec<(usize, u32)>],
    tu_counts: &MultiSampleCounts,
) -> (Vec<Vec<u64>>, Vec<Vec<u64>>) {
    let mut gene_sample_counts = vec![vec![0u64; tu_counts.samples.len()]; genes.len()];
    let mut gene_group_counts = vec![vec![0u64; tu_counts.groups.len()]; genes.len()];

    for (tu_idx, gene_hits) in overlaps.iter().enumerate() {
        let mut seen = HashSet::new();
        for &(gene_idx, _) in gene_hits {
            if !seen.insert(gene_idx) {
                continue;
            }
            for (sample_idx, count) in tu_counts.counts[tu_idx].iter().enumerate() {
                gene_sample_counts[gene_idx][sample_idx] += count;
            }
            for (group_idx, count) in tu_counts.group_counts[tu_idx].iter().enumerate() {
                gene_group_counts[gene_idx][group_idx] += count;
            }
        }
    }

    (gene_sample_counts, gene_group_counts)
}

fn write_gene_outputs(
    cli: &Cli,
    genes: &[GeneRecord],
    overlaps: &[Vec<(usize, u32)>],
    tu_counts: &[u64],
    multi_counts: Option<&MultiSampleCounts>,
) -> anyhow::Result<()> {
    let gene_ids: Vec<String> = genes.iter().map(|gene| gene.id.clone()).collect();

    if let Some(path) = cli.out_gene_count.as_ref() {
        let gene_total_counts = build_gene_total_counts(genes, overlaps, tu_counts);
        write_total_counts_csv(path, "gene_id", &gene_ids, &gene_total_counts)?;
    }

    if let Some(multi_counts) = multi_counts {
        let (gene_sample_counts, gene_group_counts) =
            build_gene_multi_counts(genes, overlaps, multi_counts);

        if let Some(path) = cli.out_gene_sample_count_matrix.as_ref() {
            write_entity_matrix(
                path,
                "gene_id",
                &gene_ids,
                &multi_counts.samples,
                &gene_sample_counts,
            )?;
        }

        if let Some(path) = cli.out_gene_group_count_matrix.as_ref() {
            if multi_counts.has_groups() {
                write_entity_matrix(
                    path,
                    "gene_id",
                    &gene_ids,
                    &multi_counts.groups,
                    &gene_group_counts,
                )?;
            }
        }
    }

    Ok(())
}

fn visit_membership_assignments<F>(path: &Path, mut visit: F) -> anyhow::Result<()>
where
    F: FnMut(&str, &str) -> anyhow::Result<()>,
{
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    for (line_idx, line_result) in reader.lines().enumerate() {
        let line_number = line_idx + 1;
        let line = line_result?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let mut fields: Vec<&str> = line.split('\t').collect();
        if fields.len() == 1 {
            fields = line.split_whitespace().collect();
        }
        if fields.len() < 2 {
            anyhow::bail!(
                "{path:?}:{line_number}: expected at least 2 columns (read_id, tu_id), got {}",
                fields.len()
            );
        }

        visit(fields[0], fields[1])?;
    }

    Ok(())
}

fn collect_membership_tu_ids(path: &Path) -> anyhow::Result<Vec<String>> {
    let mut tu_ids: BTreeSet<String> = BTreeSet::new();
    visit_membership_assignments(path, |_, tu_id| {
        tu_ids.insert(tu_id.to_owned());
        Ok(())
    })?;
    Ok(tu_ids.into_iter().collect())
}

fn write_empty_cluster_outputs(
    cli: &Cli,
    samples: Option<&[SampleManifestRecord]>,
) -> anyhow::Result<()> {
    let out_tu = cli
        .out_tu
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("--out-tu is required when clustering"))?;
    let out_membership = cli
        .out_membership
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("--out-membership is required when clustering"))?;

    File::create(out_tu)?;
    File::create(out_membership)?;

    if let Some(out_tu_count) = cli.out_tu_count.as_ref() {
        write_total_counts_csv(out_tu_count, "tu_id", &[], &[])?;
    }

    if let Some(out_tu_gene) = cli.out_tu_gene.as_ref() {
        let mut writer = std::io::BufWriter::new(File::create(out_tu_gene)?);
        writeln!(
            writer,
            "#contig\tstrand\ttu_id\ttu_start\ttu_end\tgene_id\tgene_start\tgene_end\toverlap_bp"
        )?;
    }

    if let Some(out_tu_bed12) = cli.out_tu_bed12.as_ref() {
        File::create(out_tu_bed12)?;
    }

    if let Some(out_pooled_reads) = cli.out_pooled_reads.as_ref() {
        File::create(out_pooled_reads)?;
    }

    if let Some(samples) = samples {
        let empty_counts = MultiSampleCounter::new(samples, Vec::new()).finish();
        write_multi_sample_outputs(cli, &empty_counts)?;
    }

    if want_gene_outputs(cli) {
        let genes = if let Some(annotation_path) = cli.annotation_bed.as_ref() {
            read_annotation_bed6(annotation_path)?
        } else {
            Vec::new()
        };
        let empty_overlaps: Vec<Vec<(usize, u32)>> = Vec::new();
        let empty_multi_counts =
            samples.map(|samples| MultiSampleCounter::new(samples, Vec::new()).finish());
        write_gene_outputs(
            cli,
            &genes,
            &empty_overlaps,
            &[],
            empty_multi_counts.as_ref(),
        )?;
    }

    Ok(())
}

fn build_multi_counts_from_membership(
    samples: &[SampleManifestRecord],
    membership_path: &Path,
    min_tu_count: Option<u64>,
) -> anyhow::Result<MultiSampleCounts> {
    let tu_ids = collect_membership_tu_ids(membership_path)?;
    let mut counter = MultiSampleCounter::new(samples, tu_ids);
    visit_membership_assignments(membership_path, |read_id, tu_id| {
        counter.add_assignment(read_id, tu_id)?;
        Ok(())
    })?;
    Ok(filter_multi_counts(counter.finish(), min_tu_count))
}

fn run_cluster_mode(
    cli: &Cli,
    threads: usize,
    samples: Option<&[SampleManifestRecord]>,
) -> anyhow::Result<()> {
    let total_start = Instant::now();

    let out_tu = cli
        .out_tu
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("--out-tu is required when clustering"))?;
    let out_membership = cli
        .out_membership
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("--out-membership is required when clustering"))?;

    if (cli.out_tu_gene.is_some() || cli.out_tu_bed12.is_some()) && cli.annotation_bed.is_none() {
        anyhow::bail!("--annotation-bed is required when using --out-tu-gene/--out-tu-bed12");
    }

    let t_read_start = Instant::now();
    let mut reads = match samples {
        Some(samples) => {
            let manifest_format = resolve_manifest_input_format(cli, samples)?;
            match manifest_format {
                ResolvedInputFormat::Bed6
                | ResolvedInputFormat::Bed12
                | ResolvedInputFormat::Tsv => {
                    load_manifest_reads(samples, manifest_format)?
                }
                ResolvedInputFormat::Fastq => anyhow::bail!(
                    "manifest `reads` entries cannot be FASTQ here; use `trackclustertu map` or `trackclustertu run`, or preprocess each FASTQ to BED6 before calling `trackclustertu cluster`"
                ),
            }
        }
        None => {
            let input_format = resolve_single_input_format(cli)?;
            read_reads(
                cli.input
                    .as_ref()
                    .expect("single-input mode must have --in"),
                input_format,
            )?
        }
    };
    filter_reads(&mut reads, cli.min_read_len);
    let t_read = t_read_start.elapsed();

    let t_pooled_reads = if let Some(out_pooled_reads) = cli.out_pooled_reads.as_ref() {
        let start = Instant::now();
        write_pooled_reads_bed6(out_pooled_reads, &reads)?;
        Some(start.elapsed())
    } else {
        None
    };

    if reads.is_empty() {
        let t_empty_start = Instant::now();
        write_empty_cluster_outputs(cli, samples)?;
        if cli.timings {
            eprintln!("[timings] threads={threads}");
            eprintln!("[timings] read_input={t_read:?} (reads=0)");
            if let Some(t_pooled_reads) = t_pooled_reads {
                eprintln!("[timings] write_pooled_reads={t_pooled_reads:?}");
            }
            eprintln!(
                "[timings] write_empty_outputs={:?}",
                t_empty_start.elapsed()
            );
            eprintln!("[timings] total={:?}", total_start.elapsed());
        }
        return Ok(());
    }

    let t_cluster_start = Instant::now();
    let clustering_options = TuClusteringOptions {
        attach_contained_reads: !cli.skip_score2_attachment,
        three_prime_tolerance_bp: cli.three_prime_tolerance_bp,
        max_five_prime_delta_bp: cli.max_five_prime_delta_bp,
    };
    let (result, stats) = if cli.timings {
        let (result, stats) = crate::tu::cluster_tus_with_stats_options(
            &reads,
            cli.score1_threshold,
            cli.score2_threshold,
            clustering_options,
        )?;
        (result, Some(stats))
    } else {
        (
            crate::tu::cluster_tus_with_options(
                &reads,
                cli.score1_threshold,
                cli.score2_threshold,
                clustering_options,
            )?,
            None,
        )
    };
    let t_cluster = t_cluster_start.elapsed();

    let prepared = prepare_outputs(&result, cli.min_tu_count);

    let t_write_tu_start = Instant::now();
    let tu_records: Vec<crate::io::bed::Bed6Record> = prepared
        .tus
        .iter()
        .map(|tu| crate::io::bed::Bed6Record {
            chrom: tu.contig.clone(),
            start: tu.interval.start,
            end: tu.interval.end,
            name: tu.id.clone(),
            score: 0,
            strand: tu.strand,
            extra_fields: Vec::new(),
        })
        .collect();
    crate::io::bed::write_bed6(out_tu, tu_records.iter())?;
    let t_write_tu = t_write_tu_start.elapsed();

    let mut read_indices: Vec<usize> = (0..reads.len()).collect();
    read_indices
        .par_sort_unstable_by(|&a, &b| cmp_read(&reads[a], &reads[b]).then_with(|| a.cmp(&b)));

    let want_tu_members = cli.out_tu_bed12.is_some();
    let mut tu_members: Option<Vec<Vec<usize>>> = if want_tu_members {
        Some(vec![Vec::new(); prepared.tus.len()])
    } else {
        None
    };

    let mut multi_counter = samples.and_then(|samples| {
        if want_any_multi_sample_outputs(cli) {
            Some(MultiSampleCounter::new(
                samples,
                prepared.tus.iter().map(|tu| tu.id.clone()).collect(),
            ))
        } else {
            None
        }
    });

    let t_membership_start = Instant::now();
    let mut writer = std::io::BufWriter::new(File::create(out_membership)?);
    for read_idx in read_indices {
        let read = &reads[read_idx];
        let old_tu_idx = result.read_to_tu[read_idx];
        let Some(tu_idx) = prepared.old_to_new[old_tu_idx] else {
            continue;
        };
        let tu = &prepared.tus[tu_idx];

        if let Some(members) = tu_members.as_mut() {
            members[tu_idx].push(read_idx);
        }
        if let Some(counter) = multi_counter.as_mut() {
            counter.add_assignment(&read.id, &tu.id)?;
        }

        let score1 = crate::score::score1_interval(read.interval, tu.interval);
        let score2 = crate::score::score2_interval(read.interval, tu.interval);
        writeln!(writer, "{}\t{}\t{score1:.6}\t{score2:.6}", read.id, tu.id)?;
    }
    let t_membership = t_membership_start.elapsed();

    if let Some(out_tu_count) = cli.out_tu_count.as_ref() {
        let t_tu_count_start = Instant::now();
        let tu_ids: Vec<String> = prepared.tus.iter().map(|tu| tu.id.clone()).collect();
        write_total_counts_csv(out_tu_count, "tu_id", &tu_ids, &prepared.counts_kept)?;
        if cli.timings {
            eprintln!("[timings] write_tu_count={:?}", t_tu_count_start.elapsed());
        }
    }

    let multi_counts = if let Some(counter) = multi_counter {
        let t_multi_count_start = Instant::now();
        let counts = counter.finish();
        if want_multi_sample_outputs(cli) {
            write_multi_sample_outputs(cli, &counts)?;
        }
        if cli.timings {
            eprintln!(
                "[timings] write_multi_sample_counts={:?}",
                t_multi_count_start.elapsed()
            );
        }
        Some(counts)
    } else {
        None
    };

    let needs_annotation = cli.annotation_bed.is_some()
        && (cli.out_tu_gene.is_some() || cli.out_tu_bed12.is_some() || want_gene_outputs(cli));
    if needs_annotation {
        let t_annotation_start = Instant::now();
        let annotation_path = cli.annotation_bed.as_ref().expect("checked");
        let genes = read_annotation_bed6(annotation_path)?;
        let overlaps = tu_gene_overlaps(&prepared.tus, &genes);

        if let Some(out_tu_gene) = cli.out_tu_gene.as_ref() {
            let mut writer = std::io::BufWriter::new(File::create(out_tu_gene)?);
            writeln!(
                writer,
                "#contig\tstrand\ttu_id\ttu_start\ttu_end\tgene_id\tgene_start\tgene_end\toverlap_bp"
            )?;

            for (tu_idx, tu) in prepared.tus.iter().enumerate() {
                for &(gene_idx, overlap_bp) in &overlaps[tu_idx] {
                    let gene = &genes[gene_idx];
                    writeln!(
                        writer,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        tu.contig,
                        tu.strand.as_char(),
                        tu.id,
                        tu.interval.start.get(),
                        tu.interval.end.get(),
                        gene.id,
                        gene.interval.start.get(),
                        gene.interval.end.get(),
                        overlap_bp
                    )?;
                }
            }
        }

        if let Some(out_tu_bed12) = cli.out_tu_bed12.as_ref() {
            let tu_members = tu_members
                .as_ref()
                .expect("tu_members must be built when --out-tu-bed12 is set");
            let mut transcripts: Vec<Transcript> = Vec::with_capacity(prepared.tus.len());
            for (tu_idx, tu) in prepared.tus.iter().enumerate() {
                let read_count = tu_members[tu_idx].len();
                let name2 = format_name2(&tu_members[tu_idx], &reads, read_count);

                let gene_names: Vec<String> = overlaps[tu_idx]
                    .iter()
                    .map(|(gene_idx, _)| genes[*gene_idx].id.clone())
                    .collect();

                let exons: Vec<Interval> = if gene_names.is_empty() {
                    vec![tu.interval]
                } else {
                    overlaps[tu_idx]
                        .iter()
                        .filter_map(|(gene_idx, _)| {
                            tu.interval.intersection(genes[*gene_idx].interval)
                        })
                        .collect()
                };

                let gene_field = if gene_names.is_empty() {
                    ".".to_owned()
                } else {
                    gene_names.join(",")
                };

                let transcript = Transcript::new(
                    tu.contig.clone(),
                    tu.strand,
                    tu.interval.start,
                    tu.interval.end,
                    tu.id.clone(),
                    exons,
                    Bed12Attrs {
                        score: 0,
                        thick_start: tu.interval.start,
                        thick_end: tu.interval.end,
                        item_rgb: "0".to_owned(),
                        extra_fields: vec![name2, gene_field],
                    },
                )?;
                transcripts.push(transcript);
            }

            crate::io::bed::write_bed12(out_tu_bed12, transcripts.iter())?;
        }

        if want_gene_outputs(cli) {
            write_gene_outputs(
                cli,
                &genes,
                &overlaps,
                &prepared.counts_kept,
                multi_counts.as_ref(),
            )?;
        }

        if cli.timings {
            eprintln!(
                "[timings] annotation_and_outputs={:?}",
                t_annotation_start.elapsed()
            );
        }
    }

    if cli.timings {
        eprintln!("[timings] threads={threads}");
        eprintln!(
            "[timings] read_input={t_read:?} (reads={}, tus={})",
            reads.len(),
            prepared.tus.len()
        );
        if let Some(t_pooled_reads) = t_pooled_reads {
            eprintln!("[timings] write_pooled_reads={t_pooled_reads:?}");
        }
        if let Some(stats) = stats.as_ref() {
            eprintln!(
                "[timings] partitions={} max_partition_reads={}",
                stats.partition_count, stats.max_partition_reads
            );
            eprintln!(
                "[timings] regions={} max_region_reads={}",
                stats.region_count, stats.max_region_reads
            );
        }
        eprintln!("[timings] cluster_tus={t_cluster:?}");
        eprintln!("[timings] write_tu_bed={t_write_tu:?}");
        eprintln!("[timings] write_membership={t_membership:?}");
        eprintln!("[timings] total={:?}", total_start.elapsed());
    }

    Ok(())
}

fn run_counts_only_mode(
    cli: &Cli,
    threads: usize,
    samples: &[SampleManifestRecord],
    membership: &Path,
) -> anyhow::Result<()> {
    if cli.out_tu.is_some()
        || cli.out_membership.is_some()
        || cli.out_pooled_reads.is_some()
        || cli.annotation_bed.is_some()
        || cli.out_tu_gene.is_some()
        || cli.out_tu_bed12.is_some()
        || cli.out_gene_count.is_some()
        || cli.out_gene_sample_count_matrix.is_some()
        || cli.out_gene_group_count_matrix.is_some()
    {
        anyhow::bail!(
            "--pooled-membership mode only writes TU count tables; do not combine it with TU, annotation, or gene outputs"
        );
    }
    if !want_any_tu_count_outputs(cli) {
        anyhow::bail!(
            "--pooled-membership mode requires at least one of --out-tu-count, --out-tu-sample-count-long, --out-tu-sample-count-matrix, or --out-tu-group-count-matrix"
        );
    }

    let total_start = Instant::now();
    let t_count_start = Instant::now();
    let counts = build_multi_counts_from_membership(samples, membership, cli.min_tu_count)?;
    let t_count = t_count_start.elapsed();

    if let Some(out_tu_count) = cli.out_tu_count.as_ref() {
        let totals = total_counts_from_multi(&counts);
        write_total_counts_csv(out_tu_count, "tu_id", &counts.tu_ids, &totals)?;
    }
    write_multi_sample_outputs(cli, &counts)?;

    if cli.timings {
        eprintln!("[timings] threads={threads}");
        eprintln!(
            "[timings] count_from_membership={t_count:?} (tus={})",
            counts.tu_ids.len()
        );
        eprintln!("[timings] total={:?}", total_start.elapsed());
    }

    Ok(())
}

fn run_cli(mut cli: Cli) -> anyhow::Result<()> {
    let mode = resolve_mode(&cli)?;
    apply_default_outputs(&mut cli, &mode);
    prepare_output_dirs(&cli)?;

    let threads = cli.threads.unwrap_or_else(default_threads);
    if threads == 0 {
        anyhow::bail!("--threads must be >= 1");
    }

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| anyhow::anyhow!("failed to build thread pool (threads={threads}): {e}"))?;
    pool.install(|| match &mode {
        InputMode::Single => run_cluster_mode(&cli, threads, None),
        InputMode::ManifestCluster { samples } => {
            run_cluster_mode(&cli, threads, Some(samples.as_slice()))
        }
        InputMode::ManifestCountOnly {
            samples,
            membership,
        } => run_counts_only_mode(&cli, threads, samples.as_slice(), membership.as_path()),
    })
}

pub(crate) fn run_cluster_from_args<I, T>(args: I) -> anyhow::Result<()>
where
    I: IntoIterator<Item = T>,
    T: Into<std::ffi::OsString> + Clone,
{
    let cli = ClusterCli::parse_from(args);
    run_cli(cli.into())
}

pub(crate) fn run_recount_from_args<I, T>(args: I) -> anyhow::Result<()>
where
    I: IntoIterator<Item = T>,
    T: Into<std::ffi::OsString> + Clone,
{
    let cli = RecountCli::parse_from(args);
    run_cli(cli.into())
}
