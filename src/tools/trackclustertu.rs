use std::ffi::OsString;
use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use anyhow::{Context, Result};
use clap::Parser;
use serde_json::{json, Value};

use crate::io::delimited::{DelimitedWriter, Delimiter};
use crate::tools::cluster_pipeline::{
    AnnotationRequest, ClusterConfig, ClusterInput, ClusterOptions, ClusterOutputPaths,
    ClusterRequest, InputFormat, TuIdStyle,
};
use crate::tools::output_transaction::OutputTransaction;
use crate::tools::run_manifest::{
    Executable, InputDescriptor, ManifestOsValue, RunManifest, RunManifestPublisher,
};
use crate::tu::multi::{parse_manifest, SampleManifestRecord};

const TOP_LEVEL_HELP: &str = "\
trackclustertu

Usage:
  trackclustertu <command> [options]

Commands:
  run         Full pipeline from FASTQ manifest to TU/gene counts
  map         FASTQ manifest to sorted BAMs plus BED manifest
  bam-to-bed  BAM input(s) to BED6 plus optional evidence sidecars
  cluster     BED input(s) to TU/gene outputs
  recount     Membership TSV to count tables
  diagnose-missed-tus  Report high-support boundary modes missing from current TU calls
  rescue-missed-tus    Promote high-support boundary modes into rescued TU calls
  gff-to-bed  GFF3 gene annotations to BED6

Use `trackclustertu <command> --help` for command-specific flags.
";

#[derive(Parser, Debug, Clone)]
#[command(
    name = "trackclustertu run",
    version,
    about = "Run the full trackclusterTU workflow from FASTQ manifest to TU/gene counts"
)]
struct RunCli {
    /// FASTQ manifest TSV with columns: sample, reads, [group].
    #[arg(long)]
    manifest: PathBuf,

    /// Reference FASTA used by minimap2.
    #[arg(long)]
    reference_fasta: PathBuf,

    /// Optional existing annotation BED6.
    #[arg(long, conflicts_with = "annotation_gff")]
    annotation_bed: Option<PathBuf>,

    /// Optional GFF3 annotation converted to BED before clustering.
    #[arg(long, conflicts_with = "annotation_bed")]
    annotation_gff: Option<PathBuf>,

    /// Output directory for mapping intermediates and clustering results.
    #[arg(long)]
    out_dir: PathBuf,

    /// Number of worker threads to use (default: all logical CPUs).
    ///
    /// Mapping totals of at least two are split between minimap2 and samtools sort.
    /// A requested total of one still assigns one thread to each of those two stages.
    #[arg(long)]
    threads: Option<usize>,

    /// Minimum MAPQ kept during BAM generation.
    #[arg(long, default_value_t = 0)]
    min_mapq: u8,

    /// One additional minimap2 argument. Repeat this flag to preserve argument boundaries.
    #[arg(
        long = "minimap2-arg",
        action = clap::ArgAction::Append,
        allow_hyphen_values = true,
        conflicts_with = "minimap2_args_compat"
    )]
    minimap2_arg: Vec<OsString>,

    /// Deprecated: whitespace-separated minimap2 arguments.
    ///
    /// Use repeated --minimap2-arg values instead. This compatibility option cannot represent
    /// an individual argument containing whitespace and conflicts with --minimap2-arg.
    #[arg(
        long = "minimap2-args",
        value_name = "ARGS",
        allow_hyphen_values = true,
        conflicts_with = "minimap2_arg"
    )]
    minimap2_args_compat: Option<String>,

    /// minimap2 executable name or explicit path.
    #[arg(long, default_value = "minimap2")]
    minimap2: PathBuf,

    /// samtools executable name or explicit path.
    #[arg(long, default_value = "samtools")]
    samtools: PathBuf,

    /// Span-Jaccard threshold (overlap / union).
    #[arg(
        long = "span-jaccard-threshold",
        visible_alias = "score1-threshold",
        default_value_t = 0.95
    )]
    score1_threshold: f64,

    /// Overlap-over-longer threshold (overlap / max_len).
    #[arg(
        long = "overlap-over-longer-threshold",
        visible_alias = "score2-threshold",
        default_value_t = 0.80
    )]
    score2_threshold: f64,

    /// Allowed strand-aware 3 prime mismatch during overlap-over-longer attachment (bp).
    #[arg(long, default_value_t = 12)]
    three_prime_tolerance_bp: u32,

    /// Optional maximum strand-aware 5 prime delta allowed during the attachment pass (bp).
    ///
    /// When set, pairs within this 5 prime cap and the 3 prime tolerance may still merge even
    /// if overlap over longer falls below its threshold.
    #[arg(long = "max-5p-delta")]
    max_five_prime_delta_bp: Option<u32>,

    /// Skip overlap-over-longer attachment and keep span-Jaccard seed clusters as final TUs.
    #[arg(
        long = "skip-overlap-over-longer-attachment",
        visible_alias = "skip-score2-attachment"
    )]
    skip_score2_attachment: bool,

    /// Optional minimum read length filter (bp).
    #[arg(long)]
    min_read_len: Option<u32>,

    /// Minimum pre-assignment clustering-family support required to emit a TU.
    #[arg(long)]
    min_tu_count: Option<u64>,

    /// TU identifier strategy used by the clustering stage.
    #[arg(long, default_value = "stable", value_parser = ["stable", "sequential"])]
    tu_id_style: String,

    /// Minimum TU/gene overlap in base pairs.
    #[arg(long, default_value_t = 1)]
    gene_min_overlap_bp: u32,

    /// Minimum fraction of a TU covered by a qualifying gene relationship.
    #[arg(long, default_value_t = 0.0)]
    gene_min_tu_fraction: f64,

    /// Minimum fraction of a gene covered by a qualifying TU relationship.
    #[arg(long, default_value_t = 0.0)]
    gene_min_gene_fraction: f64,

    /// Maximum best-vs-second assignment-score margin classified as ambiguous.
    #[arg(long, default_value_t = 0.02)]
    ambiguity_margin: f64,

    /// Split ambiguous reads equally across all candidates within the ambiguity margin.
    #[arg(long)]
    fractional_assignment: bool,

    /// Print a timing breakdown to stderr.
    #[arg(long)]
    timings: bool,
}

#[derive(Parser, Debug, Clone)]
#[command(
    name = "trackclustertu map",
    version,
    about = "Map FASTQ samples to sorted BAMs and generate BAM/BED manifests"
)]
struct MapCli {
    /// FASTQ manifest TSV with columns: sample, reads, [group].
    #[arg(long)]
    manifest: PathBuf,

    /// Reference FASTA used by minimap2.
    #[arg(long)]
    reference_fasta: PathBuf,

    /// Output directory. Writes bam/, bed/, logs/, samples.bam.tsv, and samples.bed.tsv.
    #[arg(long)]
    out_dir: PathBuf,

    /// Number of worker threads to use (default: all logical CPUs).
    ///
    /// Mapping totals of at least two are split between minimap2 and samtools sort.
    /// A requested total of one still assigns one thread to each of those two stages.
    #[arg(long)]
    threads: Option<usize>,

    /// Minimum MAPQ kept during BAM generation.
    #[arg(long, default_value_t = 0)]
    min_mapq: u8,

    /// One additional minimap2 argument. Repeat this flag to preserve argument boundaries.
    #[arg(
        long = "minimap2-arg",
        action = clap::ArgAction::Append,
        allow_hyphen_values = true,
        conflicts_with = "minimap2_args_compat"
    )]
    minimap2_arg: Vec<OsString>,

    /// Deprecated: whitespace-separated minimap2 arguments.
    ///
    /// Use repeated --minimap2-arg values instead. This compatibility option cannot represent
    /// an individual argument containing whitespace and conflicts with --minimap2-arg.
    #[arg(
        long = "minimap2-args",
        value_name = "ARGS",
        allow_hyphen_values = true,
        conflicts_with = "minimap2_arg"
    )]
    minimap2_args_compat: Option<String>,

    /// minimap2 executable name or explicit path.
    #[arg(long, default_value = "minimap2")]
    minimap2: PathBuf,

    /// samtools executable name or explicit path.
    #[arg(long, default_value = "samtools")]
    samtools: PathBuf,
}

#[derive(Parser, Debug, Clone)]
#[command(
    name = "trackclustertu bam-to-bed",
    version,
    about = "Convert BAM input(s) to BED6 with optional evidence and pre-clustering filters",
    group(
        clap::ArgGroup::new("input_mode")
            .required(true)
            .args(["input_bam", "manifest"])
    )
)]
struct BamToBedCli {
    /// Single input BAM file.
    #[arg(long = "in-bam", conflicts_with = "manifest", requires = "out_bed")]
    input_bam: Option<PathBuf>,

    /// BAM manifest TSV with sample, reads, [group], [evidence], [library_profile].
    #[arg(long, conflicts_with = "input_bam", requires = "out_dir")]
    manifest: Option<PathBuf>,

    /// Output BED6 path for single-BAM mode.
    #[arg(long = "out-bed", requires = "input_bam")]
    out_bed: Option<PathBuf>,

    /// Output directory for manifest mode.
    #[arg(long, requires = "manifest")]
    out_dir: Option<PathBuf>,

    /// Optional versioned per-record evidence TSV for single-BAM mode.
    #[arg(
        long = "out-evidence",
        requires = "input_bam",
        conflicts_with = "manifest"
    )]
    out_evidence: Option<PathBuf>,

    /// Optional read evidence TSV keyed by read_name for single-BAM mode.
    ///
    /// Recognized optional columns are poly_a_evidence, five_prime_adapter_evidence,
    /// three_prime_adapter_evidence, full_length_evidence, and library_preparation.
    #[arg(
        long = "in-evidence",
        requires = "input_bam",
        conflicts_with = "manifest"
    )]
    input_evidence: Option<PathBuf>,

    /// Emit one versioned evidence TSV per sample in manifest mode.
    #[arg(long, requires = "manifest", conflicts_with = "input_bam")]
    emit_evidence: bool,

    /// Library chemistry used to interpret optional full-length evidence.
    ///
    /// direct-rna infers full length from poly(A)+5' adapter evidence; direct-cdna and
    /// pcr-cdna infer it from 5'+3' adapter evidence. The selected profile and inferred
    /// effective_full_length state are recorded in any emitted evidence regardless of
    /// --require-full-length. Only that flag uses the state to filter BED6 records.
    #[arg(long, default_value = "direct-rna")]
    library_profile: crate::bam::LibraryProfile,

    /// Minimum mapping quality retained in BED6 (default preserves all MAPQ values).
    #[arg(long, default_value_t = 0)]
    min_mapq: u8,

    /// Retain only reads with explicit or profile-inferred full-length evidence.
    #[arg(long)]
    require_full_length: bool,

    /// Minimum exact chrom/start/end/strand support retained before clustering.
    #[arg(long, default_value_t = 1)]
    min_boundary_support: usize,
}

#[derive(Parser, Debug, Clone)]
#[command(
    name = "trackclustertu gff-to-bed",
    version,
    about = "Convert GFF3 gene annotations to BED6"
)]
struct GffToBedCli {
    /// Input GFF3 annotation file.
    #[arg(long = "annotation-gff", alias = "gff")]
    annotation_gff: PathBuf,

    /// Output BED6 path.
    #[arg(long = "out-bed")]
    out_bed: PathBuf,
}

fn default_threads() -> usize {
    std::thread::available_parallelism()
        .map(usize::from)
        .unwrap_or(1)
}

fn prepend_program<I>(program: &str, args: I) -> Vec<OsString>
where
    I: IntoIterator<Item = OsString>,
{
    let mut values = vec![OsString::from(program)];
    values.extend(args);
    values
}

fn print_top_level_help() {
    println!("{TOP_LEVEL_HELP}");
}

fn ensure_path_exists(path: &Path, label: &str) -> Result<()> {
    if path.exists() {
        Ok(())
    } else {
        anyhow::bail!("{label} not found: {}", path.display())
    }
}

fn write_manifest(manifest_path: &Path, records: &[SampleManifestRecord]) -> Result<()> {
    let mut writer = DelimitedWriter::create(manifest_path, Delimiter::Tab, &[])?;
    writer.write_record(["sample", "group", "reads"])?;
    for record in records {
        writer.write_record([
            record.sample.clone(),
            record.group.clone().unwrap_or_default(),
            record.reads.display().to_string(),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

fn ensure_success(status: std::process::ExitStatus, step: &str, log_path: &Path) -> Result<()> {
    if status.success() {
        Ok(())
    } else {
        anyhow::bail!(
            "{step} failed with status {status}; see log {}",
            log_path.display()
        )
    }
}

#[derive(Clone, Debug)]
struct RawMapConfig {
    manifest: PathBuf,
    reference_fasta: PathBuf,
    out_dir: PathBuf,
    threads: Option<usize>,
    min_mapq: u8,
    minimap2_arg: Vec<OsString>,
    minimap2_args_compat: Option<String>,
    minimap2: PathBuf,
    samtools: PathBuf,
}

/// Validated, fully effective configuration for the complete `map` command.
///
/// Unlike `SampleMapConfig`, this owns canonical command inputs, resolved tools, effective
/// defaults, and parsed sample records. Clap values cross into this type exactly once.
#[derive(Clone, Debug)]
struct MapConfig {
    manifest: PathBuf,
    reference_fasta: PathBuf,
    out_dir: PathBuf,
    threads: usize,
    minimap2_threads: usize,
    samtools_sort_threads: usize,
    min_mapq: u8,
    minimap2_args: Vec<OsString>,
    minimap2: Executable,
    samtools: Executable,
    records: Vec<SampleManifestRecord>,
    library_profile: crate::bam::LibraryProfile,
}

impl MapConfig {
    fn validate(raw: RawMapConfig) -> Result<Self> {
        let manifest = canonical_input_file(&raw.manifest, "manifest")?;
        let reference_fasta = canonical_input_file(&raw.reference_fasta, "reference FASTA")?;
        let threads = raw.threads.unwrap_or_else(default_threads);
        if threads == 0 {
            anyhow::bail!("--threads must be >= 1");
        }
        let (minimap2_threads, samtools_sort_threads) = allocate_pipeline_threads(threads);
        let minimap2_args =
            effective_minimap2_args(raw.minimap2_arg, raw.minimap2_args_compat.as_deref())?;
        if raw.minimap2_args_compat.is_some() {
            eprintln!(
                "warning: --minimap2-args is deprecated; repeat --minimap2-arg for each argument"
            );
        }

        let mut records = parse_manifest(&manifest)
            .with_context(|| format!("failed to parse FASTQ manifest {manifest:?}"))?;
        crate::tools::bam_to_bed::validate_unique_sanitized_sample_names(&records)?;
        for record in &mut records {
            record.reads = canonical_input_file(
                &record.reads,
                &format!("FASTQ for sample {:?}", record.sample),
            )?;
        }

        let out_dir = absolute_path(&raw.out_dir)?;
        let minimap2 = Executable::resolve(raw.minimap2, "minimap2")?;
        let samtools = Executable::resolve(raw.samtools, "samtools")?;

        Ok(Self {
            manifest,
            reference_fasta,
            out_dir,
            threads,
            minimap2_threads,
            samtools_sort_threads,
            min_mapq: raw.min_mapq,
            minimap2_args,
            minimap2,
            samtools,
            records,
            library_profile: crate::bam::LibraryProfile::default(),
        })
    }

    fn input_descriptors(&self) -> Vec<InputDescriptor> {
        let mut inputs = vec![
            InputDescriptor::new("fastq_manifest", &self.manifest),
            InputDescriptor::new("reference_fasta", &self.reference_fasta),
        ];
        inputs.extend(
            self.records
                .iter()
                .map(|record| InputDescriptor::for_sample("fastq", &record.sample, &record.reads)),
        );
        inputs
    }

    fn effective_configuration(&self) -> Value {
        json!({
            "manifest": ManifestOsValue::path(&self.manifest),
            "reference_fasta": ManifestOsValue::path(&self.reference_fasta),
            "out_dir": ManifestOsValue::path(&self.out_dir),
            "threads": self.threads,
            "pipeline_threads": {
                "minimap2": self.minimap2_threads,
                "samtools_sort": self.samtools_sort_threads,
            },
            "min_mapq": self.min_mapq,
            "minimap2_args": self.minimap2_args.iter()
                .map(|argument| ManifestOsValue::new(argument))
                .collect::<Vec<_>>(),
            "minimap2": {
                "requested": ManifestOsValue::path(self.minimap2.requested()),
                "resolved": ManifestOsValue::path(self.minimap2.resolved()),
            },
            "samtools": {
                "requested": ManifestOsValue::path(self.samtools.requested()),
                "resolved": ManifestOsValue::path(self.samtools.resolved()),
            },
            "library_profile": self.library_profile.to_string(),
            "sample_count": self.records.len(),
        })
    }
}

impl TryFrom<MapCli> for MapConfig {
    type Error = anyhow::Error;

    fn try_from(cli: MapCli) -> Result<Self> {
        Self::validate(RawMapConfig {
            manifest: cli.manifest,
            reference_fasta: cli.reference_fasta,
            out_dir: cli.out_dir,
            threads: cli.threads,
            min_mapq: cli.min_mapq,
            minimap2_arg: cli.minimap2_arg,
            minimap2_args_compat: cli.minimap2_args_compat,
            minimap2: cli.minimap2,
            samtools: cli.samtools,
        })
    }
}

/// Validated configuration for `run`, including its owned mapping configuration.
#[derive(Clone, Debug)]
struct RunConfig {
    map: MapConfig,
    annotation_bed: Option<PathBuf>,
    annotation_gff: Option<PathBuf>,
    score1_threshold: f64,
    score2_threshold: f64,
    three_prime_tolerance_bp: u32,
    max_five_prime_delta_bp: Option<u32>,
    skip_score2_attachment: bool,
    min_read_len: Option<u32>,
    min_tu_count: Option<u64>,
    tu_id_style: String,
    gene_min_overlap_bp: u32,
    gene_min_tu_fraction: f64,
    gene_min_gene_fraction: f64,
    ambiguity_margin: f64,
    fractional_assignment: bool,
    timings: bool,
}

impl TryFrom<RunCli> for RunConfig {
    type Error = anyhow::Error;

    fn try_from(cli: RunCli) -> Result<Self> {
        let RunCli {
            manifest,
            reference_fasta,
            annotation_bed,
            annotation_gff,
            out_dir,
            threads,
            min_mapq,
            minimap2_arg,
            minimap2_args_compat,
            minimap2,
            samtools,
            score1_threshold,
            score2_threshold,
            three_prime_tolerance_bp,
            max_five_prime_delta_bp,
            skip_score2_attachment,
            min_read_len,
            min_tu_count,
            tu_id_style,
            gene_min_overlap_bp,
            gene_min_tu_fraction,
            gene_min_gene_fraction,
            ambiguity_margin,
            fractional_assignment,
            timings,
        } = cli;

        validate_unit_interval(score1_threshold, "--span-jaccard-threshold")?;
        validate_unit_interval(score2_threshold, "--overlap-over-longer-threshold")?;
        validate_unit_interval(ambiguity_margin, "--ambiguity-margin")?;
        validate_unit_interval(gene_min_tu_fraction, "--gene-min-tu-fraction")?;
        validate_unit_interval(gene_min_gene_fraction, "--gene-min-gene-fraction")?;
        let annotation_bed = annotation_bed
            .map(|path| canonical_input_file(&path, "annotation BED"))
            .transpose()?;
        let annotation_gff = annotation_gff
            .map(|path| canonical_input_file(&path, "annotation GFF"))
            .transpose()?;
        if annotation_bed.is_some() && annotation_gff.is_some() {
            anyhow::bail!("--annotation-bed conflicts with --annotation-gff");
        }

        let map = MapConfig::validate(RawMapConfig {
            manifest,
            reference_fasta,
            out_dir,
            threads,
            min_mapq,
            minimap2_arg,
            minimap2_args_compat,
            minimap2,
            samtools,
        })?;

        Ok(Self {
            map,
            annotation_bed,
            annotation_gff,
            score1_threshold,
            score2_threshold,
            three_prime_tolerance_bp,
            max_five_prime_delta_bp,
            skip_score2_attachment,
            min_read_len,
            min_tu_count,
            tu_id_style,
            gene_min_overlap_bp,
            gene_min_tu_fraction,
            gene_min_gene_fraction,
            ambiguity_margin,
            fractional_assignment,
            timings,
        })
    }
}

impl RunConfig {
    fn input_descriptors(&self) -> Vec<InputDescriptor> {
        let mut inputs = self.map.input_descriptors();
        if let Some(path) = &self.annotation_bed {
            inputs.push(InputDescriptor::new("annotation_bed", path));
        }
        if let Some(path) = &self.annotation_gff {
            inputs.push(InputDescriptor::new("annotation_gff", path));
        }
        inputs
    }

    fn effective_configuration(&self) -> Value {
        let generated_annotation_bed = self
            .annotation_gff
            .as_ref()
            .map(|_| ManifestOsValue::path(&self.map.out_dir.join("annotation.bed")));
        json!({
            "mapping": self.map.effective_configuration(),
            "annotation": {
                "bed": self.annotation_bed.as_deref().map(ManifestOsValue::path),
                "gff": self.annotation_gff.as_deref().map(ManifestOsValue::path),
                "generated_bed": generated_annotation_bed,
            },
            "clustering": {
                "input_format": "bed6",
                "threads": self.map.threads,
                "span_jaccard_threshold": self.score1_threshold,
                "overlap_over_longer_threshold": self.score2_threshold,
                "three_prime_tolerance_bp": self.three_prime_tolerance_bp,
                "max_five_prime_delta_bp": self.max_five_prime_delta_bp,
                "skip_overlap_over_longer_attachment": self.skip_score2_attachment,
                "min_read_len": self.min_read_len,
                "min_tu_count": self.min_tu_count,
                "tu_id_style": self.tu_id_style,
                "gene_min_overlap_bp": self.gene_min_overlap_bp,
                "gene_min_tu_fraction": self.gene_min_tu_fraction,
                "gene_min_gene_fraction": self.gene_min_gene_fraction,
                "ambiguity_margin": self.ambiguity_margin,
                "fractional_assignment": self.fractional_assignment,
                "timings": self.timings,
                "out_dir": ManifestOsValue::path(&self.map.out_dir),
            },
            "outputs": {
                "bam_dir": ManifestOsValue::path(&self.map.out_dir.join("bam")),
                "bed_dir": ManifestOsValue::path(&self.map.out_dir.join("bed")),
                "logs_dir": ManifestOsValue::path(&self.map.out_dir.join("logs")),
                "bam_manifest": ManifestOsValue::path(&self.map.out_dir.join("samples.bam.tsv")),
                "bed_manifest": ManifestOsValue::path(&self.map.out_dir.join("samples.bed.tsv")),
                "run_manifest": ManifestOsValue::path(
                    &self.map.out_dir.join(crate::tools::run_manifest::RUN_MANIFEST_FILE_NAME)
                ),
            },
        })
    }
}

fn validate_unit_interval(value: f64, flag: &str) -> Result<()> {
    if value.is_finite() && (0.0..=1.0).contains(&value) {
        Ok(())
    } else {
        anyhow::bail!("{flag} must be finite and between 0 and 1, got {value}")
    }
}

fn canonical_input_file(path: &Path, label: &str) -> Result<PathBuf> {
    ensure_path_exists(path, label)?;
    let canonical = path
        .canonicalize()
        .with_context(|| format!("failed to canonicalize {label} {}", path.display()))?;
    if !canonical.is_file() {
        anyhow::bail!("{label} is not a regular file: {}", canonical.display());
    }
    Ok(canonical)
}

fn absolute_path(path: &Path) -> Result<PathBuf> {
    if path.is_absolute() {
        Ok(path.to_path_buf())
    } else {
        Ok(std::env::current_dir()
            .context("failed to resolve current directory")?
            .join(path))
    }
}

fn effective_minimap2_args(
    repeated: Vec<OsString>,
    compatibility: Option<&str>,
) -> Result<Vec<OsString>> {
    if !repeated.is_empty() && compatibility.is_some() {
        anyhow::bail!("--minimap2-arg conflicts with deprecated --minimap2-args");
    }
    let mut effective = vec![OsString::from("-ax"), OsString::from("map-ont")];
    if !repeated.is_empty() {
        effective.extend(repeated);
        return Ok(effective);
    }
    if let Some(raw) = compatibility {
        effective.extend(raw.split_whitespace().map(OsString::from));
    }
    Ok(effective)
}

struct SampleMapConfig<'a> {
    reference_fasta: &'a Path,
    logs_dir: &'a Path,
    sample_stem: &'a str,
    minimap2_threads: usize,
    samtools_sort_threads: usize,
    min_mapq: u8,
    minimap2_args: &'a [OsString],
    minimap2: &'a Path,
    samtools: &'a Path,
}

fn allocate_pipeline_threads(total: usize) -> (usize, usize) {
    // The mapper and sorter each need one allocated thread; samtools view stays single-threaded.
    // With a user budget of one, both configurable tools receive their minimum. For all larger
    // budgets the allocations sum to the requested total. Give an odd spare worker to the mapper.
    if total <= 1 {
        return (1, 1);
    }
    let samtools_sort_threads = total / 2;
    let minimap2_threads = total - samtools_sort_threads;
    (minimap2_threads, samtools_sort_threads)
}

fn samtools_sort_additional_threads(allocated_threads: usize) -> usize {
    // `samtools sort -@` excludes its main thread from the requested count.
    allocated_threads.saturating_sub(1)
}

fn map_fastq_to_bam(
    config: &SampleMapConfig<'_>,
    fastq_path: &Path,
    bam_path: &Path,
) -> Result<()> {
    let minimap2_log = config
        .logs_dir
        .join(format!("{}.minimap2.log", config.sample_stem));
    let samtools_view_log = config
        .logs_dir
        .join(format!("{}.samtools_view.log", config.sample_stem));
    let samtools_sort_log = config
        .logs_dir
        .join(format!("{}.samtools_sort.log", config.sample_stem));
    let samtools_index_log = config
        .logs_dir
        .join(format!("{}.samtools_index.log", config.sample_stem));

    let mut minimap2 = Command::new(config.minimap2);
    minimap2.args(config.minimap2_args);
    minimap2
        .arg("-t")
        .arg(config.minimap2_threads.to_string())
        .arg(config.reference_fasta)
        .arg(fastq_path)
        .stdout(Stdio::piped())
        .stderr(Stdio::from(File::create(&minimap2_log).with_context(
            || format!("failed to create log {:?}", minimap2_log),
        )?));
    let mut minimap2_child = minimap2.spawn().with_context(|| {
        format!(
            "failed to start minimap2 for FASTQ {}",
            fastq_path.display()
        )
    })?;
    let minimap2_stdout = minimap2_child
        .stdout
        .take()
        .context("failed to capture minimap2 stdout")?;

    let mut samtools_view = Command::new(config.samtools);
    samtools_view
        .arg("view")
        .arg("-b")
        .arg("-F")
        .arg("260")
        .arg("-F")
        .arg("2048")
        .arg("-q")
        .arg(config.min_mapq.to_string())
        .arg("-")
        .stdin(Stdio::from(minimap2_stdout))
        .stdout(Stdio::piped())
        .stderr(Stdio::from(File::create(&samtools_view_log).with_context(
            || format!("failed to create log {:?}", samtools_view_log),
        )?));
    let mut samtools_view_child = samtools_view.spawn().with_context(|| {
        format!(
            "failed to start samtools view for FASTQ {}",
            fastq_path.display()
        )
    })?;
    let samtools_view_stdout = samtools_view_child
        .stdout
        .take()
        .context("failed to capture samtools view stdout")?;

    let mut samtools_sort = Command::new(config.samtools);
    samtools_sort
        .arg("sort")
        .arg("-@")
        .arg(samtools_sort_additional_threads(config.samtools_sort_threads).to_string())
        .arg("-o")
        .arg(bam_path)
        .arg("-")
        .stdin(Stdio::from(samtools_view_stdout))
        .stderr(Stdio::from(File::create(&samtools_sort_log).with_context(
            || format!("failed to create log {:?}", samtools_sort_log),
        )?));
    let mut samtools_sort_child = samtools_sort.spawn().with_context(|| {
        format!(
            "failed to start samtools sort for FASTQ {}",
            fastq_path.display()
        )
    })?;

    let sort_status = samtools_sort_child.wait()?;
    let view_status = samtools_view_child.wait()?;
    let minimap2_status = minimap2_child.wait()?;

    ensure_success(sort_status, "samtools sort", &samtools_sort_log)?;
    ensure_success(view_status, "samtools view", &samtools_view_log)?;
    ensure_success(minimap2_status, "minimap2", &minimap2_log)?;

    let index_status = Command::new(config.samtools)
        .arg("index")
        .arg(bam_path)
        .stderr(Stdio::from(
            File::create(&samtools_index_log)
                .with_context(|| format!("failed to create log {:?}", samtools_index_log))?,
        ))
        .status()
        .with_context(|| format!("failed to start samtools index for {}", bam_path.display()))?;
    ensure_success(index_status, "samtools index", &samtools_index_log)?;

    Ok(())
}

fn execute_map(config: &MapConfig) -> Result<PathBuf> {
    fs::create_dir_all(&config.out_dir)
        .with_context(|| format!("failed to create output directory {:?}", config.out_dir))?;
    let bam_dir = config.out_dir.join("bam");
    let logs_dir = config.out_dir.join("logs");
    fs::create_dir_all(&bam_dir)
        .with_context(|| format!("failed to create BAM output directory {:?}", bam_dir))?;
    fs::create_dir_all(&logs_dir)
        .with_context(|| format!("failed to create log output directory {:?}", logs_dir))?;

    let mut bam_records: Vec<SampleManifestRecord> = Vec::with_capacity(config.records.len());
    for record in &config.records {
        let sample_stem = crate::tools::bam_to_bed::sanitize_sample_name(&record.sample);
        let bam_path = bam_dir.join(format!("{sample_stem}.sorted.bam"));
        let sample_config = SampleMapConfig {
            reference_fasta: &config.reference_fasta,
            logs_dir: &logs_dir,
            sample_stem: &sample_stem,
            minimap2_threads: config.minimap2_threads,
            samtools_sort_threads: config.samtools_sort_threads,
            min_mapq: config.min_mapq,
            minimap2_args: &config.minimap2_args,
            minimap2: config.minimap2.resolved(),
            samtools: config.samtools.resolved(),
        };
        map_fastq_to_bam(&sample_config, &record.reads, &bam_path)?;
        let bam_path = bam_path
            .canonicalize()
            .with_context(|| format!("failed to canonicalize BAM path {:?}", bam_path))?;

        bam_records.push(SampleManifestRecord {
            sample: record.sample.clone(),
            reads: bam_path,
            group: record.group.clone(),
            evidence: None,
        });
    }

    let bam_manifest = config.out_dir.join("samples.bam.tsv");
    write_manifest(&bam_manifest, &bam_records)?;
    let bam_manifest = bam_manifest
        .canonicalize()
        .with_context(|| format!("failed to canonicalize BAM manifest {:?}", bam_manifest))?;

    let bed_manifest =
        crate::tools::bam_to_bed::convert_records_to_bed_manifest(&bam_records, &config.out_dir)?;
    let bed_manifest = bed_manifest
        .canonicalize()
        .with_context(|| format!("failed to canonicalize BED manifest {:?}", bed_manifest))?;

    println!("bam_manifest={}", bam_manifest.display());
    println!("bed_manifest={}", bed_manifest.display());

    Ok(bed_manifest)
}

fn run_map(config: &MapConfig) -> Result<PathBuf> {
    let inputs = config.input_descriptors();
    let publisher = RunManifestPublisher::new(&config.out_dir, &inputs)?;
    let timestamp = crate::tools::run_manifest::timestamp_now()?;
    let manifest = RunManifest::capture(
        "map",
        config.effective_configuration(),
        config.library_profile.to_string(),
        &inputs,
        &[
            ("minimap2", &config.minimap2),
            ("samtools", &config.samtools),
        ],
        timestamp,
    )?;

    let bed_manifest = execute_map(config)?;
    let manifest_path = publisher.final_path().to_path_buf();
    publisher.publish(&manifest)?;
    println!("run_manifest={}", manifest_path.display());
    Ok(bed_manifest)
}

fn run_map_from_args<I>(args: I) -> Result<()>
where
    I: IntoIterator<Item = OsString>,
{
    let cli = MapCli::parse_from(args);
    let config = MapConfig::try_from(cli)?;
    run_map(&config).map(|_| ())
}

fn run_bam_to_bed_from_args<I>(args: I) -> Result<()>
where
    I: IntoIterator<Item = OsString>,
{
    let cli = BamToBedCli::parse_from(args);
    let config = crate::bam::BamConversionConfig {
        min_mapq: cli.min_mapq,
        require_full_length: cli.require_full_length,
        min_boundary_support: cli.min_boundary_support,
        library_profile: cli.library_profile,
    };
    match (cli.input_bam.as_deref(), cli.manifest.as_deref()) {
        (Some(input_bam), None) => {
            let out_bed = cli.out_bed.as_deref().expect("required by clap");
            let mut inputs = vec![input_bam.to_path_buf()];
            if let Some(path) = &cli.input_evidence {
                inputs.push(path.clone());
            }
            let mut outputs = vec![out_bed.to_path_buf()];
            if let Some(path) = &cli.out_evidence {
                outputs.push(path.clone());
            }
            let transaction = OutputTransaction::new(inputs, outputs)?;
            let staged_bed = transaction.staged_path(out_bed)?;
            let staged_evidence = cli
                .out_evidence
                .as_deref()
                .map(|path| transaction.staged_path(path))
                .transpose()?;
            let summary = crate::tools::bam_to_bed::convert_single_bam_to_bed_with_config(
                input_bam,
                &staged_bed,
                staged_evidence.as_deref(),
                cli.input_evidence.as_deref(),
                &config,
            )?;
            transaction.commit()?;
            println!("out_bed={}", out_bed.display());
            if let Some(path) = &cli.out_evidence {
                println!("out_evidence={}", path.display());
            }
            print_bam_conversion_summary(&summary);
            Ok(())
        }
        (None, Some(manifest)) => {
            let out_dir = cli.out_dir.as_deref().expect("required by clap");
            let result = crate::tools::bam_to_bed::convert_manifest_to_bed_manifest_with_config(
                manifest,
                out_dir,
                &config,
                cli.emit_evidence,
            )?;
            println!("bed_manifest={}", result.bed_manifest.display());
            print_bam_conversion_summary(&result.summary);
            Ok(())
        }
        _ => anyhow::bail!("pass either --in-bam/--out-bed or --manifest/--out-dir"),
    }
}

fn print_bam_conversion_summary(summary: &crate::bam::BamConversionSummary) {
    eprint!("bam_to_bed_counts");
    for (name, value) in summary.fields() {
        eprint!("\t{name}={value}");
    }
    eprintln!();
}

fn run_gff_to_bed_from_args<I>(args: I) -> Result<()>
where
    I: IntoIterator<Item = OsString>,
{
    let cli = GffToBedCli::parse_from(args);
    let transaction = OutputTransaction::new([&cli.annotation_gff], [&cli.out_bed])?;
    let staged_bed = transaction.staged_path(&cli.out_bed)?;
    let gene_count =
        crate::tools::gff_to_bed::convert_gff_to_bed(&cli.annotation_gff, &staged_bed)?;
    transaction.commit()?;
    println!("genes={gene_count}");
    println!("out_bed={}", cli.out_bed.display());
    Ok(())
}

fn run_full_pipeline(config: &RunConfig) -> Result<()> {
    let inputs = config.input_descriptors();
    let publisher = RunManifestPublisher::new(&config.map.out_dir, &inputs)?;
    let timestamp = crate::tools::run_manifest::timestamp_now()?;
    let manifest = RunManifest::capture(
        "run",
        config.effective_configuration(),
        config.map.library_profile.to_string(),
        &inputs,
        &[
            ("minimap2", &config.map.minimap2),
            ("samtools", &config.map.samtools),
        ],
        timestamp,
    )?;

    // `run` owns publication of the final run manifest. The mapping stage therefore executes
    // directly instead of publishing an intermediate map-only manifest at the same path.
    let bed_manifest = execute_map(&config.map)?;

    let annotation_bed = match (&config.annotation_bed, &config.annotation_gff) {
        (Some(path), None) => Some(path.clone()),
        (None, Some(gff)) => {
            let out_bed = config.map.out_dir.join("annotation.bed");
            let gene_count = crate::tools::gff_to_bed::convert_gff_to_bed(gff, &out_bed)?;
            println!("genes={gene_count}");
            println!("out_bed={}", out_bed.display());
            Some(
                out_bed
                    .canonicalize()
                    .with_context(|| format!("failed to canonicalize BED path {:?}", out_bed))?,
            )
        }
        (None, None) => None,
        (Some(_), Some(_)) => unreachable!("clap enforces conflict"),
    };

    let tu_id_style = match config.tu_id_style.as_str() {
        "stable" => TuIdStyle::Stable,
        "sequential" => TuIdStyle::Sequential,
        _ => unreachable!("RunCli constrains --tu-id-style"),
    };
    let annotation = annotation_bed.map(|bed| AnnotationRequest {
        bed,
        min_overlap_bp: config.gene_min_overlap_bp,
        min_tu_fraction: config.gene_min_tu_fraction,
        min_gene_fraction: config.gene_min_gene_fraction,
    });
    let cluster_config = ClusterConfig::validate(ClusterRequest {
        input: ClusterInput::Manifest(bed_manifest),
        format: InputFormat::Bed6,
        outputs: ClusterOutputPaths {
            out_dir: Some(config.map.out_dir.clone()),
            ..ClusterOutputPaths::default()
        },
        options: ClusterOptions {
            span_jaccard_threshold: config.score1_threshold,
            overlap_over_longer_threshold: config.score2_threshold,
            three_prime_tolerance_bp: config.three_prime_tolerance_bp,
            max_five_prime_delta_bp: config.max_five_prime_delta_bp,
            attach_contained_reads: !config.skip_score2_attachment,
            ambiguity_margin: config.ambiguity_margin,
            fractional_assignment: config.fractional_assignment,
            tu_id_style,
            strict_read_errors: false,
        },
        min_read_len: config.min_read_len,
        min_tu_count: config.min_tu_count,
        annotation,
        threads: Some(config.map.threads),
        timings: config.timings,
    })?;
    crate::tools::cluster_pipeline::run_cluster(cluster_config)?;
    let manifest_path = publisher.final_path().to_path_buf();
    publisher.publish(&manifest)?;
    println!("run_manifest={}", manifest_path.display());
    Ok(())
}

fn run_full_from_args<I>(args: I) -> Result<()>
where
    I: IntoIterator<Item = OsString>,
{
    let cli = RunCli::parse_from(args);
    let config = RunConfig::try_from(cli)?;
    run_full_pipeline(&config)
}

pub fn entrypoint() -> Result<()> {
    let mut args = std::env::args_os();
    let _program = args.next();

    let Some(command) = args.next() else {
        print_top_level_help();
        return Ok(());
    };

    match command.to_str() {
        Some("-h") | Some("--help") | Some("help") => {
            print_top_level_help();
            Ok(())
        }
        Some("-V") | Some("--version") => {
            println!("trackclustertu {}", env!("CARGO_PKG_VERSION"));
            Ok(())
        }
        Some("run") => run_full_from_args(prepend_program("trackclustertu run", args)),
        Some("map") => run_map_from_args(prepend_program("trackclustertu map", args)),
        Some("bam-to-bed") => {
            run_bam_to_bed_from_args(prepend_program("trackclustertu bam-to-bed", args))
        }
        Some("cluster") => crate::tools::cluster_pipeline::run_cluster_from_args(prepend_program(
            "trackclustertu cluster",
            args,
        )),
        Some("recount") => crate::tools::cluster_pipeline::run_recount_from_args(prepend_program(
            "trackclustertu recount",
            args,
        )),
        Some("diagnose-missed-tus") => crate::tools::diagnose_missed_tus::run_from_args(
            prepend_program("trackclustertu diagnose-missed-tus", args),
        ),
        Some("rescue-missed-tus") => crate::tools::diagnose_missed_tus::run_rescue_from_args(
            prepend_program("trackclustertu rescue-missed-tus", args),
        ),
        Some("gff-to-bed") => {
            run_gff_to_bed_from_args(prepend_program("trackclustertu gff-to-bed", args))
        }
        Some(other) if other.starts_with('-') => {
            anyhow::bail!(
                "subcommands are required; use `trackclustertu cluster ...` or `trackclustertu run ...`\n\n{TOP_LEVEL_HELP}"
            )
        }
        Some(other) => anyhow::bail!("unknown subcommand {other:?}\n\n{TOP_LEVEL_HELP}"),
        None => unreachable!("command already checked"),
    }
}

#[cfg(test)]
mod tests {
    use std::ffi::OsString;
    use std::path::PathBuf;

    use clap::error::ErrorKind;
    use clap::Parser;

    use super::{
        allocate_pipeline_threads, effective_minimap2_args, samtools_sort_additional_threads,
        MapCli,
    };

    #[test]
    fn concurrent_mapper_and_sorter_split_the_thread_budget() {
        assert_eq!(allocate_pipeline_threads(1), (1, 1));
        assert_eq!(allocate_pipeline_threads(2), (1, 1));
        assert_eq!(allocate_pipeline_threads(3), (2, 1));
        assert_eq!(allocate_pipeline_threads(8), (4, 4));

        for total in 2..=32 {
            let (mapper, sorter) = allocate_pipeline_threads(total);
            assert!(mapper >= 1);
            assert!(sorter >= 1);
            assert_eq!(mapper + sorter, total);
        }

        assert_eq!(samtools_sort_additional_threads(1), 0);
        assert_eq!(samtools_sort_additional_threads(4), 3);
    }

    #[test]
    fn repeated_minimap2_arguments_preserve_spaces_and_boundaries() {
        let cli = MapCli::try_parse_from([
            OsString::from("trackclustertu map"),
            OsString::from("--manifest"),
            OsString::from("manifest with spaces.tsv"),
            OsString::from("--reference-fasta"),
            OsString::from("reference with spaces.fa"),
            OsString::from("--out-dir"),
            OsString::from("output with spaces"),
            OsString::from("--minimap2-arg"),
            OsString::from("-ax"),
            OsString::from("--minimap2-arg"),
            OsString::from("one argument with spaces"),
            OsString::from("--minimap2"),
            OsString::from("tools with spaces/minimap2"),
            OsString::from("--samtools"),
            OsString::from("tools with spaces/samtools"),
        ])
        .unwrap();

        assert_eq!(
            cli.minimap2_arg,
            vec![
                OsString::from("-ax"),
                OsString::from("one argument with spaces")
            ]
        );
        assert_eq!(cli.minimap2, PathBuf::from("tools with spaces/minimap2"));
        assert_eq!(cli.samtools, PathBuf::from("tools with spaces/samtools"));
    }

    #[test]
    fn deprecated_minimap2_args_conflicts_with_lossless_form() {
        let error = MapCli::try_parse_from([
            "trackclustertu map",
            "--manifest",
            "samples.tsv",
            "--reference-fasta",
            "reference.fa",
            "--out-dir",
            "output",
            "--minimap2-arg",
            "-ax",
            "--minimap2-args",
            "-x map-ont",
        ])
        .unwrap_err();
        assert_eq!(error.kind(), ErrorKind::ArgumentConflict);
    }

    #[test]
    fn minimap2_defaults_precede_additional_arguments() {
        assert_eq!(
            effective_minimap2_args(Vec::new(), None).unwrap(),
            vec![OsString::from("-ax"), OsString::from("map-ont")]
        );
        assert_eq!(
            effective_minimap2_args(vec![OsString::from("--secondary=no")], None).unwrap(),
            vec![
                OsString::from("-ax"),
                OsString::from("map-ont"),
                OsString::from("--secondary=no")
            ]
        );
        assert_eq!(
            effective_minimap2_args(Vec::new(), Some("-ax map-pb --secondary=no")).unwrap(),
            vec![
                OsString::from("-ax"),
                OsString::from("map-ont"),
                OsString::from("-ax"),
                OsString::from("map-pb"),
                OsString::from("--secondary=no")
            ]
        );

        let cli = MapCli::try_parse_from([
            "trackclustertu map",
            "--manifest",
            "samples.tsv",
            "--reference-fasta",
            "reference.fa",
            "--out-dir",
            "output",
            "--minimap2-args",
            "-ax map-ont",
        ])
        .unwrap();
        assert_eq!(cli.minimap2_args_compat.as_deref(), Some("-ax map-ont"));
    }

    #[cfg(unix)]
    #[test]
    fn repeated_minimap2_arguments_accept_non_utf8_os_values() {
        use std::os::unix::ffi::OsStringExt;

        let raw = OsString::from_vec(vec![b'a', 0x80, b' ', b'b']);
        let cli = MapCli::try_parse_from(vec![
            OsString::from("trackclustertu map"),
            OsString::from("--manifest"),
            OsString::from("samples.tsv"),
            OsString::from("--reference-fasta"),
            OsString::from("reference.fa"),
            OsString::from("--out-dir"),
            OsString::from("output"),
            OsString::from("--minimap2-arg"),
            raw.clone(),
        ])
        .unwrap();
        assert_eq!(cli.minimap2_arg, vec![raw]);
    }
}
