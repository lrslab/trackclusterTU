use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;

use crate::io::bed::Bed6Record;
use crate::model::{Coord, Interval, Strand};
use crate::score::{score1_interval, score2_interval};
use crate::tu::ReadRecord;

#[derive(Parser, Debug, Clone)]
#[command(
    name = "trackclustertu diagnose-missed-tus",
    version,
    about = "Report high-support boundary modes that are not represented by existing TU calls"
)]
struct DiagnoseMissedTusCli {
    /// Input BED6 reads file used for clustering.
    #[arg(long = "in")]
    input: PathBuf,

    /// Existing TU BED6 output (for example, `results/tus.bed`).
    #[arg(long = "existing-tu")]
    existing_tu: PathBuf,

    /// Optional gene annotation BED6 used to label candidates.
    #[arg(long)]
    annotation_bed: Option<PathBuf>,

    /// Output TSV report of candidate missed TUs.
    #[arg(long = "out-tsv")]
    out_tsv: PathBuf,

    /// Optional BED6 file for loading candidate missed TUs into IGV.
    #[arg(long = "out-bed")]
    out_bed: Option<PathBuf>,

    /// Strand-aware 3 prime window used to form end families and TU matches (bp).
    #[arg(long, default_value_t = 12)]
    three_prime_window_bp: u32,

    /// Strand-aware 5 prime window used to call start modes and TU matches (bp).
    #[arg(long, default_value_t = 10)]
    five_prime_window_bp: u32,

    /// Minimum reads required in a 3 prime family before diagnosing it.
    #[arg(long, default_value_t = 20)]
    min_family_support: usize,

    /// Minimum reads required for a candidate 5 prime mode.
    #[arg(long, default_value_t = 20)]
    min_mode_support: usize,

    /// Minimum fraction of the 3 prime family required for a candidate 5 prime mode.
    #[arg(long, default_value_t = 0.02)]
    min_mode_fraction: f64,

    /// Maximum number of candidate modes to emit per 3 prime family.
    #[arg(long, default_value_t = 3)]
    max_candidates_per_family: usize,

    /// Optional minimum read length filter (bp).
    #[arg(long)]
    min_read_len: Option<u32>,
}

#[derive(Parser, Debug, Clone)]
#[command(
    name = "trackclustertu rescue-missed-tus",
    version,
    about = "Promote high-support boundary modes into rescued TU calls"
)]
struct RescueMissedTusCli {
    /// Input BED6 reads file used for clustering.
    #[arg(long = "in")]
    input: PathBuf,

    /// Existing TU BED6 output (for example, `results/tus.bed`).
    #[arg(long = "existing-tu")]
    existing_tu: PathBuf,

    /// Existing membership TSV from the original cluster run.
    #[arg(long = "existing-membership")]
    existing_membership: PathBuf,

    /// Optional gene annotation BED6 used to label rescued candidates.
    #[arg(long)]
    annotation_bed: Option<PathBuf>,

    /// Output rescued TU BED6.
    #[arg(long = "out-tu")]
    out_tu: PathBuf,

    /// Output rescued membership TSV.
    #[arg(long = "out-membership")]
    out_membership: PathBuf,

    /// Optional rescued TU counts CSV (`tu_id,count`).
    #[arg(long = "out-tu-count")]
    out_tu_count: Option<PathBuf>,

    /// Optional TSV report of the rescued candidate modes.
    #[arg(long = "out-candidates-tsv")]
    out_candidates_tsv: Option<PathBuf>,

    /// Optional BED6 track of the rescued candidate modes.
    #[arg(long = "out-candidates-bed")]
    out_candidates_bed: Option<PathBuf>,

    /// Prefix used when naming rescued TU IDs.
    #[arg(long, default_value = "RESC")]
    rescue_prefix: String,

    /// Strand-aware 3 prime window used to form end families and TU matches (bp).
    #[arg(long, default_value_t = 12)]
    three_prime_window_bp: u32,

    /// Strand-aware 5 prime window used to call start modes and TU matches (bp).
    #[arg(long, default_value_t = 10)]
    five_prime_window_bp: u32,

    /// Minimum reads required in a 3 prime family before diagnosing it.
    #[arg(long, default_value_t = 20)]
    min_family_support: usize,

    /// Minimum reads required for a candidate 5 prime mode.
    #[arg(long, default_value_t = 20)]
    min_mode_support: usize,

    /// Minimum fraction of the 3 prime family required for a candidate 5 prime mode.
    #[arg(long, default_value_t = 0.02)]
    min_mode_fraction: f64,

    /// Maximum number of candidate modes to emit per 3 prime family.
    #[arg(long, default_value_t = 3)]
    max_candidates_per_family: usize,

    /// Optional minimum read length filter (bp).
    #[arg(long)]
    min_read_len: Option<u32>,
}

#[derive(Clone, Debug)]
struct ExistingTu {
    id: String,
    contig: String,
    strand: Strand,
    interval: Interval,
}

#[derive(Clone, Debug)]
struct GeneRecord {
    id: String,
    contig: String,
    strand: Strand,
    interval: Interval,
}

#[derive(Clone, Debug)]
struct ThreePrimeFamily {
    contig: String,
    strand: Strand,
    three_prime_coord: u32,
    read_indices: Vec<usize>,
}

#[derive(Clone, Debug)]
struct BoundaryPeak {
    five_prime_coord: u32,
    support: usize,
}

#[derive(Clone, Debug)]
struct BestExistingTu {
    id: String,
    interval: Interval,
    three_prime_delta_bp: u32,
    five_prime_delta_bp: u32,
    score1: f64,
    score2: f64,
}

#[derive(Clone, Debug)]
struct CandidateMode {
    report_id: String,
    contig: String,
    strand: Strand,
    interval: Interval,
    three_prime_coord: u32,
    five_prime_coord: u32,
    family_support: usize,
    mode_support: usize,
    mode_fraction: f64,
    reason: String,
    nearest_tu: Option<BestExistingTu>,
    gene_hint: String,
    support_read_indices: Vec<usize>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct BoundaryQuality {
    three_prime_delta_bp: u32,
    five_prime_delta_bp: u32,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BoundaryKind {
    ThreePrime,
    FivePrime,
}

#[derive(Clone, Copy, Debug)]
struct DetectionParams {
    three_prime_window_bp: u32,
    five_prime_window_bp: u32,
    min_family_support: usize,
    min_mode_support: usize,
    min_mode_fraction: f64,
    max_candidates_per_family: usize,
}

pub(crate) fn run_from_args<I, T>(args: I) -> Result<()>
where
    I: IntoIterator<Item = T>,
    T: Into<std::ffi::OsString> + Clone,
{
    let cli = DiagnoseMissedTusCli::parse_from(args);
    run(&cli)
}

pub(crate) fn run_rescue_from_args<I, T>(args: I) -> Result<()>
where
    I: IntoIterator<Item = T>,
    T: Into<std::ffi::OsString> + Clone,
{
    let cli = RescueMissedTusCli::parse_from(args);
    run_rescue(&cli)
}

fn run(cli: &DiagnoseMissedTusCli) -> Result<()> {
    let params = detection_params(
        cli.three_prime_window_bp,
        cli.five_prime_window_bp,
        cli.min_family_support,
        cli.min_mode_support,
        cli.min_mode_fraction,
        cli.max_candidates_per_family,
    )?;
    let (reads, existing_tus, genes) = load_detection_inputs(
        &cli.input,
        &cli.existing_tu,
        cli.annotation_bed.as_deref(),
        cli.min_read_len,
    )?;
    let mut candidates = detect_candidate_modes(&reads, &existing_tus, &genes, params);
    assign_candidate_ids(&mut candidates, "MISS");

    write_tsv_report(&cli.out_tsv, &candidates)?;
    if let Some(out_bed) = cli.out_bed.as_ref() {
        write_bed_report(out_bed, &candidates)?;
    }

    println!("candidate_count={}", candidates.len());
    println!("out_tsv={}", cli.out_tsv.display());
    if let Some(out_bed) = cli.out_bed.as_ref() {
        println!("out_bed={}", out_bed.display());
    }

    Ok(())
}

fn run_rescue(cli: &RescueMissedTusCli) -> Result<()> {
    let params = detection_params(
        cli.three_prime_window_bp,
        cli.five_prime_window_bp,
        cli.min_family_support,
        cli.min_mode_support,
        cli.min_mode_fraction,
        cli.max_candidates_per_family,
    )?;
    let (reads, existing_tus, genes) = load_detection_inputs(
        &cli.input,
        &cli.existing_tu,
        cli.annotation_bed.as_deref(),
        cli.min_read_len,
    )?;
    let candidate_modes = detect_candidate_modes(&reads, &existing_tus, &genes, params);
    let candidates = deduplicate_candidates(candidate_modes, params);

    let read_to_existing = read_existing_membership(&cli.existing_membership)?;
    let read_lookup: HashMap<&str, usize> = reads
        .iter()
        .enumerate()
        .map(|(idx, read)| (read.id.as_str(), idx))
        .collect();
    let tu_lookup: HashMap<&str, &ExistingTu> =
        existing_tus.iter().map(|tu| (tu.id.as_str(), tu)).collect();

    let mut best_candidate_for_read: Vec<Option<(usize, BoundaryQuality)>> =
        vec![None; reads.len()];
    for (candidate_idx, candidate) in candidates.iter().enumerate() {
        for &read_idx in &candidate.support_read_indices {
            let quality = boundary_quality_for_read(&reads[read_idx], candidate.interval);
            let replace = match best_candidate_for_read[read_idx] {
                None => true,
                Some((current_idx, current_quality)) => better_candidate_choice(
                    candidate_idx,
                    quality,
                    candidate.mode_support,
                    current_idx,
                    current_quality,
                    candidates[current_idx].mode_support,
                ),
            };
            if replace {
                best_candidate_for_read[read_idx] = Some((candidate_idx, quality));
            }
        }
    }

    let mut provisional_candidate_counts = vec![0usize; candidates.len()];
    for (candidate_idx, _) in best_candidate_for_read.iter().flatten() {
        provisional_candidate_counts[*candidate_idx] += 1;
    }
    let candidate_kept: Vec<bool> = provisional_candidate_counts
        .iter()
        .map(|&count| count >= cli.min_mode_support)
        .collect();

    let mut final_read_to_candidate: Vec<Option<usize>> = vec![None; reads.len()];
    let mut final_candidate_counts = vec![0usize; candidates.len()];
    for (read_idx, read) in reads.iter().enumerate() {
        let Some((candidate_idx, candidate_quality)) = best_candidate_for_read[read_idx] else {
            continue;
        };
        if !candidate_kept[candidate_idx] {
            continue;
        }

        let current_tu = read_to_existing
            .get(read.id.as_str())
            .and_then(|tu_id| tu_lookup.get(tu_id.as_str()).copied());
        let use_candidate = match current_tu {
            Some(tu) => {
                let current_quality = boundary_quality_for_read(read, tu.interval);
                better_boundary_quality(
                    candidate_quality,
                    candidates[candidate_idx].mode_support,
                    current_quality,
                    0,
                )
            }
            None => true,
        };

        if use_candidate {
            final_read_to_candidate[read_idx] = Some(candidate_idx);
            final_candidate_counts[candidate_idx] += 1;
        }
    }

    let final_candidate_kept: Vec<bool> = final_candidate_counts
        .iter()
        .map(|&count| count >= cli.min_mode_support)
        .collect();
    let kept_candidate_index_map = build_kept_candidate_index_map(&final_candidate_kept);
    let mut kept_candidates: Vec<CandidateMode> = Vec::new();
    for (candidate_idx, candidate) in candidates.into_iter().enumerate() {
        if final_candidate_kept[candidate_idx] {
            kept_candidates.push(candidate);
        }
    }
    let mut candidates = kept_candidates;

    let rescue_width = 4usize.max(candidates.len().to_string().len());
    let mut rescue_ids: Vec<String> = Vec::with_capacity(candidates.len());
    for idx in 0..candidates.len() {
        rescue_ids.push(format!(
            "{}{:0width$}",
            cli.rescue_prefix,
            idx + 1,
            width = rescue_width
        ));
    }

    let mut final_assignments: Vec<Option<(&str, Interval)>> = vec![None; reads.len()];
    for (read_id, tu_id) in &read_to_existing {
        if let Some(&read_idx) = read_lookup.get(read_id.as_str()) {
            if let Some(tu) = tu_lookup.get(tu_id.as_str()) {
                final_assignments[read_idx] = Some((tu.id.as_str(), tu.interval));
            }
        }
    }
    for read_idx in 0..reads.len() {
        let Some(original_candidate_idx) = final_read_to_candidate[read_idx] else {
            continue;
        };
        let Some(new_candidate_idx) = kept_candidate_index_map[original_candidate_idx] else {
            continue;
        };
        final_assignments[read_idx] = Some((
            rescue_ids[new_candidate_idx].as_str(),
            candidates[new_candidate_idx].interval,
        ));
    }

    let mut existing_final_counts: HashMap<String, u64> = HashMap::new();
    let mut rescue_final_counts: Vec<u64> = vec![0; candidates.len()];
    let rescue_lookup: HashMap<&str, usize> = rescue_ids
        .iter()
        .enumerate()
        .map(|(idx, id)| (id.as_str(), idx))
        .collect();
    for assignment in &final_assignments {
        let Some((tu_id, _)) = assignment else {
            continue;
        };
        if let Some(&rescue_idx) = rescue_lookup.get(tu_id) {
            rescue_final_counts[rescue_idx] += 1;
        } else {
            *existing_final_counts
                .entry((*tu_id).to_owned())
                .or_default() += 1;
        }
    }

    let mut final_tus: Vec<(String, String, Strand, Interval)> = Vec::new();
    for tu in &existing_tus {
        if existing_final_counts.get(&tu.id).copied().unwrap_or(0) > 0 {
            final_tus.push((tu.id.clone(), tu.contig.clone(), tu.strand, tu.interval));
        }
    }
    for (candidate_idx, candidate) in candidates.iter().enumerate() {
        if rescue_final_counts[candidate_idx] > 0 {
            final_tus.push((
                rescue_ids[candidate_idx].clone(),
                candidate.contig.clone(),
                candidate.strand,
                candidate.interval,
            ));
        }
    }
    final_tus.sort_by(cmp_interval_record);

    write_rescued_tu_bed(&cli.out_tu, &final_tus)?;
    write_rescued_membership(&cli.out_membership, &reads, &final_assignments)?;
    if let Some(out_tu_count) = cli.out_tu_count.as_ref() {
        write_rescued_counts_csv(
            out_tu_count,
            &final_tus,
            &existing_final_counts,
            &rescue_ids,
            &rescue_final_counts,
        )?;
    }

    assign_candidate_ids(&mut candidates, "MISS");
    if let Some(out_tsv) = cli.out_candidates_tsv.as_ref() {
        write_tsv_report(out_tsv, &candidates)?;
    }
    if let Some(out_bed) = cli.out_candidates_bed.as_ref() {
        write_bed_report(out_bed, &candidates)?;
    }

    println!("rescued_candidate_count={}", candidates.len());
    println!("out_tu={}", cli.out_tu.display());
    println!("out_membership={}", cli.out_membership.display());
    if let Some(out_tu_count) = cli.out_tu_count.as_ref() {
        println!("out_tu_count={}", out_tu_count.display());
    }

    Ok(())
}

fn detection_params(
    three_prime_window_bp: u32,
    five_prime_window_bp: u32,
    min_family_support: usize,
    min_mode_support: usize,
    min_mode_fraction: f64,
    max_candidates_per_family: usize,
) -> Result<DetectionParams> {
    if !(0.0..=1.0).contains(&min_mode_fraction) {
        anyhow::bail!(
            "--min-mode-fraction must be between 0 and 1, got {}",
            min_mode_fraction
        );
    }
    if max_candidates_per_family == 0 {
        anyhow::bail!("--max-candidates-per-family must be >= 1");
    }
    Ok(DetectionParams {
        three_prime_window_bp,
        five_prime_window_bp,
        min_family_support,
        min_mode_support,
        min_mode_fraction,
        max_candidates_per_family,
    })
}

fn load_detection_inputs(
    input: &Path,
    existing_tu: &Path,
    annotation_bed: Option<&Path>,
    min_read_len: Option<u32>,
) -> Result<(Vec<ReadRecord>, Vec<ExistingTu>, Vec<GeneRecord>)> {
    let mut reads = read_reads_bed6(input)?;
    if let Some(min_len) = min_read_len {
        reads.retain(|read| read.interval.len() >= min_len);
    }
    let existing_tus = read_existing_tus_bed6(existing_tu)?;
    let genes = match annotation_bed {
        Some(path) => read_gene_bed6(path)?,
        None => Vec::new(),
    };
    Ok((reads, existing_tus, genes))
}

fn detect_candidate_modes(
    reads: &[ReadRecord],
    existing_tus: &[ExistingTu],
    genes: &[GeneRecord],
    params: DetectionParams,
) -> Vec<CandidateMode> {
    let mut tus_by_partition: HashMap<(String, Strand), Vec<ExistingTu>> = HashMap::new();
    for tu in existing_tus {
        tus_by_partition
            .entry((tu.contig.clone(), tu.strand))
            .or_default()
            .push(tu.clone());
    }

    let mut genes_by_partition: HashMap<(String, Strand), Vec<GeneRecord>> = HashMap::new();
    for gene in genes {
        genes_by_partition
            .entry((gene.contig.clone(), gene.strand))
            .or_default()
            .push(gene.clone());
    }

    let families = build_three_prime_families(reads, params.three_prime_window_bp);
    let mut candidates: Vec<CandidateMode> = Vec::new();

    for family in families {
        let family_support = family.read_indices.len();
        if family_support < params.min_family_support {
            continue;
        }

        let partition_key = (family.contig.clone(), family.strand);
        let partition_tus = tus_by_partition
            .get(&partition_key)
            .map(Vec::as_slice)
            .unwrap_or(&[]);
        let partition_genes = genes_by_partition
            .get(&partition_key)
            .map(Vec::as_slice)
            .unwrap_or(&[]);

        let family_reads: Vec<&ReadRecord> =
            family.read_indices.iter().map(|&idx| &reads[idx]).collect();
        let mut peaks =
            detect_five_prime_modes(&family_reads, family.strand, params.five_prime_window_bp);
        peaks.sort_by(|a, b| {
            b.support.cmp(&a.support).then_with(|| {
                boundary_priority(family.strand, BoundaryKind::FivePrime, a.five_prime_coord)
                    .cmp(&boundary_priority(
                        family.strand,
                        BoundaryKind::FivePrime,
                        b.five_prime_coord,
                    ))
                    .reverse()
            })
        });

        for peak in peaks.into_iter().take(params.max_candidates_per_family) {
            let support_read_indices: Vec<usize> = family
                .read_indices
                .iter()
                .copied()
                .filter(|&read_idx| {
                    let read = &reads[read_idx];
                    three_prime_coord(read).abs_diff(family.three_prime_coord)
                        <= params.three_prime_window_bp
                        && five_prime_coord(read).abs_diff(peak.five_prime_coord)
                            <= params.five_prime_window_bp
                })
                .collect();
            let mode_support = support_read_indices.len();
            if mode_support < params.min_mode_support {
                continue;
            }

            let mode_fraction = mode_support as f64 / family_support as f64;
            if mode_fraction < params.min_mode_fraction {
                continue;
            }

            let Some(interval) = interval_from_boundaries(
                family.strand,
                peak.five_prime_coord,
                family.three_prime_coord,
            ) else {
                continue;
            };

            let best_match = best_existing_tu_match(
                &family.contig,
                family.strand,
                interval,
                partition_tus,
                params.three_prime_window_bp,
                params.five_prime_window_bp,
            );
            if best_match.represented {
                continue;
            }

            candidates.push(CandidateMode {
                report_id: String::new(),
                contig: family.contig.clone(),
                strand: family.strand,
                interval,
                three_prime_coord: family.three_prime_coord,
                five_prime_coord: peak.five_prime_coord,
                family_support,
                mode_support,
                mode_fraction,
                reason: infer_reason(&best_match),
                nearest_tu: best_match.nearest,
                gene_hint: gene_hint(interval, partition_genes),
                support_read_indices,
            });
        }
    }

    candidates
}

fn assign_candidate_ids(candidates: &mut [CandidateMode], prefix: &str) {
    candidates.sort_by(|a, b| {
        a.contig
            .cmp(&b.contig)
            .then_with(|| a.strand.cmp(&b.strand))
            .then_with(|| a.interval.start.cmp(&b.interval.start))
            .then_with(|| a.interval.end.cmp(&b.interval.end))
            .then_with(|| b.mode_support.cmp(&a.mode_support))
    });

    let width = 4usize.max(candidates.len().to_string().len());
    for (idx, candidate) in candidates.iter_mut().enumerate() {
        candidate.report_id = format!("{prefix}{:0width$}", idx + 1, width = width);
    }
}

fn candidate_report_id(candidate: &CandidateMode) -> &str {
    candidate.report_id.as_str()
}

fn candidate_reason(candidate: &CandidateMode) -> &str {
    candidate.reason.as_str()
}

fn deduplicate_candidates(
    mut candidates: Vec<CandidateMode>,
    params: DetectionParams,
) -> Vec<CandidateMode> {
    candidates.sort_by(|a, b| {
        b.mode_support
            .cmp(&a.mode_support)
            .then_with(|| b.mode_fraction.total_cmp(&a.mode_fraction))
            .then_with(|| a.contig.cmp(&b.contig))
            .then_with(|| a.strand.cmp(&b.strand))
            .then_with(|| a.interval.start.cmp(&b.interval.start))
            .then_with(|| a.interval.end.cmp(&b.interval.end))
    });

    let mut kept: Vec<CandidateMode> = Vec::new();
    'candidate: for candidate in candidates {
        for existing in &kept {
            if existing.contig != candidate.contig || existing.strand != candidate.strand {
                continue;
            }
            let existing_quality = boundary_quality_for_interval(
                candidate.interval,
                existing.interval,
                candidate.strand,
            );
            if existing_quality.three_prime_delta_bp <= params.three_prime_window_bp
                && existing_quality.five_prime_delta_bp <= params.five_prime_window_bp
            {
                continue 'candidate;
            }
        }
        kept.push(candidate);
    }

    kept.sort_by(|a, b| {
        a.contig
            .cmp(&b.contig)
            .then_with(|| a.strand.cmp(&b.strand))
            .then_with(|| a.interval.start.cmp(&b.interval.start))
            .then_with(|| a.interval.end.cmp(&b.interval.end))
            .then_with(|| b.mode_support.cmp(&a.mode_support))
    });
    kept
}

fn read_existing_membership(path: &Path) -> Result<HashMap<String, String>> {
    let mut assignments = HashMap::new();
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read membership {:?}", path))?;
    for (line_idx, line) in text.lines().enumerate() {
        let line_number = line_idx + 1;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 2 {
            anyhow::bail!(
                "{path:?}:{line_number}: expected at least 2 columns (read_id, tu_id), got {}",
                fields.len()
            );
        }
        assignments.insert(fields[0].to_owned(), fields[1].to_owned());
    }
    Ok(assignments)
}

fn boundary_quality_for_read(read: &ReadRecord, interval: Interval) -> BoundaryQuality {
    BoundaryQuality {
        three_prime_delta_bp: three_prime_coord(read)
            .abs_diff(three_prime_coord_interval(interval, read.strand)),
        five_prime_delta_bp: five_prime_coord(read)
            .abs_diff(five_prime_coord_interval(interval, read.strand)),
    }
}

fn boundary_quality_for_interval(a: Interval, b: Interval, strand: Strand) -> BoundaryQuality {
    BoundaryQuality {
        three_prime_delta_bp: three_prime_coord_interval(a, strand)
            .abs_diff(three_prime_coord_interval(b, strand)),
        five_prime_delta_bp: five_prime_coord_interval(a, strand)
            .abs_diff(five_prime_coord_interval(b, strand)),
    }
}

fn better_boundary_quality(
    candidate_quality: BoundaryQuality,
    candidate_support: usize,
    current_quality: BoundaryQuality,
    current_support: usize,
) -> bool {
    (
        candidate_quality.three_prime_delta_bp,
        candidate_quality.five_prime_delta_bp,
    ) < (
        current_quality.three_prime_delta_bp,
        current_quality.five_prime_delta_bp,
    ) || ((
        candidate_quality.three_prime_delta_bp,
        candidate_quality.five_prime_delta_bp,
    ) == (
        current_quality.three_prime_delta_bp,
        current_quality.five_prime_delta_bp,
    ) && candidate_support > current_support)
}

fn better_candidate_choice(
    candidate_idx: usize,
    candidate_quality: BoundaryQuality,
    candidate_support: usize,
    current_idx: usize,
    current_quality: BoundaryQuality,
    current_support: usize,
) -> bool {
    better_boundary_quality(
        candidate_quality,
        candidate_support,
        current_quality,
        current_support,
    ) || (candidate_quality == current_quality
        && candidate_support == current_support
        && candidate_idx < current_idx)
}

fn build_kept_candidate_index_map(final_candidate_kept: &[bool]) -> Vec<Option<usize>> {
    let mut kept_index_map = vec![None; final_candidate_kept.len()];
    let mut next_idx = 0usize;
    for (idx, kept) in final_candidate_kept.iter().copied().enumerate() {
        if kept {
            kept_index_map[idx] = Some(next_idx);
            next_idx += 1;
        }
    }
    kept_index_map
}

fn cmp_interval_record(
    a: &(String, String, Strand, Interval),
    b: &(String, String, Strand, Interval),
) -> Ordering {
    a.1.cmp(&b.1)
        .then_with(|| a.2.cmp(&b.2))
        .then_with(|| a.3.start.cmp(&b.3.start))
        .then_with(|| a.3.end.cmp(&b.3.end))
        .then_with(|| a.0.cmp(&b.0))
}

fn cmp_read(a: &ReadRecord, b: &ReadRecord) -> Ordering {
    a.contig
        .cmp(&b.contig)
        .then_with(|| a.strand.cmp(&b.strand))
        .then_with(|| a.interval.start.cmp(&b.interval.start))
        .then_with(|| a.interval.end.cmp(&b.interval.end))
        .then_with(|| a.id.cmp(&b.id))
}

fn write_rescued_tu_bed(
    path: &Path,
    final_tus: &[(String, String, Strand, Interval)],
) -> Result<()> {
    let records: Vec<Bed6Record> = final_tus
        .iter()
        .map(|(id, contig, strand, interval)| Bed6Record {
            chrom: contig.clone(),
            start: interval.start,
            end: interval.end,
            name: id.clone(),
            score: 0,
            strand: *strand,
            extra_fields: Vec::new(),
        })
        .collect();
    crate::io::bed::write_bed6(path, records.iter())
        .with_context(|| format!("failed to write rescued TU BED {:?}", path))?;
    Ok(())
}

fn write_rescued_membership(
    path: &Path,
    reads: &[ReadRecord],
    assignments: &[Option<(&str, Interval)>],
) -> Result<()> {
    let mut read_indices: Vec<usize> = (0..reads.len()).collect();
    read_indices.sort_unstable_by(|&a, &b| cmp_read(&reads[a], &reads[b]).then_with(|| a.cmp(&b)));

    let mut writer = BufWriter::new(
        File::create(path).with_context(|| format!("failed to create membership {:?}", path))?,
    );
    for read_idx in read_indices {
        let Some((tu_id, interval)) = assignments[read_idx] else {
            continue;
        };
        let read = &reads[read_idx];
        let score1 = score1_interval(read.interval, interval);
        let score2 = score2_interval(read.interval, interval);
        writeln!(writer, "{}\t{}\t{score1:.6}\t{score2:.6}", read.id, tu_id)?;
    }
    Ok(())
}

fn write_rescued_counts_csv(
    path: &Path,
    final_tus: &[(String, String, Strand, Interval)],
    existing_final_counts: &HashMap<String, u64>,
    rescue_ids: &[String],
    rescue_final_counts: &[u64],
) -> Result<()> {
    let rescue_count_map: HashMap<&str, u64> = rescue_ids
        .iter()
        .zip(rescue_final_counts.iter().copied())
        .map(|(id, count)| (id.as_str(), count))
        .collect();
    let mut writer = BufWriter::new(
        File::create(path).with_context(|| format!("failed to create TU counts {:?}", path))?,
    );
    writeln!(writer, "tu_id,count")?;
    for (tu_id, _, _, _) in final_tus {
        let count = existing_final_counts
            .get(tu_id)
            .copied()
            .or_else(|| rescue_count_map.get(tu_id.as_str()).copied())
            .unwrap_or(0);
        writeln!(writer, "{tu_id},{count}")?;
    }
    Ok(())
}

fn read_reads_bed6(path: &Path) -> Result<Vec<ReadRecord>> {
    let reader = crate::io::bed::read_bed6(path)
        .with_context(|| format!("failed to open BED6 reads {:?}", path))?;
    reader
        .map(|record| {
            let record =
                record.with_context(|| format!("failed to parse BED6 reads {:?}", path))?;
            let interval = Interval::new(record.start, record.end)
                .with_context(|| format!("invalid BED6 interval in {:?}", path))?;
            Ok(ReadRecord {
                contig: record.chrom,
                strand: record.strand,
                interval,
                id: record.name,
            })
        })
        .collect()
}

fn read_existing_tus_bed6(path: &Path) -> Result<Vec<ExistingTu>> {
    let reader = crate::io::bed::read_bed6(path)
        .with_context(|| format!("failed to open TU BED6 {:?}", path))?;
    reader
        .map(|record| {
            let record = record.with_context(|| format!("failed to parse TU BED6 {:?}", path))?;
            let interval = Interval::new(record.start, record.end)
                .with_context(|| format!("invalid TU interval in {:?}", path))?;
            Ok(ExistingTu {
                id: record.name,
                contig: record.chrom,
                strand: record.strand,
                interval,
            })
        })
        .collect()
}

fn read_gene_bed6(path: &Path) -> Result<Vec<GeneRecord>> {
    let reader = crate::io::bed::read_bed6(path)
        .with_context(|| format!("failed to open annotation BED6 {:?}", path))?;
    reader
        .map(|record| {
            let record =
                record.with_context(|| format!("failed to parse annotation BED6 {:?}", path))?;
            let interval = Interval::new(record.start, record.end)
                .with_context(|| format!("invalid gene interval in {:?}", path))?;
            Ok(GeneRecord {
                id: record.name,
                contig: record.chrom,
                strand: record.strand,
                interval,
            })
        })
        .collect()
}

fn build_three_prime_families(reads: &[ReadRecord], window_bp: u32) -> Vec<ThreePrimeFamily> {
    let mut indices: Vec<usize> = (0..reads.len()).collect();
    indices.sort_unstable_by(|&a, &b| {
        reads[a]
            .contig
            .cmp(&reads[b].contig)
            .then_with(|| reads[a].strand.cmp(&reads[b].strand))
            .then_with(|| three_prime_coord(&reads[a]).cmp(&three_prime_coord(&reads[b])))
            .then_with(|| five_prime_coord(&reads[a]).cmp(&five_prime_coord(&reads[b])))
            .then_with(|| reads[a].id.cmp(&reads[b].id))
    });

    let mut families: Vec<ThreePrimeFamily> = Vec::new();
    let mut family_indices: Vec<usize> = Vec::new();
    let mut family_contig: Option<String> = None;
    let mut family_strand = Strand::Unknown;
    let mut prev_three_prime: Option<u32> = None;

    for read_idx in indices {
        let read = &reads[read_idx];
        let coord = three_prime_coord(read);
        let same_partition = family_contig
            .as_deref()
            .map(|contig| contig == read.contig.as_str() && family_strand == read.strand)
            .unwrap_or(false);
        let same_family = same_partition
            && prev_three_prime
                .map(|prev| coord.saturating_sub(prev) <= window_bp)
                .unwrap_or(false);

        if !family_indices.is_empty() && !same_family {
            families.push(finalize_three_prime_family(
                reads,
                std::mem::take(&mut family_indices),
            ));
        }

        if family_indices.is_empty() {
            family_contig = Some(read.contig.clone());
            family_strand = read.strand;
        }

        family_indices.push(read_idx);
        prev_three_prime = Some(coord);
    }

    if !family_indices.is_empty() {
        families.push(finalize_three_prime_family(reads, family_indices));
    }

    families
}

fn finalize_three_prime_family(reads: &[ReadRecord], read_indices: Vec<usize>) -> ThreePrimeFamily {
    let first = &reads[read_indices[0]];
    let three_prime_coord = choose_mode_coordinate(
        read_indices
            .iter()
            .map(|&idx| three_prime_coord(&reads[idx])),
        first.strand,
        BoundaryKind::ThreePrime,
    );

    ThreePrimeFamily {
        contig: first.contig.clone(),
        strand: first.strand,
        three_prime_coord,
        read_indices,
    }
}

fn detect_five_prime_modes(
    reads: &[&ReadRecord],
    strand: Strand,
    window_bp: u32,
) -> Vec<BoundaryPeak> {
    let mut counts: BTreeMap<u32, usize> = BTreeMap::new();
    for read in reads {
        *counts.entry(five_prime_coord(read)).or_default() += 1;
    }

    let mut remaining: Vec<(u32, usize)> = counts.into_iter().collect();
    let mut peaks: Vec<BoundaryPeak> = Vec::new();

    while !remaining.is_empty() {
        let mut best_idx = 0usize;
        let mut best_support = 0usize;
        let mut best_count = 0usize;

        for i in 0..remaining.len() {
            let center = remaining[i].0;
            let support = remaining
                .iter()
                .filter(|(coord, _)| center.abs_diff(*coord) <= window_bp)
                .map(|(_, count)| *count)
                .sum::<usize>();
            let count = remaining[i].1;

            let better = support > best_support
                || (support == best_support && count > best_count)
                || (support == best_support
                    && count == best_count
                    && boundary_priority(strand, BoundaryKind::FivePrime, center)
                        > boundary_priority(
                            strand,
                            BoundaryKind::FivePrime,
                            remaining[best_idx].0,
                        ));
            if better {
                best_idx = i;
                best_support = support;
                best_count = count;
            }
        }

        let center = remaining[best_idx].0;
        let mut next_remaining: Vec<(u32, usize)> = Vec::new();
        let mut support = 0usize;
        for (coord, count) in remaining {
            if center.abs_diff(coord) <= window_bp {
                support += count;
            } else {
                next_remaining.push((coord, count));
            }
        }
        peaks.push(BoundaryPeak {
            five_prime_coord: center,
            support,
        });
        remaining = next_remaining;
    }

    peaks
}

#[derive(Default)]
struct MatchSummary {
    represented: bool,
    any_same_three_prime: bool,
    any_same_five_prime: bool,
    nearest: Option<BestExistingTu>,
}

fn best_existing_tu_match(
    contig: &str,
    strand: Strand,
    candidate: Interval,
    existing_tus: &[ExistingTu],
    three_prime_window_bp: u32,
    five_prime_window_bp: u32,
) -> MatchSummary {
    let candidate_three_prime = three_prime_coord_interval(candidate, strand);
    let candidate_five_prime = five_prime_coord_interval(candidate, strand);
    let mut summary = MatchSummary::default();

    for tu in existing_tus {
        if tu.contig != contig || tu.strand != strand {
            continue;
        }

        let three_prime_delta_bp =
            candidate_three_prime.abs_diff(three_prime_coord_interval(tu.interval, strand));
        let five_prime_delta_bp =
            candidate_five_prime.abs_diff(five_prime_coord_interval(tu.interval, strand));
        let score1 = score1_interval(candidate, tu.interval);
        let score2 = score2_interval(candidate, tu.interval);

        if three_prime_delta_bp <= three_prime_window_bp {
            summary.any_same_three_prime = true;
        }
        if five_prime_delta_bp <= five_prime_window_bp {
            summary.any_same_five_prime = true;
        }
        if three_prime_delta_bp <= three_prime_window_bp
            && five_prime_delta_bp <= five_prime_window_bp
        {
            summary.represented = true;
        }

        let candidate_match = BestExistingTu {
            id: tu.id.clone(),
            interval: tu.interval,
            three_prime_delta_bp,
            five_prime_delta_bp,
            score1,
            score2,
        };

        let replace = match summary.nearest.as_ref() {
            None => true,
            Some(current) => {
                (
                    candidate_match.three_prime_delta_bp,
                    candidate_match.five_prime_delta_bp,
                ) < (current.three_prime_delta_bp, current.five_prime_delta_bp)
                    || ((
                        candidate_match.three_prime_delta_bp,
                        candidate_match.five_prime_delta_bp,
                    ) == (current.three_prime_delta_bp, current.five_prime_delta_bp)
                        && candidate_match.score2.total_cmp(&current.score2) == Ordering::Greater)
            }
        };
        if replace {
            summary.nearest = Some(candidate_match);
        }
    }

    summary
}

fn infer_reason(match_summary: &MatchSummary) -> String {
    if match_summary.any_same_three_prime && !match_summary.any_same_five_prime {
        "same_3p_missing_5p_mode".to_owned()
    } else if match_summary.any_same_five_prime && !match_summary.any_same_three_prime {
        "same_5p_missing_3p_mode".to_owned()
    } else if match_summary.any_same_three_prime && match_summary.any_same_five_prime {
        "boundary_mode_shifted".to_owned()
    } else {
        "no_matching_tu".to_owned()
    }
}

fn gene_hint(candidate: Interval, genes: &[GeneRecord]) -> String {
    let overlaps: Vec<&str> = genes
        .iter()
        .filter(|gene| gene.interval.overlaps(candidate))
        .map(|gene| gene.id.as_str())
        .collect();
    if !overlaps.is_empty() {
        return overlaps.join(",");
    }

    let mut nearest: Option<(&str, u32)> = None;
    for gene in genes {
        let distance = interval_distance(candidate, gene.interval);
        let replace = match nearest {
            None => true,
            Some((_, best_distance)) => distance < best_distance,
        };
        if replace {
            nearest = Some((gene.id.as_str(), distance));
        }
    }

    nearest
        .map(|(id, _)| id.to_owned())
        .unwrap_or_else(|| ".".to_owned())
}

fn interval_distance(a: Interval, b: Interval) -> u32 {
    if a.overlaps(b) {
        0
    } else if a.end <= b.start {
        b.start.get().saturating_sub(a.end.get())
    } else {
        a.start.get().saturating_sub(b.end.get())
    }
}

fn interval_from_boundaries(strand: Strand, five_prime: u32, three_prime: u32) -> Option<Interval> {
    let start = match strand {
        Strand::Plus | Strand::Unknown => five_prime,
        Strand::Minus => three_prime,
    };
    let end = match strand {
        Strand::Plus | Strand::Unknown => three_prime,
        Strand::Minus => five_prime,
    };
    Interval::new(Coord::new(start), Coord::new(end)).ok()
}

fn three_prime_coord(read: &ReadRecord) -> u32 {
    three_prime_coord_interval(read.interval, read.strand)
}

fn five_prime_coord(read: &ReadRecord) -> u32 {
    five_prime_coord_interval(read.interval, read.strand)
}

fn three_prime_coord_interval(interval: Interval, strand: Strand) -> u32 {
    match strand {
        Strand::Plus | Strand::Unknown => interval.end.get(),
        Strand::Minus => interval.start.get(),
    }
}

fn five_prime_coord_interval(interval: Interval, strand: Strand) -> u32 {
    match strand {
        Strand::Plus | Strand::Unknown => interval.start.get(),
        Strand::Minus => interval.end.get(),
    }
}

fn choose_mode_coordinate<I>(coords: I, strand: Strand, boundary: BoundaryKind) -> u32
where
    I: IntoIterator<Item = u32>,
{
    let mut counts: BTreeMap<u32, usize> = BTreeMap::new();
    for coord in coords {
        *counts.entry(coord).or_default() += 1;
    }

    counts
        .into_iter()
        .max_by(|(coord_a, count_a), (coord_b, count_b)| {
            count_a.cmp(count_b).then_with(|| {
                boundary_priority(strand, boundary, *coord_a)
                    .cmp(&boundary_priority(strand, boundary, *coord_b))
            })
        })
        .map(|(coord, _)| coord)
        .unwrap_or(0)
}

fn boundary_priority(strand: Strand, boundary: BoundaryKind, coord: u32) -> i64 {
    let coord = coord as i64;
    match (strand, boundary) {
        (Strand::Plus, BoundaryKind::ThreePrime) => coord,
        (Strand::Minus, BoundaryKind::ThreePrime) => -coord,
        (Strand::Unknown, BoundaryKind::ThreePrime) => coord,
        (Strand::Plus, BoundaryKind::FivePrime) => -coord,
        (Strand::Minus, BoundaryKind::FivePrime) => coord,
        (Strand::Unknown, BoundaryKind::FivePrime) => -coord,
    }
}

fn write_tsv_report(path: &Path, candidates: &[CandidateMode]) -> Result<()> {
    let mut writer = BufWriter::new(
        File::create(path).with_context(|| format!("failed to create report {:?}", path))?,
    );
    writeln!(
        writer,
        "candidate_id\tcontig\tstrand\tcandidate_start\tcandidate_end\tthree_prime_coord\tfive_prime_coord\tfamily_support\tmode_support\tmode_fraction\treason\tnearest_tu_id\tnearest_tu_start\tnearest_tu_end\tnearest_tu_three_prime_delta_bp\tnearest_tu_five_prime_delta_bp\tnearest_tu_score1\tnearest_tu_score2\tgene_hint"
    )?;

    for candidate in candidates {
        let nearest = candidate.nearest_tu.as_ref();
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            candidate_report_id(candidate),
            candidate.contig,
            candidate.strand.as_char(),
            candidate.interval.start.get(),
            candidate.interval.end.get(),
            candidate.three_prime_coord,
            candidate.five_prime_coord,
            candidate.family_support,
            candidate.mode_support,
            candidate.mode_fraction,
            candidate_reason(candidate),
            nearest.map(|tu| tu.id.as_str()).unwrap_or("."),
            nearest
                .map(|tu| tu.interval.start.get().to_string())
                .unwrap_or_else(|| ".".to_owned()),
            nearest
                .map(|tu| tu.interval.end.get().to_string())
                .unwrap_or_else(|| ".".to_owned()),
            nearest
                .map(|tu| tu.three_prime_delta_bp.to_string())
                .unwrap_or_else(|| ".".to_owned()),
            nearest
                .map(|tu| tu.five_prime_delta_bp.to_string())
                .unwrap_or_else(|| ".".to_owned()),
            nearest
                .map(|tu| format!("{:.6}", tu.score1))
                .unwrap_or_else(|| ".".to_owned()),
            nearest
                .map(|tu| format!("{:.6}", tu.score2))
                .unwrap_or_else(|| ".".to_owned()),
            candidate.gene_hint,
        )?;
    }

    Ok(())
}

fn write_bed_report(path: &Path, candidates: &[CandidateMode]) -> Result<()> {
    let records: Vec<Bed6Record> = candidates
        .iter()
        .map(|candidate| Bed6Record {
            chrom: candidate.contig.clone(),
            start: candidate.interval.start,
            end: candidate.interval.end,
            name: format!(
                "{}|{}|{}",
                candidate_report_id(candidate),
                candidate_reason(candidate),
                candidate.gene_hint
            ),
            score: candidate.mode_support.min(u32::MAX as usize) as u32,
            strand: candidate.strand,
            extra_fields: Vec::new(),
        })
        .collect();
    crate::io::bed::write_bed6(path, records.iter())
        .with_context(|| format!("failed to write BED report {:?}", path))?;
    Ok(())
}
