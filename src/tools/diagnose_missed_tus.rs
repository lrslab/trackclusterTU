use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;

use crate::tools::output_transaction::OutputTransaction;
use crate::tools::read_rejections::{
    emit_read_input_summary, write_read_rejections, ReadRejection,
};
use crate::tu::ReadRecord;

mod detection;
mod parsing;
mod rescue;
mod writing;

use detection::{detect_candidate_modes, CandidateMode, DetectionParams};
use parsing::{ExistingMembershipRow, ExistingTu, GeneRecord};

#[derive(Parser, Debug)]
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

    /// Read records excluded during BED parsing, validation, deduplication, or filtering.
    ///
    /// Defaults to `diagnose_read_rejections.tsv` beside `--out-tsv` and is written even
    /// when no reads are rejected.
    #[arg(long = "out-read-rejections")]
    out_read_rejections: Option<PathBuf>,

    /// Stop when any individual read is rejected instead of recording and continuing.
    #[arg(long = "strict-read-errors")]
    strict_read_errors: bool,

    /// Strand-aware 3 prime window used to form end families and TU matches (bp).
    #[arg(long, default_value_t = 12)]
    three_prime_window_bp: u32,

    /// Maximum 3 prime coordinate diameter within one end family (bp).
    ///
    /// Defaults to `--three-prime-window-bp`, which prevents single-linkage chains.
    #[arg(long)]
    max_three_prime_family_diameter_bp: Option<u32>,

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

#[derive(Parser, Debug)]
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

    /// Optional TSV report mapping candidate IDs to rescued TU IDs and final hard counts.
    #[arg(long = "out-candidates-tsv")]
    out_candidates_tsv: Option<PathBuf>,

    /// Optional BED6 track whose name field contains both candidate and rescued TU IDs.
    #[arg(long = "out-candidates-bed")]
    out_candidates_bed: Option<PathBuf>,

    /// Read records excluded during BED parsing, validation, deduplication, or filtering.
    ///
    /// Defaults to `rescue_read_rejections.tsv` beside `--out-tu` and is written even when
    /// no reads are rejected.
    #[arg(long = "out-read-rejections")]
    out_read_rejections: Option<PathBuf>,

    /// Stop when any individual read is rejected instead of recording and continuing.
    #[arg(long = "strict-read-errors")]
    strict_read_errors: bool,

    /// Prefix used when naming rescued TU IDs.
    #[arg(long, default_value = "RESC")]
    rescue_prefix: String,

    /// Strand-aware 3 prime window used to form end families and TU matches (bp).
    #[arg(long, default_value_t = 12)]
    three_prime_window_bp: u32,

    /// Maximum 3 prime coordinate diameter within one end family (bp).
    ///
    /// Defaults to `--three-prime-window-bp`, which prevents single-linkage chains.
    #[arg(long)]
    max_three_prime_family_diameter_bp: Option<u32>,

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

/// Validated application configuration for missed-TU diagnosis.
#[derive(Clone, Debug)]
struct DiagnoseConfig {
    input: PathBuf,
    existing_tu: PathBuf,
    annotation_bed: Option<PathBuf>,
    out_tsv: PathBuf,
    out_bed: Option<PathBuf>,
    out_read_rejections: PathBuf,
    strict_read_errors: bool,
    detection: DetectionParams,
    min_read_len: Option<u32>,
}

impl TryFrom<DiagnoseMissedTusCli> for DiagnoseConfig {
    type Error = anyhow::Error;

    fn try_from(cli: DiagnoseMissedTusCli) -> Result<Self> {
        let out_read_rejections = cli
            .out_read_rejections
            .clone()
            .unwrap_or_else(|| default_sidecar_path(&cli.out_tsv, "diagnose_read_rejections.tsv"));
        let detection = DetectionParams::try_new(
            cli.three_prime_window_bp,
            cli.max_three_prime_family_diameter_bp,
            cli.five_prime_window_bp,
            cli.min_family_support,
            cli.min_mode_support,
            cli.min_mode_fraction,
            cli.max_candidates_per_family,
        )?;
        validate_input_file(&cli.input, "input BED6")?;
        validate_input_file(&cli.existing_tu, "existing TU BED6")?;
        if let Some(path) = &cli.annotation_bed {
            validate_input_file(path, "annotation BED6")?;
        }
        validate_output_file(&cli.out_tsv, "output TSV")?;
        if let Some(path) = &cli.out_bed {
            validate_output_file(path, "output BED6")?;
        }
        validate_output_file(&out_read_rejections, "read rejection TSV")?;

        Ok(Self {
            input: cli.input,
            existing_tu: cli.existing_tu,
            annotation_bed: cli.annotation_bed,
            out_tsv: cli.out_tsv,
            out_bed: cli.out_bed,
            out_read_rejections,
            strict_read_errors: cli.strict_read_errors,
            detection,
            min_read_len: cli.min_read_len,
        })
    }
}

impl DiagnoseConfig {
    fn input_paths(&self) -> Vec<&Path> {
        [
            Some(self.input.as_path()),
            Some(self.existing_tu.as_path()),
            self.annotation_bed.as_deref(),
        ]
        .into_iter()
        .flatten()
        .collect()
    }

    fn output_paths(&self) -> Vec<&Path> {
        [
            Some(self.out_tsv.as_path()),
            self.out_bed.as_deref(),
            Some(self.out_read_rejections.as_path()),
        ]
        .into_iter()
        .flatten()
        .collect()
    }

    fn with_staged_outputs(&self, transaction: &OutputTransaction) -> Result<Self> {
        let mut staged = self.clone();
        staged.out_tsv = transaction.staged_path(&self.out_tsv)?;
        staged.out_bed = self
            .out_bed
            .as_deref()
            .map(|path| transaction.staged_path(path))
            .transpose()?;
        staged.out_read_rejections = transaction.staged_path(&self.out_read_rejections)?;
        Ok(staged)
    }
}

/// Validated application configuration for rescue and its complete output set.
#[derive(Clone, Debug)]
struct RescueConfig {
    input: PathBuf,
    existing_tu: PathBuf,
    existing_membership: PathBuf,
    annotation_bed: Option<PathBuf>,
    out_tu: PathBuf,
    out_membership: PathBuf,
    out_tu_count: Option<PathBuf>,
    out_candidates_tsv: Option<PathBuf>,
    out_candidates_bed: Option<PathBuf>,
    out_read_rejections: PathBuf,
    strict_read_errors: bool,
    rescue_prefix: String,
    detection: DetectionParams,
    min_read_len: Option<u32>,
}

impl TryFrom<RescueMissedTusCli> for RescueConfig {
    type Error = anyhow::Error;

    fn try_from(cli: RescueMissedTusCli) -> Result<Self> {
        let out_read_rejections = cli
            .out_read_rejections
            .clone()
            .unwrap_or_else(|| default_sidecar_path(&cli.out_tu, "rescue_read_rejections.tsv"));
        validate_rescue_prefix(&cli.rescue_prefix)?;
        let detection = DetectionParams::try_new(
            cli.three_prime_window_bp,
            cli.max_three_prime_family_diameter_bp,
            cli.five_prime_window_bp,
            cli.min_family_support,
            cli.min_mode_support,
            cli.min_mode_fraction,
            cli.max_candidates_per_family,
        )?;
        validate_input_file(&cli.input, "input BED6")?;
        validate_input_file(&cli.existing_tu, "existing TU BED6")?;
        validate_input_file(&cli.existing_membership, "existing membership TSV")?;
        if let Some(path) = &cli.annotation_bed {
            validate_input_file(path, "annotation BED6")?;
        }
        validate_output_file(&cli.out_tu, "rescued TU BED6")?;
        validate_output_file(&cli.out_membership, "rescued membership TSV")?;
        for (path, label) in [
            (cli.out_tu_count.as_ref(), "rescued TU count CSV"),
            (cli.out_candidates_tsv.as_ref(), "candidate TSV"),
            (cli.out_candidates_bed.as_ref(), "candidate BED6"),
        ] {
            if let Some(path) = path {
                validate_output_file(path, label)?;
            }
        }
        validate_output_file(&out_read_rejections, "read rejection TSV")?;

        Ok(Self {
            input: cli.input,
            existing_tu: cli.existing_tu,
            existing_membership: cli.existing_membership,
            annotation_bed: cli.annotation_bed,
            out_tu: cli.out_tu,
            out_membership: cli.out_membership,
            out_tu_count: cli.out_tu_count,
            out_candidates_tsv: cli.out_candidates_tsv,
            out_candidates_bed: cli.out_candidates_bed,
            out_read_rejections,
            strict_read_errors: cli.strict_read_errors,
            rescue_prefix: cli.rescue_prefix,
            detection,
            min_read_len: cli.min_read_len,
        })
    }
}

impl RescueConfig {
    fn input_paths(&self) -> Vec<&Path> {
        [
            Some(self.input.as_path()),
            Some(self.existing_tu.as_path()),
            Some(self.existing_membership.as_path()),
            self.annotation_bed.as_deref(),
        ]
        .into_iter()
        .flatten()
        .collect()
    }

    fn output_paths(&self) -> Vec<&Path> {
        [
            Some(self.out_tu.as_path()),
            Some(self.out_membership.as_path()),
            self.out_tu_count.as_deref(),
            self.out_candidates_tsv.as_deref(),
            self.out_candidates_bed.as_deref(),
            Some(self.out_read_rejections.as_path()),
        ]
        .into_iter()
        .flatten()
        .collect()
    }

    fn with_staged_outputs(&self, transaction: &OutputTransaction) -> Result<Self> {
        let mut staged = self.clone();
        staged.out_tu = transaction.staged_path(&self.out_tu)?;
        staged.out_membership = transaction.staged_path(&self.out_membership)?;
        staged.out_tu_count = staged_path(transaction, self.out_tu_count.as_deref())?;
        staged.out_candidates_tsv = staged_path(transaction, self.out_candidates_tsv.as_deref())?;
        staged.out_candidates_bed = staged_path(transaction, self.out_candidates_bed.as_deref())?;
        staged.out_read_rejections = transaction.staged_path(&self.out_read_rejections)?;
        Ok(staged)
    }
}

struct DiagnosisInputs {
    reads: Vec<ReadRecord>,
    existing_tus: Vec<ExistingTu>,
    genes: Vec<GeneRecord>,
    read_rejections: Vec<ReadRejection>,
    total_read_records: u64,
}

impl DiagnosisInputs {
    fn load(config: &DiagnoseConfig) -> Result<Self> {
        let (read_load, existing_tus, genes) = parsing::load_detection_inputs(
            &config.input,
            &config.existing_tu,
            config.annotation_bed.as_deref(),
            config.min_read_len,
        )?;
        enforce_strict_read_errors(config.strict_read_errors, &read_load.outcome.rejections)?;
        let total_read_records = read_load.outcome.total_records;
        let read_rejections = read_load.outcome.rejections;
        let reads = read_load
            .outcome
            .reads
            .into_iter()
            .map(|located| located.read)
            .collect();
        Ok(Self {
            reads,
            existing_tus,
            genes,
            read_rejections,
            total_read_records,
        })
    }
}

struct RescueInputs {
    reads: Vec<ReadRecord>,
    existing_tus: Vec<ExistingTu>,
    genes: Vec<GeneRecord>,
    read_to_existing: Vec<Option<usize>>,
    existing_memberships: Vec<Option<ExistingMembershipRow>>,
    membership_schema: rescue::MembershipSchema,
    membership_metadata: Vec<String>,
    read_rejections: Vec<ReadRejection>,
    total_read_records: u64,
}

impl RescueInputs {
    fn load(config: &RescueConfig) -> Result<Self> {
        let read_load = parsing::read_reads_bed6(&config.input, config.min_read_len)?;
        enforce_strict_read_errors(config.strict_read_errors, &read_load.outcome.rejections)?;
        let input_read_ids = read_load.input_read_ids;
        let total_read_records = read_load.outcome.total_records;
        let read_rejections = read_load.outcome.rejections;
        let reads: Vec<ReadRecord> = read_load
            .outcome
            .reads
            .into_iter()
            .map(|located| located.read)
            .collect();
        let existing_tus = parsing::read_existing_tus_bed6(&config.existing_tu)?;
        let genes = match config.annotation_bed.as_deref() {
            Some(path) => parsing::read_gene_bed6(path)?,
            None => Vec::new(),
        };
        let membership_rows: Vec<ExistingMembershipRow> =
            parsing::read_existing_membership(&config.existing_membership)?;
        let membership_metadata =
            crate::tools::membership::read_membership_metadata(&config.existing_membership)?;
        let validated_memberships = rescue::validate_existing_memberships(
            membership_rows,
            &membership_metadata,
            &input_read_ids,
            &reads,
            &existing_tus,
        )
        .map_err(|error| anyhow::anyhow!("{:?}:{error}", config.existing_membership))?;

        Ok(Self {
            reads,
            existing_tus,
            genes,
            read_to_existing: validated_memberships.read_to_existing,
            existing_memberships: validated_memberships.rows_by_read,
            membership_schema: validated_memberships.schema,
            membership_metadata,
            read_rejections,
            total_read_records,
        })
    }
}

pub(crate) fn run_from_args<I, T>(args: I) -> Result<()>
where
    I: IntoIterator<Item = T>,
    T: Into<std::ffi::OsString> + Clone,
{
    let cli = DiagnoseMissedTusCli::parse_from(args);
    let config = DiagnoseConfig::try_from(cli)?;
    // Validate and parse every input before the transaction creates any output directory.
    let inputs = DiagnosisInputs::load(&config)?;
    let total_read_records = inputs.total_read_records;
    let retained_reads = inputs.reads.len();
    let transaction = OutputTransaction::new(config.input_paths(), config.output_paths())?;
    let staged = config.with_staged_outputs(&transaction)?;
    let read_rejections = execute_diagnosis(&staged, inputs)?;
    transaction.commit()?;
    emit_read_input_summary(total_read_records, retained_reads, &read_rejections);

    println!("out_tsv={}", config.out_tsv.display());
    if let Some(out_bed) = config.out_bed.as_ref() {
        println!("out_bed={}", out_bed.display());
    }
    println!(
        "out_read_rejections={}",
        config.out_read_rejections.display()
    );
    Ok(())
}

pub(crate) fn run_rescue_from_args<I, T>(args: I) -> Result<()>
where
    I: IntoIterator<Item = T>,
    T: Into<std::ffi::OsString> + Clone,
{
    let cli = RescueMissedTusCli::parse_from(args);
    let config = RescueConfig::try_from(cli)?;
    // Structural inputs and membership foreign keys validate before output setup; recoverable
    // read-record problems have already been quarantined in `inputs`.
    let inputs = RescueInputs::load(&config)?;
    let total_read_records = inputs.total_read_records;
    let retained_reads = inputs.reads.len();
    let transaction = OutputTransaction::new(config.input_paths(), config.output_paths())?;
    let staged = config.with_staged_outputs(&transaction)?;
    let read_rejections = execute_rescue(&staged, inputs)?;
    transaction.commit()?;
    emit_read_input_summary(total_read_records, retained_reads, &read_rejections);

    println!("out_tu={}", config.out_tu.display());
    println!("out_membership={}", config.out_membership.display());
    if let Some(out_tu_count) = config.out_tu_count.as_ref() {
        println!("out_tu_count={}", out_tu_count.display());
    }
    println!(
        "out_read_rejections={}",
        config.out_read_rejections.display()
    );
    Ok(())
}

fn execute_diagnosis(
    config: &DiagnoseConfig,
    inputs: DiagnosisInputs,
) -> Result<Vec<ReadRejection>> {
    let mut candidates = detect_candidate_modes(
        &inputs.reads,
        &inputs.existing_tus,
        &inputs.genes,
        config.detection,
    );
    writing::assign_candidate_ids(&mut candidates, "MISS");
    writing::write_tsv_report(&config.out_tsv, &candidates)?;
    if let Some(out_bed) = config.out_bed.as_ref() {
        writing::write_bed_report(out_bed, &candidates)?;
    }
    write_read_rejections(&config.out_read_rejections, &inputs.read_rejections)?;
    println!("candidate_count={}", candidates.len());
    Ok(inputs.read_rejections)
}

fn execute_rescue(config: &RescueConfig, inputs: RescueInputs) -> Result<Vec<ReadRejection>> {
    let candidate_modes: Vec<CandidateMode> = detect_candidate_modes(
        &inputs.reads,
        &inputs.existing_tus,
        &inputs.genes,
        config.detection,
    );
    let mut plan = rescue::build_rescue_plan(
        &inputs.reads,
        &inputs.existing_tus,
        candidate_modes,
        inputs.read_to_existing,
        config.detection,
        &config.rescue_prefix,
    )?;

    writing::write_rescued_tu_bed(&config.out_tu, &plan.final_tus)?;
    writing::write_rescued_membership(
        &config.out_membership,
        &inputs.reads,
        &inputs.existing_memberships,
        inputs.membership_schema,
        &inputs.membership_metadata,
        &inputs.existing_tus,
        &plan,
    )?;
    if let Some(out_tu_count) = config.out_tu_count.as_ref() {
        writing::write_rescued_counts_csv(out_tu_count, &plan.final_tus)?;
    }

    // Sorting for report IDs happens only after all index-aligned rescue outputs are complete.
    // Candidate and rescued-TU IDs move together so report mappings remain explicit.
    writing::assign_rescued_candidate_ids(&mut plan.candidates, &mut plan.rescue_ids, "MISS")?;
    if let Some(out_tsv) = config.out_candidates_tsv.as_ref() {
        writing::write_rescued_tsv_report(
            out_tsv,
            &plan.candidates,
            &plan.rescue_ids,
            &plan.final_tus,
        )?;
    }
    if let Some(out_bed) = config.out_candidates_bed.as_ref() {
        writing::write_rescued_bed_report(out_bed, &plan.candidates, &plan.rescue_ids)?;
    }
    write_read_rejections(&config.out_read_rejections, &inputs.read_rejections)?;
    println!("rescued_candidate_count={}", plan.candidates.len());
    Ok(inputs.read_rejections)
}

fn enforce_strict_read_errors(strict: bool, rejections: &[ReadRejection]) -> Result<()> {
    if strict && !rejections.is_empty() {
        anyhow::bail!(
            "--strict-read-errors rejected {} read record(s); first rejection: {}",
            rejections.len(),
            rejections[0].context()
        );
    }
    Ok(())
}

fn validate_rescue_prefix(prefix: &str) -> Result<()> {
    if prefix.trim().is_empty() {
        anyhow::bail!("--rescue-prefix must not be empty");
    }
    if prefix.contains(['\t', '\n', '\r', ',']) {
        anyhow::bail!(
            "--rescue-prefix must not contain tabs, newlines, carriage returns, or commas"
        );
    }
    Ok(())
}

fn validate_input_file(path: &Path, label: &str) -> Result<()> {
    let metadata = std::fs::metadata(path)
        .with_context(|| format!("{label} not found or unreadable: {}", path.display()))?;
    if !metadata.is_file() {
        anyhow::bail!("{label} is not a regular file: {}", path.display());
    }
    Ok(())
}

fn validate_output_file(path: &Path, label: &str) -> Result<()> {
    if path.as_os_str().is_empty() || path.file_name().is_none() {
        anyhow::bail!("{label} must name a file");
    }
    if path.is_dir() {
        anyhow::bail!("{label} is an existing directory: {}", path.display());
    }
    Ok(())
}

fn default_sidecar_path(primary_output: &Path, filename: &str) -> PathBuf {
    primary_output
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .join(filename)
}

fn staged_path(transaction: &OutputTransaction, path: Option<&Path>) -> Result<Option<PathBuf>> {
    path.map(|path| transaction.staged_path(path)).transpose()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn invalid_cli_values_fail_during_config_translation() {
        let temp = std::env::temp_dir();
        let cli = RescueMissedTusCli {
            input: temp.join("not-read-because-prefix-fails.bed"),
            existing_tu: temp.join("not-read-because-prefix-fails.tus.bed"),
            existing_membership: temp.join("not-read-because-prefix-fails.tsv"),
            annotation_bed: None,
            out_tu: temp.join("not-created.tus.bed"),
            out_membership: temp.join("not-created.membership.tsv"),
            out_tu_count: None,
            out_candidates_tsv: None,
            out_candidates_bed: None,
            out_read_rejections: None,
            strict_read_errors: false,
            rescue_prefix: "bad\tprefix".to_owned(),
            three_prime_window_bp: 12,
            max_three_prime_family_diameter_bp: None,
            five_prime_window_bp: 10,
            min_family_support: 20,
            min_mode_support: 20,
            min_mode_fraction: 0.02,
            max_candidates_per_family: 3,
            min_read_len: None,
        };
        let error = RescueConfig::try_from(cli).unwrap_err().to_string();
        assert!(error.contains("--rescue-prefix"));
    }
}
