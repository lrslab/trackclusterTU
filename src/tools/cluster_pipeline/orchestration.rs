use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::time::Instant;

use clap::Parser;
use rayon::prelude::*;
use thiserror::Error;

use super::annotation::{
    classify_tus, format_name2, merge_overlapping_intervals, read_annotation_bed6,
    tu_gene_contexts, write_tu_gene_relationships, write_tu_gff3, write_tu_semantics,
    GeneOverlapPolicy, GeneRecord,
};
use super::config::{
    AnnotationRequest, ClusterConfig, ClusterInput, ClusterOptions, ClusterOutputPaths,
    ClusterRequest, ConfigError, InputFormat, RecountConfig, RecountOutputPaths, RecountRequest,
    TuIdStyle,
};
use super::counting::{
    build_gene_multi_counts, build_gene_total_counts, build_multi_counts_from_membership,
    total_counts_from_multi, write_count_metrics_csv, write_entity_metric_long_counts,
    write_entity_metric_matrix, write_gene_count_metrics_csv, CountMetrics, GENE_COUNT_SEMANTICS,
};
use super::input::{
    cmp_read, filter_reads, load_manifest_reads, read_reads, resolve_manifest_input_format,
    resolve_single_input_format, ResolvedInputFormat,
};
use super::output::{
    create_membership_writer, hard_assignment_tu_index, prepare_outputs, write_membership_row,
    write_pooled_reads_bed6, write_tu_endpoint_stats, write_tu_id_map, PreparedOutputs,
};

use crate::model::{Bed12Attrs, Interval, Transcript};
use crate::tools::output_transaction::OutputTransaction;
use crate::tools::read_rejections::{
    emit_read_input_summary, write_read_rejections, ReadLoadOutcome,
};
use crate::tu::multi::{
    parse_manifest, MultiSampleCounter, MultiSampleCounts, SampleManifestRecord,
};
use crate::tu::{AssignmentStatus, TuClusteringOptions};

#[derive(Parser, Debug, Clone)]
#[command(
    name = "trackclustertu cluster",
    version,
    about = "Cluster bacterial directRNA reads into transcript units (TUs)",
    group(
        clap::ArgGroup::new("input_mode")
            .required(true)
            .args(["input", "manifest"])
    )
)]
struct ClusterCli {
    /// Single input reads file path.
    #[arg(long = "in", conflicts_with = "manifest")]
    input: Option<PathBuf>,

    /// Multi-sample manifest TSV with columns: sample, reads, [group], [evidence].
    #[arg(long, conflicts_with = "input")]
    manifest: Option<PathBuf>,

    /// Input format. Defaults to auto-detecting from file extensions.
    #[arg(long, value_enum, default_value_t = InputFormat::Auto)]
    format: InputFormat,

    /// Output directory. Missing output paths default to files inside this directory.
    #[arg(long)]
    out_dir: Option<PathBuf>,

    /// Span Jaccard threshold (overlap / union).
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

    /// Output shared TU BED6.
    #[arg(long)]
    out_tu: Option<PathBuf>,

    /// Output pooled membership TSV.
    #[arg(long)]
    out_membership: Option<PathBuf>,

    /// Optional pooled BED6 used for clustering after sample-tagging read IDs.
    #[arg(long, requires = "manifest")]
    out_pooled_reads: Option<PathBuf>,

    /// Read records excluded during parsing, validation, deduplication, or filtering.
    ///
    /// Defaults to `read_rejections.tsv` inside `--out-dir` and is written even when empty.
    #[arg(long)]
    out_read_rejections: Option<PathBuf>,

    /// Fail atomically if any read is rejected; no output set is published.
    #[arg(long)]
    strict_read_errors: bool,

    /// Endpoint consensus/support audit table. Defaults to `tu_endpoint_stats.tsv`.
    #[arg(long)]
    out_tu_endpoint_stats: Option<PathBuf>,

    /// TU identifier strategy. Stable coordinate-derived IDs are the default.
    #[arg(long, value_enum, default_value_t = TuIdStyle::Stable)]
    tu_id_style: TuIdStyle,

    /// Emitted/stable/sequential TU identifier mapping table.
    #[arg(long)]
    out_tu_id_map: Option<PathBuf>,

    /// Maximum overlap-over-longer difference marking the best and runner-up as ambiguous.
    #[arg(long, default_value_t = 0.02)]
    ambiguity_margin: f64,

    /// Split an ambiguous read equally across every candidate inside `--ambiguity-margin`.
    #[arg(long)]
    fractional_assignment: bool,

    /// TU count CSV with legacy `count` plus unique, evidence, total, and fractional columns.
    #[arg(long)]
    out_tu_count: Option<PathBuf>,

    /// Per-sample TU long table with explicit unique/evidence/total/fractional columns.
    #[arg(long, requires = "manifest")]
    out_tu_sample_count_long: Option<PathBuf>,

    /// Per-sample TU matrix with `<sample>.<count_semantics>` columns.
    #[arg(long, requires = "manifest")]
    out_tu_sample_count_matrix: Option<PathBuf>,

    /// Per-group TU matrix with `<group>.<count_semantics>` columns.
    ///
    /// This file is only written when the manifest has a non-empty `group` column.
    #[arg(long, requires = "manifest")]
    out_tu_group_count_matrix: Option<PathBuf>,

    /// Minimum pre-assignment clustering-family support required to emit a TU.
    #[arg(long)]
    min_tu_count: Option<u64>,

    /// Optional annotation BED6 (e.g., genes) used to anchor TU results.
    #[arg(long)]
    annotation_bed: Option<PathBuf>,

    /// Minimum overlap in base pairs for a TU/gene relationship.
    #[arg(long, default_value_t = 1, requires = "annotation_bed")]
    gene_min_overlap_bp: u32,

    /// Minimum fraction of the TU covered by a gene relationship.
    #[arg(long, default_value_t = 0.0, requires = "annotation_bed")]
    gene_min_tu_fraction: f64,

    /// Minimum fraction of the gene covered by a TU relationship.
    #[arg(long, default_value_t = 0.0, requires = "annotation_bed")]
    gene_min_gene_fraction: f64,

    /// Optional TU-to-gene overlap TSV (one line per TU×gene overlap).
    #[arg(long, requires = "annotation_bed")]
    out_tu_gene: Option<PathBuf>,

    /// Versioned TU semantics table with ordered gene-context signatures and classifications.
    #[arg(long, requires = "annotation_bed")]
    out_tu_semantics: Option<PathBuf>,

    /// Standard GFF3 transcript-unit features and TU-to-gene relationship features.
    #[arg(long, requires = "annotation_bed")]
    out_tu_gff3: Option<PathBuf>,

    /// Optional TU BED12 output anchored to qualifying same-strand genes.
    ///
    /// Blocks are clipped to same-strand gene intersections. When none qualify, the whole TU
    /// span is emitted as one block.
    ///
    /// Adds extra columns:
    /// - name2: comma-separated subreads, with `|<read_count>` suffix (TrackCluster-style)
    /// - gene_list: comma-separated qualifying same-strand genes (or '.')
    #[arg(long, requires = "annotation_bed")]
    out_tu_bed12: Option<PathBuf>,

    /// Gene count CSV with legacy `count` plus explicit assignment-semantics columns.
    #[arg(long, requires = "annotation_bed")]
    out_gene_count: Option<PathBuf>,

    /// Per-sample gene matrix with `<sample>.<count_semantics>` columns.
    #[arg(long, requires_all = ["annotation_bed", "manifest"])]
    out_gene_sample_count_matrix: Option<PathBuf>,

    /// Per-group gene matrix with `<group>.<count_semantics>` columns.
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

    /// TU count CSV with legacy `count` plus unique, evidence, total, and fractional columns.
    #[arg(long)]
    out_tu_count: Option<PathBuf>,

    /// Per-sample TU long table with explicit unique/evidence/total/fractional columns.
    #[arg(long)]
    out_tu_sample_count_long: Option<PathBuf>,

    /// Per-sample TU matrix with `<sample>.<count_semantics>` columns.
    #[arg(long)]
    out_tu_sample_count_matrix: Option<PathBuf>,

    /// Per-group TU matrix with `<group>.<count_semantics>` columns.
    ///
    /// This file is only written when the manifest has a non-empty `group` column.
    #[arg(long)]
    out_tu_group_count_matrix: Option<PathBuf>,

    /// Minimum total hard-assignment count required to retain a TU row.
    ///
    /// Fractional-only contributions do not satisfy this threshold.
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
struct ExecutionConfig {
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
    out_read_rejections: Option<PathBuf>,
    strict_read_errors: bool,
    out_tu_endpoint_stats: Option<PathBuf>,
    tu_id_style: TuIdStyle,
    out_tu_id_map: Option<PathBuf>,
    ambiguity_margin: f64,
    fractional_assignment: bool,
    out_tu_count: Option<PathBuf>,
    out_tu_sample_count_long: Option<PathBuf>,
    out_tu_sample_count_matrix: Option<PathBuf>,
    out_tu_group_count_matrix: Option<PathBuf>,
    min_tu_count: Option<u64>,
    annotation_bed: Option<PathBuf>,
    gene_min_overlap_bp: u32,
    gene_min_tu_fraction: f64,
    gene_min_gene_fraction: f64,
    out_tu_gene: Option<PathBuf>,
    out_tu_semantics: Option<PathBuf>,
    out_tu_gff3: Option<PathBuf>,
    out_tu_bed12: Option<PathBuf>,
    out_gene_count: Option<PathBuf>,
    out_gene_sample_count_matrix: Option<PathBuf>,
    out_gene_group_count_matrix: Option<PathBuf>,
    threads: Option<usize>,
    timings: bool,
}

impl TryFrom<ClusterCli> for ClusterConfig {
    type Error = ConfigError;

    fn try_from(cli: ClusterCli) -> Result<Self, Self::Error> {
        let input = match (cli.input, cli.manifest) {
            (Some(path), None) => ClusterInput::Single(path),
            (None, Some(path)) => ClusterInput::Manifest(path),
            (None, None) => return Err(ConfigError::MissingInput),
            (Some(_), Some(_)) => return Err(ConfigError::ConflictingInput),
        };
        let annotation = cli.annotation_bed.map(|bed| AnnotationRequest {
            bed,
            min_overlap_bp: cli.gene_min_overlap_bp,
            min_tu_fraction: cli.gene_min_tu_fraction,
            min_gene_fraction: cli.gene_min_gene_fraction,
        });
        ClusterConfig::validate(ClusterRequest {
            input,
            format: cli.format,
            outputs: ClusterOutputPaths {
                out_dir: cli.out_dir,
                tu: cli.out_tu,
                membership: cli.out_membership,
                pooled_reads: cli.out_pooled_reads,
                read_rejections: cli.out_read_rejections,
                endpoint_stats: cli.out_tu_endpoint_stats,
                id_map: cli.out_tu_id_map,
                tu_count: cli.out_tu_count,
                tu_sample_long: cli.out_tu_sample_count_long,
                tu_sample_matrix: cli.out_tu_sample_count_matrix,
                tu_group_matrix: cli.out_tu_group_count_matrix,
                tu_gene: cli.out_tu_gene,
                tu_semantics: cli.out_tu_semantics,
                tu_gff3: cli.out_tu_gff3,
                tu_bed12: cli.out_tu_bed12,
                gene_count: cli.out_gene_count,
                gene_sample_matrix: cli.out_gene_sample_count_matrix,
                gene_group_matrix: cli.out_gene_group_count_matrix,
            },
            options: ClusterOptions {
                span_jaccard_threshold: cli.score1_threshold,
                overlap_over_longer_threshold: cli.score2_threshold,
                three_prime_tolerance_bp: cli.three_prime_tolerance_bp,
                max_five_prime_delta_bp: cli.max_five_prime_delta_bp,
                attach_contained_reads: !cli.skip_score2_attachment,
                ambiguity_margin: cli.ambiguity_margin,
                fractional_assignment: cli.fractional_assignment,
                tu_id_style: cli.tu_id_style,
                strict_read_errors: cli.strict_read_errors,
            },
            min_read_len: cli.min_read_len,
            min_tu_count: cli.min_tu_count,
            annotation,
            threads: cli.threads,
            timings: cli.timings,
        })
    }
}

impl TryFrom<RecountCli> for RecountConfig {
    type Error = ConfigError;

    fn try_from(cli: RecountCli) -> Result<Self, Self::Error> {
        RecountConfig::validate(RecountRequest {
            manifest: cli.manifest,
            membership: cli.pooled_membership,
            outputs: RecountOutputPaths {
                out_dir: cli.out_dir,
                tu_count: cli.out_tu_count,
                tu_sample_long: cli.out_tu_sample_count_long,
                tu_sample_matrix: cli.out_tu_sample_count_matrix,
                tu_group_matrix: cli.out_tu_group_count_matrix,
            },
            min_tu_count: cli.min_tu_count,
            threads: cli.threads,
            timings: cli.timings,
        })
    }
}

impl From<ClusterConfig> for ExecutionConfig {
    fn from(config: ClusterConfig) -> Self {
        let (input, manifest) = match config.input {
            ClusterInput::Single(path) => (Some(path), None),
            ClusterInput::Manifest(path) => (None, Some(path)),
        };
        let (annotation_bed, min_overlap, min_tu_fraction, min_gene_fraction) = config
            .annotation
            .map(|annotation| {
                (
                    Some(annotation.bed),
                    annotation.min_overlap_bp,
                    annotation.min_tu_fraction,
                    annotation.min_gene_fraction,
                )
            })
            .unwrap_or((None, 1, 0.0, 0.0));
        let outputs = config.outputs;
        Self {
            input,
            manifest,
            pooled_membership: None,
            format: config.format,
            out_dir: outputs.out_dir,
            score1_threshold: config.options.span_jaccard_threshold,
            score2_threshold: config.options.overlap_over_longer_threshold,
            three_prime_tolerance_bp: config.options.three_prime_tolerance_bp,
            max_five_prime_delta_bp: config.options.max_five_prime_delta_bp,
            skip_score2_attachment: !config.options.attach_contained_reads,
            min_read_len: config.min_read_len,
            out_tu: outputs.tu,
            out_membership: outputs.membership,
            out_pooled_reads: outputs.pooled_reads,
            out_read_rejections: outputs.read_rejections,
            strict_read_errors: config.options.strict_read_errors,
            out_tu_endpoint_stats: outputs.endpoint_stats,
            tu_id_style: config.options.tu_id_style,
            out_tu_id_map: outputs.id_map,
            ambiguity_margin: config.options.ambiguity_margin,
            fractional_assignment: config.options.fractional_assignment,
            out_tu_count: outputs.tu_count,
            out_tu_sample_count_long: outputs.tu_sample_long,
            out_tu_sample_count_matrix: outputs.tu_sample_matrix,
            out_tu_group_count_matrix: outputs.tu_group_matrix,
            min_tu_count: config.min_tu_count,
            annotation_bed,
            gene_min_overlap_bp: min_overlap,
            gene_min_tu_fraction: min_tu_fraction,
            gene_min_gene_fraction: min_gene_fraction,
            out_tu_gene: outputs.tu_gene,
            out_tu_semantics: outputs.tu_semantics,
            out_tu_gff3: outputs.tu_gff3,
            out_tu_bed12: outputs.tu_bed12,
            out_gene_count: outputs.gene_count,
            out_gene_sample_count_matrix: outputs.gene_sample_matrix,
            out_gene_group_count_matrix: outputs.gene_group_matrix,
            threads: Some(config.threads),
            timings: config.timings,
        }
    }
}

impl From<RecountConfig> for ExecutionConfig {
    fn from(config: RecountConfig) -> Self {
        let outputs = config.outputs;
        Self {
            input: None,
            manifest: Some(config.manifest),
            pooled_membership: Some(config.membership),
            format: InputFormat::Auto,
            out_dir: outputs.out_dir,
            score1_threshold: 0.95,
            score2_threshold: 0.80,
            three_prime_tolerance_bp: 12,
            max_five_prime_delta_bp: None,
            skip_score2_attachment: false,
            min_read_len: None,
            out_tu: None,
            out_membership: None,
            out_pooled_reads: None,
            out_read_rejections: None,
            strict_read_errors: false,
            out_tu_endpoint_stats: None,
            tu_id_style: TuIdStyle::Stable,
            out_tu_id_map: None,
            ambiguity_margin: 0.02,
            fractional_assignment: false,
            out_tu_count: outputs.tu_count,
            out_tu_sample_count_long: outputs.tu_sample_long,
            out_tu_sample_count_matrix: outputs.tu_sample_matrix,
            out_tu_group_count_matrix: outputs.tu_group_matrix,
            min_tu_count: config.min_tu_count,
            annotation_bed: None,
            gene_min_overlap_bp: 1,
            gene_min_tu_fraction: 0.0,
            gene_min_gene_fraction: 0.0,
            out_tu_gene: None,
            out_tu_semantics: None,
            out_tu_gff3: None,
            out_tu_bed12: None,
            out_gene_count: None,
            out_gene_sample_count_matrix: None,
            out_gene_group_count_matrix: None,
            threads: Some(config.threads),
            timings: config.timings,
        }
    }
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

fn default_output_dir(cli: &ExecutionConfig, mode: &InputMode) -> PathBuf {
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

fn apply_default_outputs(cli: &mut ExecutionConfig, mode: &InputMode) {
    let out_dir = default_output_dir(cli, mode);
    cli.out_dir = Some(out_dir.clone());

    match mode {
        InputMode::Single => {
            set_default_path(&mut cli.out_tu, &out_dir, "tus.bed");
            set_default_path(&mut cli.out_membership, &out_dir, "membership.tsv");
            set_default_path(
                &mut cli.out_tu_endpoint_stats,
                &out_dir,
                "tu_endpoint_stats.tsv",
            );
            set_default_path(&mut cli.out_tu_id_map, &out_dir, "tu_id_map.tsv");
            set_default_path(&mut cli.out_tu_count, &out_dir, "tu_count.csv");
            set_default_path(
                &mut cli.out_read_rejections,
                &out_dir,
                "read_rejections.tsv",
            );

            if cli.annotation_bed.is_some() {
                set_default_path(&mut cli.out_tu_gene, &out_dir, "tu_gene.tsv");
                set_default_path(&mut cli.out_tu_semantics, &out_dir, "tu_semantics.tsv");
                set_default_path(&mut cli.out_tu_gff3, &out_dir, "tus.gff3");
                set_default_path(&mut cli.out_tu_bed12, &out_dir, "tus.anchored.bed12");
                set_default_path(&mut cli.out_gene_count, &out_dir, "gene_count.csv");
            }
        }
        InputMode::ManifestCluster { .. } => {
            set_default_path(&mut cli.out_tu, &out_dir, "tus.bed");
            set_default_path(&mut cli.out_membership, &out_dir, "membership.tsv");
            set_default_path(&mut cli.out_pooled_reads, &out_dir, "pooled.bed");
            set_default_path(
                &mut cli.out_tu_endpoint_stats,
                &out_dir,
                "tu_endpoint_stats.tsv",
            );
            set_default_path(&mut cli.out_tu_id_map, &out_dir, "tu_id_map.tsv");
            set_default_path(&mut cli.out_tu_count, &out_dir, "tu_count.csv");
            set_default_path(
                &mut cli.out_read_rejections,
                &out_dir,
                "read_rejections.tsv",
            );
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
                set_default_path(&mut cli.out_tu_semantics, &out_dir, "tu_semantics.tsv");
                set_default_path(&mut cli.out_tu_gff3, &out_dir, "tus.gff3");
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

fn prepare_output_dirs(cli: &ExecutionConfig) -> anyhow::Result<()> {
    if let Some(out_dir) = cli.out_dir.as_ref() {
        fs::create_dir_all(out_dir)?;
    }

    let file_slots = [
        cli.out_tu.as_ref(),
        cli.out_membership.as_ref(),
        cli.out_pooled_reads.as_ref(),
        cli.out_read_rejections.as_ref(),
        cli.out_tu_endpoint_stats.as_ref(),
        cli.out_tu_id_map.as_ref(),
        cli.out_tu_count.as_ref(),
        cli.out_tu_sample_count_long.as_ref(),
        cli.out_tu_sample_count_matrix.as_ref(),
        cli.out_tu_group_count_matrix.as_ref(),
        cli.out_tu_gene.as_ref(),
        cli.out_tu_semantics.as_ref(),
        cli.out_tu_gff3.as_ref(),
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

fn command_input_paths(cli: &ExecutionConfig, mode: &InputMode) -> Vec<PathBuf> {
    let mut paths: Vec<PathBuf> = [
        cli.input.as_ref(),
        cli.manifest.as_ref(),
        cli.pooled_membership.as_ref(),
        cli.annotation_bed.as_ref(),
    ]
    .into_iter()
    .flatten()
    .cloned()
    .collect();

    match mode {
        InputMode::Single => {}
        InputMode::ManifestCluster { samples } => {
            paths.extend(samples.iter().map(|sample| sample.reads.clone()));
            paths.extend(samples.iter().filter_map(|sample| sample.evidence.clone()));
        }
        InputMode::ManifestCountOnly { samples, .. } => {
            paths.extend(samples.iter().map(|sample| sample.reads.clone()));
        }
    }
    paths
}

fn command_output_paths(cli: &ExecutionConfig) -> Vec<PathBuf> {
    [
        cli.out_tu.as_ref(),
        cli.out_membership.as_ref(),
        cli.out_pooled_reads.as_ref(),
        cli.out_read_rejections.as_ref(),
        cli.out_tu_endpoint_stats.as_ref(),
        cli.out_tu_id_map.as_ref(),
        cli.out_tu_count.as_ref(),
        cli.out_tu_sample_count_long.as_ref(),
        cli.out_tu_sample_count_matrix.as_ref(),
        cli.out_tu_group_count_matrix.as_ref(),
        cli.out_tu_gene.as_ref(),
        cli.out_tu_semantics.as_ref(),
        cli.out_tu_gff3.as_ref(),
        cli.out_tu_bed12.as_ref(),
        cli.out_gene_count.as_ref(),
        cli.out_gene_sample_count_matrix.as_ref(),
        cli.out_gene_group_count_matrix.as_ref(),
    ]
    .into_iter()
    .flatten()
    .cloned()
    .collect()
}

fn stage_output_path(
    slot: &mut Option<PathBuf>,
    transaction: &OutputTransaction,
) -> anyhow::Result<()> {
    if let Some(final_path) = slot.as_ref() {
        *slot = Some(transaction.staged_path(final_path)?);
    }
    Ok(())
}

fn staged_execution_config(
    cli: &ExecutionConfig,
    transaction: &OutputTransaction,
) -> anyhow::Result<ExecutionConfig> {
    let mut staged = cli.clone();
    stage_output_path(&mut staged.out_tu, transaction)?;
    stage_output_path(&mut staged.out_membership, transaction)?;
    stage_output_path(&mut staged.out_pooled_reads, transaction)?;
    stage_output_path(&mut staged.out_read_rejections, transaction)?;
    stage_output_path(&mut staged.out_tu_endpoint_stats, transaction)?;
    stage_output_path(&mut staged.out_tu_id_map, transaction)?;
    stage_output_path(&mut staged.out_tu_count, transaction)?;
    stage_output_path(&mut staged.out_tu_sample_count_long, transaction)?;
    stage_output_path(&mut staged.out_tu_sample_count_matrix, transaction)?;
    stage_output_path(&mut staged.out_tu_group_count_matrix, transaction)?;
    stage_output_path(&mut staged.out_tu_gene, transaction)?;
    stage_output_path(&mut staged.out_tu_semantics, transaction)?;
    stage_output_path(&mut staged.out_tu_gff3, transaction)?;
    stage_output_path(&mut staged.out_tu_bed12, transaction)?;
    stage_output_path(&mut staged.out_gene_count, transaction)?;
    stage_output_path(&mut staged.out_gene_sample_count_matrix, transaction)?;
    stage_output_path(&mut staged.out_gene_group_count_matrix, transaction)?;
    Ok(staged)
}

fn default_threads() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

fn resolve_mode(cli: &ExecutionConfig) -> anyhow::Result<InputMode> {
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

fn want_multi_sample_outputs(cli: &ExecutionConfig) -> bool {
    cli.out_tu_sample_count_long.is_some()
        || cli.out_tu_sample_count_matrix.is_some()
        || cli.out_tu_group_count_matrix.is_some()
}

fn want_gene_multi_sample_outputs(cli: &ExecutionConfig) -> bool {
    cli.out_gene_sample_count_matrix.is_some() || cli.out_gene_group_count_matrix.is_some()
}

fn want_any_multi_sample_outputs(cli: &ExecutionConfig) -> bool {
    want_multi_sample_outputs(cli) || want_gene_multi_sample_outputs(cli)
}

fn want_gene_outputs(cli: &ExecutionConfig) -> bool {
    cli.out_gene_count.is_some() || want_gene_multi_sample_outputs(cli)
}

fn want_annotation_outputs(cli: &ExecutionConfig) -> bool {
    cli.out_tu_gene.is_some()
        || cli.out_tu_semantics.is_some()
        || cli.out_tu_gff3.is_some()
        || cli.out_tu_bed12.is_some()
        || want_gene_outputs(cli)
}

fn gene_overlap_policy(cli: &ExecutionConfig) -> GeneOverlapPolicy {
    GeneOverlapPolicy {
        min_overlap_bp: cli.gene_min_overlap_bp,
        min_tu_fraction: cli.gene_min_tu_fraction,
        min_gene_fraction: cli.gene_min_gene_fraction,
    }
}

fn want_any_tu_count_outputs(cli: &ExecutionConfig) -> bool {
    cli.out_tu_count.is_some() || want_multi_sample_outputs(cli)
}

fn write_multi_sample_outputs(
    cli: &ExecutionConfig,
    counts: &MultiSampleCounts,
) -> anyhow::Result<()> {
    if let Some(path) = cli.out_tu_sample_count_long.as_ref() {
        write_entity_metric_long_counts(
            path,
            "tu_id",
            "sample",
            counts.tu_ids(),
            counts.samples(),
            counts.unique_counts(),
            counts.full_length_evidence_counts(),
            counts.counts(),
            counts.fractional_counts(),
        )?;
    }

    if let Some(path) = cli.out_tu_sample_count_matrix.as_ref() {
        write_entity_metric_matrix(
            path,
            None,
            "tu_id",
            counts.tu_ids(),
            counts.samples(),
            counts.unique_counts(),
            counts.full_length_evidence_counts(),
            counts.counts(),
            counts.fractional_counts(),
        )?;
    }

    if let Some(path) = cli.out_tu_group_count_matrix.as_ref() {
        if counts.has_groups() {
            write_entity_metric_matrix(
                path,
                None,
                "tu_id",
                counts.tu_ids(),
                counts.groups(),
                counts.group_unique_counts(),
                counts.group_full_length_evidence_counts(),
                counts.group_counts(),
                counts.group_fractional_counts(),
            )?;
        }
    }

    Ok(())
}

fn write_gene_outputs(
    cli: &ExecutionConfig,
    genes: &[GeneRecord],
    overlaps: &[Vec<(usize, u32)>],
    tu_counts: &CountMetrics,
    multi_counts: Option<&MultiSampleCounts>,
) -> anyhow::Result<()> {
    let gene_ids: Vec<String> = genes.iter().map(|gene| gene.id.clone()).collect();

    if let Some(path) = cli.out_gene_count.as_ref() {
        let gene_total_counts = build_gene_total_counts(genes, overlaps, tu_counts);
        write_gene_count_metrics_csv(path, &gene_ids, &gene_total_counts)?;
    }

    if let Some(multi_counts) = multi_counts {
        let (gene_sample_counts, gene_group_counts) =
            build_gene_multi_counts(genes, overlaps, multi_counts);

        if let Some(path) = cli.out_gene_sample_count_matrix.as_ref() {
            write_entity_metric_matrix(
                path,
                Some(GENE_COUNT_SEMANTICS),
                "gene_id",
                &gene_ids,
                multi_counts.samples(),
                &gene_sample_counts.unique,
                &gene_sample_counts.full_length_evidence,
                &gene_sample_counts.total,
                &gene_sample_counts.fractional,
            )?;
        }

        if let Some(path) = cli.out_gene_group_count_matrix.as_ref() {
            if multi_counts.has_groups() {
                write_entity_metric_matrix(
                    path,
                    Some(GENE_COUNT_SEMANTICS),
                    "gene_id",
                    &gene_ids,
                    multi_counts.groups(),
                    &gene_group_counts.unique,
                    &gene_group_counts.full_length_evidence,
                    &gene_group_counts.total,
                    &gene_group_counts.fractional,
                )?;
            }
        }
    }

    Ok(())
}

fn write_empty_cluster_outputs(
    cli: &ExecutionConfig,
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
    let mut membership_writer = create_membership_writer(
        out_membership,
        cli.ambiguity_margin,
        cli.fractional_assignment,
    )?;
    membership_writer.flush()?;

    if let Some(out_endpoint_stats) = cli.out_tu_endpoint_stats.as_ref() {
        write_tu_endpoint_stats(
            out_endpoint_stats,
            &PreparedOutputs {
                tus: Vec::new(),
                endpoint_stats: Vec::new(),
                id_mappings: Vec::new(),
            },
            &[],
        )?;
    }
    if let Some(out_tu_id_map) = cli.out_tu_id_map.as_ref() {
        write_tu_id_map(
            out_tu_id_map,
            &PreparedOutputs {
                tus: Vec::new(),
                endpoint_stats: Vec::new(),
                id_mappings: Vec::new(),
            },
        )?;
    }

    if let Some(out_tu_count) = cli.out_tu_count.as_ref() {
        write_count_metrics_csv(out_tu_count, "tu_id", &[], &CountMetrics::zeros(0))?;
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

    if cli.annotation_bed.is_some() && want_annotation_outputs(cli) {
        let genes = if let Some(annotation_path) = cli.annotation_bed.as_ref() {
            read_annotation_bed6(annotation_path)?
        } else {
            Vec::new()
        };
        let policy = gene_overlap_policy(cli);
        let contexts = tu_gene_contexts(&[], &genes, policy);
        let classifications = classify_tus(&[], &contexts);
        if let Some(path) = cli.out_tu_gene.as_ref() {
            write_tu_gene_relationships(path, &[], &genes, &contexts, policy)?;
        }
        if let Some(path) = cli.out_tu_semantics.as_ref() {
            write_tu_semantics(path, &[], &genes, &contexts, &classifications, policy)?;
        }
        if let Some(path) = cli.out_tu_gff3.as_ref() {
            write_tu_gff3(path, &[], &genes, &contexts, &classifications, &[])?;
        }
        let empty_multi_counts =
            samples.map(|samples| MultiSampleCounter::new(samples, Vec::new()).finish());
        if want_gene_outputs(cli) {
            write_gene_outputs(
                cli,
                &genes,
                &contexts.same_strand,
                &CountMetrics::zeros(0),
                empty_multi_counts.as_ref(),
            )?;
        }
    }

    Ok(())
}

fn run_cluster_mode(
    cli: &ExecutionConfig,
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

    if want_annotation_outputs(cli) && cli.annotation_bed.is_none() {
        anyhow::bail!("--annotation-bed is required for TU semantics, GFF3, and gene outputs");
    }

    let t_read_start = Instant::now();
    let mut input_outcome: ReadLoadOutcome = match samples {
        Some(samples) => {
            let manifest_format = resolve_manifest_input_format(cli.format, samples)?;
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
            let input = cli
                .input
                .as_ref()
                .expect("single-input mode must have --in");
            let input_format = resolve_single_input_format(cli.format, input);
            read_reads(input, input_format)?
        }
    };
    filter_reads(&mut input_outcome, cli.min_read_len);
    if cli.strict_read_errors && !input_outcome.rejections.is_empty() {
        let first = input_outcome
            .rejections
            .first()
            .expect("non-empty rejection list checked above");
        anyhow::bail!(
            "strict read-error policy rejected {} read record(s); first rejection: {}",
            input_outcome.rejections.len(),
            first.context()
        );
    }
    let total_input_records = input_outcome.total_records;
    let rejections = std::mem::take(&mut input_outcome.rejections);
    let (reads, full_length_evidence_by_read): (Vec<crate::tu::ReadRecord>, Vec<bool>) =
        std::mem::take(&mut input_outcome.reads)
            .into_iter()
            .map(|located| (located.read, located.full_length_evidence))
            .unzip();
    let rejection_path = cli
        .out_read_rejections
        .as_deref()
        .expect("cluster output defaults always include a read rejection report");
    write_read_rejections(rejection_path, &rejections)?;
    emit_read_input_summary(total_input_records, reads.len(), &rejections);
    if reads.is_empty() && !rejections.is_empty() {
        eprintln!(
            "warning: every input read was rejected; publishing an empty clustering result and the complete rejection report"
        );
    }
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

    eprintln!(
        "pipeline_stage\tstage=clustering\tstatus=started\treads={}\tthreads={threads}",
        reads.len()
    );
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
    eprintln!(
        "pipeline_stage\tstage=clustering\tstatus=completed\treads={}\ttus={}\telapsed_seconds={:.3}",
        reads.len(),
        result.tus().len(),
        t_cluster.as_secs_f64()
    );

    let prepared = prepare_outputs(&result, cli.min_tu_count, cli.tu_id_style)?;
    eprintln!(
        "pipeline_stage\tstage=assignment\tstatus=started\treads={}\ttus={}",
        reads.len(),
        prepared.tus.len()
    );
    let t_assignment_start = Instant::now();
    let assignments = crate::tu::assign_reads_to_tus(
        &reads,
        &prepared.tus,
        cli.score1_threshold,
        cli.score2_threshold,
        clustering_options,
        cli.ambiguity_margin,
        cli.fractional_assignment,
    )?;
    let t_assignment = t_assignment_start.elapsed();
    eprintln!(
        "pipeline_stage\tstage=assignment\tstatus=completed\treads={}\ttus={}\telapsed_seconds={:.3}",
        reads.len(),
        prepared.tus.len(),
        t_assignment.as_secs_f64()
    );
    let mut count_metrics = CountMetrics::zeros(prepared.tus.len());

    let t_write_tu_start = Instant::now();
    let tu_records: Vec<crate::io::bed::Bed6Record> = prepared
        .tus
        .iter()
        .map(|tu| crate::io::bed::Bed6Record {
            chrom: tu.contig.clone(),
            start: tu.interval.start(),
            end: tu.interval.end(),
            name: tu.id.clone(),
            score: 0,
            strand: tu.strand,
            extra_fields: Vec::new(),
        })
        .collect();
    crate::io::bed::write_bed6(out_tu, tu_records.iter())?;
    if let Some(out_endpoint_stats) = cli.out_tu_endpoint_stats.as_ref() {
        write_tu_endpoint_stats(out_endpoint_stats, &prepared, &reads)?;
    }
    if let Some(out_tu_id_map) = cli.out_tu_id_map.as_ref() {
        write_tu_id_map(out_tu_id_map, &prepared)?;
    }
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
    let mut writer = create_membership_writer(
        out_membership,
        cli.ambiguity_margin,
        cli.fractional_assignment,
    )?;
    for read_idx in read_indices {
        let read = &reads[read_idx];
        let assignment = &assignments[read_idx];
        let full_length_evidence = full_length_evidence_by_read[read_idx];

        if let Some(tu_idx) = hard_assignment_tu_index(assignment) {
            count_metrics.total[tu_idx] += 1;
            if assignment.status == AssignmentStatus::Unique {
                count_metrics.unique[tu_idx] += 1;
            }
            if full_length_evidence {
                count_metrics.full_length_evidence[tu_idx] += 1;
            }
            if let Some(members) = tu_members.as_mut() {
                members[tu_idx].push(read_idx);
            }
        }
        for fraction in &assignment.fractional_assignments {
            count_metrics.fractional[fraction.tu_index] += fraction.weight;
        }

        if let Some(counter) = multi_counter.as_mut() {
            if let Some(tu_idx) = hard_assignment_tu_index(assignment) {
                let weight = assignment
                    .fractional_assignments
                    .iter()
                    .find(|fraction| fraction.tu_index == tu_idx)
                    .map(|fraction| fraction.weight)
                    .unwrap_or(1.0);
                counter.add_assignment_metrics(
                    &read.id,
                    &prepared.tus[tu_idx].id,
                    assignment.status == AssignmentStatus::Unique,
                    full_length_evidence,
                    true,
                    weight,
                )?;
            } else {
                for fraction in &assignment.fractional_assignments {
                    counter.add_assignment_metrics(
                        &read.id,
                        &prepared.tus[fraction.tu_index].id,
                        false,
                        false,
                        false,
                        fraction.weight,
                    )?;
                }
            }
        }

        write_membership_row(
            &mut writer,
            read,
            assignment,
            &prepared.tus,
            full_length_evidence,
        )?;
    }
    writer.flush()?;
    let t_membership = t_membership_start.elapsed();

    if let Some(out_tu_count) = cli.out_tu_count.as_ref() {
        let t_tu_count_start = Instant::now();
        let tu_ids: Vec<String> = prepared.tus.iter().map(|tu| tu.id.clone()).collect();
        write_count_metrics_csv(out_tu_count, "tu_id", &tu_ids, &count_metrics)?;
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

    let needs_annotation = cli.annotation_bed.is_some() && want_annotation_outputs(cli);
    if needs_annotation {
        let t_annotation_start = Instant::now();
        let annotation_path = cli.annotation_bed.as_ref().expect("checked");
        let genes = read_annotation_bed6(annotation_path)?;
        let policy = gene_overlap_policy(cli);
        let contexts = tu_gene_contexts(&prepared.tus, &genes, policy);
        let classifications = classify_tus(&prepared.tus, &contexts);

        if let Some(out_tu_gene) = cli.out_tu_gene.as_ref() {
            write_tu_gene_relationships(out_tu_gene, &prepared.tus, &genes, &contexts, policy)?;
        }
        if let Some(out_tu_semantics) = cli.out_tu_semantics.as_ref() {
            write_tu_semantics(
                out_tu_semantics,
                &prepared.tus,
                &genes,
                &contexts,
                &classifications,
                policy,
            )?;
        }
        if let Some(out_tu_gff3) = cli.out_tu_gff3.as_ref() {
            write_tu_gff3(
                out_tu_gff3,
                &prepared.tus,
                &genes,
                &contexts,
                &classifications,
                &prepared.endpoint_stats,
            )?;
        }

        if let Some(out_tu_bed12) = cli.out_tu_bed12.as_ref() {
            let tu_members = tu_members
                .as_ref()
                .expect("tu_members must be built when --out-tu-bed12 is set");
            let mut transcripts: Vec<Transcript> = Vec::with_capacity(prepared.tus.len());
            for (tu_idx, tu) in prepared.tus.iter().enumerate() {
                let read_count = tu_members[tu_idx].len();
                let name2 = format_name2(&tu_members[tu_idx], &reads, read_count);

                let gene_names: Vec<String> = contexts.same_strand[tu_idx]
                    .iter()
                    .map(|(gene_idx, _)| genes[*gene_idx].id.clone())
                    .collect();

                let exons: Vec<Interval> = if gene_names.is_empty() {
                    vec![tu.interval]
                } else {
                    merge_overlapping_intervals(
                        contexts.same_strand[tu_idx]
                            .iter()
                            .filter_map(|(gene_idx, _)| {
                                tu.interval.intersection(genes[*gene_idx].interval)
                            })
                            .collect(),
                    )
                };

                let gene_field = if gene_names.is_empty() {
                    ".".to_owned()
                } else {
                    gene_names.join(",")
                };

                let transcript = Transcript::new(
                    tu.contig.clone(),
                    tu.strand,
                    tu.interval.start(),
                    tu.interval.end(),
                    tu.id.clone(),
                    exons,
                    Bed12Attrs {
                        score: 0,
                        thick_start: tu.interval.start(),
                        thick_end: tu.interval.end(),
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
                &contexts.same_strand,
                &count_metrics,
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
        eprintln!("[timings] assign_reads_to_tus={t_assignment:?}");
        eprintln!("[timings] write_tu_bed={t_write_tu:?}");
        eprintln!("[timings] write_membership={t_membership:?}");
        eprintln!("[timings] total={:?}", total_start.elapsed());
    }

    Ok(())
}

fn run_counts_only_mode(
    cli: &ExecutionConfig,
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
        write_count_metrics_csv(out_tu_count, "tu_id", counts.tu_ids(), &totals)?;
    }
    write_multi_sample_outputs(cli, &counts)?;

    if cli.timings {
        eprintln!("[timings] threads={threads}");
        eprintln!(
            "[timings] count_from_membership={t_count:?} (tus={})",
            counts.tu_ids().len()
        );
        eprintln!("[timings] total={:?}", total_start.elapsed());
    }

    Ok(())
}

fn execute_config(mut cli: ExecutionConfig) -> anyhow::Result<()> {
    let mode = resolve_mode(&cli)?;
    apply_default_outputs(&mut cli, &mode);
    let transaction =
        OutputTransaction::new(command_input_paths(&cli, &mode), command_output_paths(&cli))?;
    let staged_config = staged_execution_config(&cli, &transaction)?;
    prepare_output_dirs(&staged_config)?;

    let threads = cli.threads.unwrap_or_else(default_threads);
    if threads == 0 {
        anyhow::bail!("--threads must be >= 1");
    }

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| anyhow::anyhow!("failed to build thread pool (threads={threads}): {e}"))?;
    pool.install(|| match &mode {
        InputMode::Single => run_cluster_mode(&staged_config, threads, None),
        InputMode::ManifestCluster { samples } => {
            run_cluster_mode(&staged_config, threads, Some(samples.as_slice()))
        }
        InputMode::ManifestCountOnly {
            samples,
            membership,
        } => run_counts_only_mode(
            &staged_config,
            threads,
            samples.as_slice(),
            membership.as_path(),
        ),
    })?;
    transaction.commit()
}

#[derive(Error, Debug)]
pub(crate) enum PipelineError {
    #[error(transparent)]
    Config(#[from] ConfigError),

    #[error("cluster/recount execution failed: {0}")]
    Execution(#[source] Box<dyn std::error::Error + Send + Sync>),
}

pub(crate) fn run_cluster(config: ClusterConfig) -> Result<(), PipelineError> {
    execute_config(config.into())
        .map_err(|error| PipelineError::Execution(error.into_boxed_dyn_error()))
}

pub(crate) fn run_recount(config: RecountConfig) -> Result<(), PipelineError> {
    execute_config(config.into())
        .map_err(|error| PipelineError::Execution(error.into_boxed_dyn_error()))
}

pub(crate) fn run_cluster_from_args<I, T>(args: I) -> anyhow::Result<()>
where
    I: IntoIterator<Item = T>,
    T: Into<std::ffi::OsString> + Clone,
{
    let cli = ClusterCli::parse_from(args);
    let config = ClusterConfig::try_from(cli)?;
    run_cluster(config).map_err(anyhow::Error::new)
}

pub(crate) fn run_recount_from_args<I, T>(args: I) -> anyhow::Result<()>
where
    I: IntoIterator<Item = T>,
    T: Into<std::ffi::OsString> + Clone,
{
    let cli = RecountCli::parse_from(args);
    let config = RecountConfig::try_from(cli)?;
    run_recount(config).map_err(anyhow::Error::new)
}

#[cfg(test)]
mod tests {
    use std::sync::atomic::{AtomicU64, Ordering};

    use super::*;

    static NEXT_TMP_ID: AtomicU64 = AtomicU64::new(0);

    fn unique_tmp_dir(label: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "trackclustertu_{label}_{}_{}",
            std::process::id(),
            NEXT_TMP_ID.fetch_add(1, Ordering::Relaxed)
        ))
    }

    #[test]
    fn typed_cluster_api_matches_cli_boundary_byte_for_byte() {
        let root = unique_tmp_dir("typed_cli_parity");
        let typed_out = root.join("typed");
        let cli_out = root.join("cli");
        let input = root.join("reads.bed");
        fs::create_dir_all(&root).unwrap();
        fs::write(
            &input,
            concat!(
                "chr1\t100\t200\tr1\t0\t+\n",
                "chr1\t101\t201\tr2\t0\t+\n",
                "chr1\t400\t500\tr3\t0\t-\n"
            ),
        )
        .unwrap();

        let typed = ClusterConfig::validate(ClusterRequest {
            input: ClusterInput::Single(input.clone()),
            format: InputFormat::Bed6,
            outputs: ClusterOutputPaths {
                out_dir: Some(typed_out.clone()),
                ..ClusterOutputPaths::default()
            },
            options: ClusterOptions {
                tu_id_style: TuIdStyle::Sequential,
                ..ClusterOptions::default()
            },
            min_read_len: None,
            min_tu_count: None,
            annotation: None,
            threads: Some(1),
            timings: false,
        })
        .unwrap();
        run_cluster(typed).unwrap();

        run_cluster_from_args([
            std::ffi::OsString::from("trackclustertu cluster"),
            std::ffi::OsString::from("--in"),
            input.as_os_str().to_owned(),
            std::ffi::OsString::from("--format"),
            std::ffi::OsString::from("bed6"),
            std::ffi::OsString::from("--out-dir"),
            cli_out.as_os_str().to_owned(),
            std::ffi::OsString::from("--tu-id-style"),
            std::ffi::OsString::from("sequential"),
            std::ffi::OsString::from("--threads"),
            std::ffi::OsString::from("1"),
        ])
        .unwrap();

        for filename in [
            "tus.bed",
            "membership.tsv",
            "tu_endpoint_stats.tsv",
            "tu_id_map.tsv",
            "tu_count.csv",
            "read_rejections.tsv",
        ] {
            assert_eq!(
                fs::read(typed_out.join(filename)).unwrap(),
                fs::read(cli_out.join(filename)).unwrap(),
                "typed and CLI output differ for {filename}"
            );
        }

        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn typed_recount_api_matches_cli_boundary_byte_for_byte() {
        let root = unique_tmp_dir("typed_recount_cli_parity");
        let typed_out = root.join("typed");
        let cli_out = root.join("cli");
        let manifest = root.join("samples.tsv");
        let membership = root.join("membership.tsv");
        fs::create_dir_all(&root).unwrap();
        fs::write(
            &manifest,
            concat!(
                "sample\treads\tgroup\n",
                "sampleA\ta.bed\tcontrol\n",
                "sampleB\tb.bed\tcontrol\n"
            ),
        )
        .unwrap();
        fs::write(
            &membership,
            concat!(
                "sampleA::r1\tTU1\t1.0\t1.0\n",
                "sampleB::r2\tTU1\t1.0\t1.0\n",
                "sampleB::r3\tTU2\t1.0\t1.0\n"
            ),
        )
        .unwrap();

        let typed = RecountConfig::validate(RecountRequest {
            manifest: manifest.clone(),
            membership: membership.clone(),
            outputs: RecountOutputPaths {
                out_dir: Some(typed_out.clone()),
                ..RecountOutputPaths::default()
            },
            min_tu_count: None,
            threads: Some(1),
            timings: false,
        })
        .unwrap();
        run_recount(typed).unwrap();

        run_recount_from_args([
            std::ffi::OsString::from("trackclustertu recount"),
            std::ffi::OsString::from("--manifest"),
            manifest.as_os_str().to_owned(),
            std::ffi::OsString::from("--pooled-membership"),
            membership.as_os_str().to_owned(),
            std::ffi::OsString::from("--out-dir"),
            cli_out.as_os_str().to_owned(),
            std::ffi::OsString::from("--threads"),
            std::ffi::OsString::from("1"),
        ])
        .unwrap();

        for filename in [
            "tu_count.csv",
            "tu_sample_long.tsv",
            "tu_sample_matrix.tsv",
            "tu_group_matrix.tsv",
        ] {
            assert_eq!(
                fs::read(typed_out.join(filename)).unwrap(),
                fs::read(cli_out.join(filename)).unwrap(),
                "typed and CLI recount output differ for {filename}"
            );
        }

        fs::remove_dir_all(root).unwrap();
    }
}
