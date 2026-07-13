use std::path::PathBuf;

use clap::ValueEnum;
use thiserror::Error;

#[derive(ValueEnum, Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum InputFormat {
    #[default]
    Auto,
    Bed6,
    Bed12,
    Tsv,
}

#[derive(ValueEnum, Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum TuIdStyle {
    /// Coordinate-derived IDs that remain unchanged when unrelated TUs are filtered.
    #[default]
    Stable,
    /// Historical TU000001-style IDs assigned after filtering.
    Sequential,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum ClusterInput {
    Single(PathBuf),
    Manifest(PathBuf),
}

impl ClusterInput {
    pub(crate) fn is_manifest(&self) -> bool {
        matches!(self, Self::Manifest(_))
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ClusterOptions {
    pub(crate) span_jaccard_threshold: f64,
    pub(crate) overlap_over_longer_threshold: f64,
    pub(crate) three_prime_tolerance_bp: u32,
    pub(crate) max_five_prime_delta_bp: Option<u32>,
    pub(crate) attach_contained_reads: bool,
    pub(crate) ambiguity_margin: f64,
    pub(crate) fractional_assignment: bool,
    pub(crate) tu_id_style: TuIdStyle,
    pub(crate) strict_read_errors: bool,
}

impl Default for ClusterOptions {
    fn default() -> Self {
        Self {
            span_jaccard_threshold: 0.95,
            overlap_over_longer_threshold: 0.80,
            three_prime_tolerance_bp: 12,
            max_five_prime_delta_bp: None,
            attach_contained_reads: true,
            ambiguity_margin: 0.02,
            fractional_assignment: false,
            tu_id_style: TuIdStyle::Stable,
            strict_read_errors: false,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct AnnotationRequest {
    pub(crate) bed: PathBuf,
    pub(crate) min_overlap_bp: u32,
    pub(crate) min_tu_fraction: f64,
    pub(crate) min_gene_fraction: f64,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct ClusterOutputPaths {
    pub(crate) out_dir: Option<PathBuf>,
    pub(crate) tu: Option<PathBuf>,
    pub(crate) membership: Option<PathBuf>,
    pub(crate) pooled_reads: Option<PathBuf>,
    pub(crate) read_rejections: Option<PathBuf>,
    pub(crate) endpoint_stats: Option<PathBuf>,
    pub(crate) id_map: Option<PathBuf>,
    pub(crate) tu_count: Option<PathBuf>,
    pub(crate) tu_sample_long: Option<PathBuf>,
    pub(crate) tu_sample_matrix: Option<PathBuf>,
    pub(crate) tu_group_matrix: Option<PathBuf>,
    pub(crate) tu_gene: Option<PathBuf>,
    pub(crate) tu_semantics: Option<PathBuf>,
    pub(crate) tu_gff3: Option<PathBuf>,
    pub(crate) tu_bed12: Option<PathBuf>,
    pub(crate) gene_count: Option<PathBuf>,
    pub(crate) gene_sample_matrix: Option<PathBuf>,
    pub(crate) gene_group_matrix: Option<PathBuf>,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ClusterRequest {
    pub(crate) input: ClusterInput,
    pub(crate) format: InputFormat,
    pub(crate) outputs: ClusterOutputPaths,
    pub(crate) options: ClusterOptions,
    pub(crate) min_read_len: Option<u32>,
    pub(crate) min_tu_count: Option<u64>,
    pub(crate) annotation: Option<AnnotationRequest>,
    pub(crate) threads: Option<usize>,
    pub(crate) timings: bool,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ClusterConfig {
    pub(super) input: ClusterInput,
    pub(super) format: InputFormat,
    pub(super) outputs: ClusterOutputPaths,
    pub(super) options: ClusterOptions,
    pub(super) min_read_len: Option<u32>,
    pub(super) min_tu_count: Option<u64>,
    pub(super) annotation: Option<AnnotationRequest>,
    pub(super) threads: usize,
    pub(super) timings: bool,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct RecountOutputPaths {
    pub(crate) out_dir: Option<PathBuf>,
    pub(crate) tu_count: Option<PathBuf>,
    pub(crate) tu_sample_long: Option<PathBuf>,
    pub(crate) tu_sample_matrix: Option<PathBuf>,
    pub(crate) tu_group_matrix: Option<PathBuf>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RecountRequest {
    pub(crate) manifest: PathBuf,
    pub(crate) membership: PathBuf,
    pub(crate) outputs: RecountOutputPaths,
    pub(crate) min_tu_count: Option<u64>,
    pub(crate) threads: Option<usize>,
    pub(crate) timings: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RecountConfig {
    pub(super) manifest: PathBuf,
    pub(super) membership: PathBuf,
    pub(super) outputs: RecountOutputPaths,
    pub(super) min_tu_count: Option<u64>,
    pub(super) threads: usize,
    pub(super) timings: bool,
}

#[derive(Error, Debug, PartialEq)]
pub(crate) enum ConfigError {
    #[error("exactly one of --in or --manifest is required")]
    MissingInput,

    #[error("--in and --manifest are mutually exclusive")]
    ConflictingInput,

    #[error("{field} must be finite and between 0 and 1, got {value}")]
    InvalidUnitInterval { field: &'static str, value: f64 },

    #[error("--threads must be >= 1")]
    ZeroThreads,

    #[error("{output} requires manifest input")]
    ManifestOutputWithSingleInput { output: &'static str },

    #[error("{output} requires annotation input")]
    AnnotationOutputWithoutAnnotation { output: &'static str },
}

fn default_threads() -> usize {
    std::thread::available_parallelism()
        .map(|threads| threads.get())
        .unwrap_or(1)
}

fn validate_unit_interval(field: &'static str, value: f64) -> Result<(), ConfigError> {
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        return Err(ConfigError::InvalidUnitInterval { field, value });
    }
    Ok(())
}

impl ClusterConfig {
    pub(crate) fn validate(request: ClusterRequest) -> Result<Self, ConfigError> {
        validate_unit_interval(
            "--span-jaccard-threshold",
            request.options.span_jaccard_threshold,
        )?;
        validate_unit_interval(
            "--overlap-over-longer-threshold",
            request.options.overlap_over_longer_threshold,
        )?;
        validate_unit_interval("--ambiguity-margin", request.options.ambiguity_margin)?;
        if let Some(annotation) = request.annotation.as_ref() {
            validate_unit_interval("--gene-min-tu-fraction", annotation.min_tu_fraction)?;
            validate_unit_interval("--gene-min-gene-fraction", annotation.min_gene_fraction)?;
        }

        let threads = request.threads.unwrap_or_else(default_threads);
        if threads == 0 {
            return Err(ConfigError::ZeroThreads);
        }

        if !request.input.is_manifest() {
            for (path, output) in [
                (&request.outputs.pooled_reads, "--out-pooled-reads"),
                (
                    &request.outputs.tu_sample_long,
                    "--out-tu-sample-count-long",
                ),
                (
                    &request.outputs.tu_sample_matrix,
                    "--out-tu-sample-count-matrix",
                ),
                (
                    &request.outputs.tu_group_matrix,
                    "--out-tu-group-count-matrix",
                ),
                (
                    &request.outputs.gene_sample_matrix,
                    "--out-gene-sample-count-matrix",
                ),
                (
                    &request.outputs.gene_group_matrix,
                    "--out-gene-group-count-matrix",
                ),
            ] {
                if path.is_some() {
                    return Err(ConfigError::ManifestOutputWithSingleInput { output });
                }
            }
        }

        if request.annotation.is_none() {
            for (path, output) in [
                (&request.outputs.tu_gene, "--out-tu-gene"),
                (&request.outputs.tu_semantics, "--out-tu-semantics"),
                (&request.outputs.tu_gff3, "--out-tu-gff3"),
                (&request.outputs.tu_bed12, "--out-tu-bed12"),
                (&request.outputs.gene_count, "--out-gene-count"),
                (
                    &request.outputs.gene_sample_matrix,
                    "--out-gene-sample-count-matrix",
                ),
                (
                    &request.outputs.gene_group_matrix,
                    "--out-gene-group-count-matrix",
                ),
            ] {
                if path.is_some() {
                    return Err(ConfigError::AnnotationOutputWithoutAnnotation { output });
                }
            }
        }

        Ok(Self {
            input: request.input,
            format: request.format,
            outputs: request.outputs,
            options: request.options,
            min_read_len: request.min_read_len,
            min_tu_count: request.min_tu_count,
            annotation: request.annotation,
            threads,
            timings: request.timings,
        })
    }
}

impl RecountConfig {
    pub(crate) fn validate(request: RecountRequest) -> Result<Self, ConfigError> {
        let threads = request.threads.unwrap_or_else(default_threads);
        if threads == 0 {
            return Err(ConfigError::ZeroThreads);
        }
        Ok(Self {
            manifest: request.manifest,
            membership: request.membership,
            outputs: request.outputs,
            min_tu_count: request.min_tu_count,
            threads,
            timings: request.timings,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base_request() -> ClusterRequest {
        ClusterRequest {
            input: ClusterInput::Single(PathBuf::from("reads.bed")),
            format: InputFormat::Bed6,
            outputs: ClusterOutputPaths::default(),
            options: ClusterOptions::default(),
            min_read_len: None,
            min_tu_count: None,
            annotation: None,
            threads: Some(1),
            timings: false,
        }
    }

    #[test]
    fn cluster_config_rejects_nonfinite_and_out_of_range_values() {
        for value in [f64::NAN, -0.1, 1.1] {
            let mut request = base_request();
            request.options.ambiguity_margin = value;
            assert!(matches!(
                ClusterConfig::validate(request),
                Err(ConfigError::InvalidUnitInterval {
                    field: "--ambiguity-margin",
                    ..
                })
            ));
        }
    }

    #[test]
    fn cluster_config_reports_canonical_score_field_names() {
        let mut request = base_request();
        request.options.span_jaccard_threshold = 1.01;
        assert_eq!(
            ClusterConfig::validate(request).unwrap_err(),
            ConfigError::InvalidUnitInterval {
                field: "--span-jaccard-threshold",
                value: 1.01,
            }
        );

        let mut request = base_request();
        request.options.overlap_over_longer_threshold = -0.01;
        assert_eq!(
            ClusterConfig::validate(request).unwrap_err(),
            ConfigError::InvalidUnitInterval {
                field: "--overlap-over-longer-threshold",
                value: -0.01,
            }
        );
    }

    #[test]
    fn cluster_config_enforces_output_requirements_without_clap() {
        let mut request = base_request();
        request.outputs.tu_sample_matrix = Some(PathBuf::from("matrix.tsv"));
        assert_eq!(
            ClusterConfig::validate(request).unwrap_err(),
            ConfigError::ManifestOutputWithSingleInput {
                output: "--out-tu-sample-count-matrix"
            }
        );

        let mut request = base_request();
        request.outputs.tu_gff3 = Some(PathBuf::from("tus.gff3"));
        assert_eq!(
            ClusterConfig::validate(request).unwrap_err(),
            ConfigError::AnnotationOutputWithoutAnnotation {
                output: "--out-tu-gff3"
            }
        );
    }

    #[test]
    fn typed_configs_reject_zero_threads() {
        let mut cluster = base_request();
        cluster.threads = Some(0);
        assert_eq!(
            ClusterConfig::validate(cluster).unwrap_err(),
            ConfigError::ZeroThreads
        );
        let recount = RecountRequest {
            manifest: PathBuf::from("samples.tsv"),
            membership: PathBuf::from("membership.tsv"),
            outputs: RecountOutputPaths::default(),
            min_tu_count: None,
            threads: Some(0),
            timings: false,
        };
        assert_eq!(
            RecountConfig::validate(recount).unwrap_err(),
            ConfigError::ZeroThreads
        );
    }
}
