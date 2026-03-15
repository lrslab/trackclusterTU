use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use anyhow::{Context, Result};
use clap::Parser;

use crate::tu::multi::{parse_manifest, SampleManifestRecord};

const TOP_LEVEL_HELP: &str = "\
trackclustertu

Usage:
  trackclustertu <command> [options]

Commands:
  run         Full pipeline from FASTQ manifest to TU/gene counts
  map         FASTQ manifest to sorted BAMs plus BED manifest
  bam-to-bed  BAM input(s) to BED6
  cluster     BED input(s) to TU/gene outputs
  recount     Membership TSV to count tables
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
    #[arg(long)]
    threads: Option<usize>,

    /// Minimum MAPQ kept during BAM generation.
    #[arg(long, default_value_t = 0)]
    min_mapq: u8,

    /// Additional minimap2 arguments provided as one whitespace-separated string.
    #[arg(long, default_value = "-ax map-ont")]
    minimap2_args: String,

    /// score1 threshold (overlap / union).
    #[arg(long, default_value_t = 0.95)]
    score1_threshold: f64,

    /// score2 threshold (overlap / min_len).
    #[arg(long, default_value_t = 0.99)]
    score2_threshold: f64,

    /// Skip the second-pass containment attachment and keep score1 seed clusters as final TUs.
    #[arg(long)]
    skip_score2_attachment: bool,

    /// Optional minimum read length filter (bp).
    #[arg(long)]
    min_read_len: Option<u32>,

    /// Optional minimum reads per TU (filters outputs).
    #[arg(long)]
    min_tu_count: Option<u64>,

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
    #[arg(long)]
    threads: Option<usize>,

    /// Minimum MAPQ kept during BAM generation.
    #[arg(long, default_value_t = 0)]
    min_mapq: u8,

    /// Additional minimap2 arguments provided as one whitespace-separated string.
    #[arg(long, default_value = "-ax map-ont")]
    minimap2_args: String,
}

#[derive(Parser, Debug, Clone)]
#[command(
    name = "trackclustertu bam-to-bed",
    version,
    about = "Convert BAM input(s) to BED6"
)]
struct BamToBedCli {
    /// Single input BAM file.
    #[arg(long = "in-bam", conflicts_with = "manifest", requires = "out_bed")]
    input_bam: Option<PathBuf>,

    /// BAM manifest TSV with columns: sample, reads, [group].
    #[arg(long, conflicts_with = "input_bam", requires = "out_dir")]
    manifest: Option<PathBuf>,

    /// Output BED6 path for single-BAM mode.
    #[arg(long = "out-bed", requires = "input_bam")]
    out_bed: Option<PathBuf>,

    /// Output directory for manifest mode.
    #[arg(long, requires = "manifest")]
    out_dir: Option<PathBuf>,
}

#[derive(Parser, Debug, Clone)]
#[command(
    name = "trackclustertu gff-to-bed",
    version,
    about = "Convert GFF3 gene annotations to BED6"
)]
struct GffToBedCli {
    #[arg(long = "annotation-gff", alias = "gff")]
    annotation_gff: PathBuf,

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
    let mut writer = BufWriter::new(
        File::create(manifest_path)
            .with_context(|| format!("failed to create manifest {:?}", manifest_path))?,
    );
    writeln!(writer, "sample\tgroup\treads")?;
    for record in records {
        writeln!(
            writer,
            "{}\t{}\t{}",
            record.sample,
            record.group.as_deref().unwrap_or(""),
            record.reads.display()
        )?;
    }
    Ok(())
}

fn split_whitespace_args(raw: &str) -> Vec<String> {
    raw.split_whitespace().map(str::to_owned).collect()
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

struct MapConfig<'a> {
    reference_fasta: &'a Path,
    logs_dir: &'a Path,
    sample_stem: &'a str,
    threads: usize,
    min_mapq: u8,
    minimap2_args: &'a str,
}

fn map_fastq_to_bam(config: &MapConfig<'_>, fastq_path: &Path, bam_path: &Path) -> Result<()> {
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

    let mut minimap2 = Command::new("minimap2");
    for arg in split_whitespace_args(config.minimap2_args) {
        minimap2.arg(arg);
    }
    minimap2
        .arg("-t")
        .arg(config.threads.to_string())
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

    let mut samtools_view = Command::new("samtools");
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

    let mut samtools_sort = Command::new("samtools");
    samtools_sort
        .arg("sort")
        .arg("-@")
        .arg(config.threads.to_string())
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

    let index_status = Command::new("samtools")
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

fn run_map(cli: &MapCli) -> Result<PathBuf> {
    ensure_path_exists(&cli.manifest, "manifest")?;
    ensure_path_exists(&cli.reference_fasta, "reference FASTA")?;

    let threads = cli.threads.unwrap_or_else(default_threads);
    if threads == 0 {
        anyhow::bail!("--threads must be >= 1");
    }

    let records = parse_manifest(&cli.manifest)
        .with_context(|| format!("failed to parse FASTQ manifest {:?}", cli.manifest))?;

    fs::create_dir_all(&cli.out_dir)
        .with_context(|| format!("failed to create output directory {:?}", cli.out_dir))?;
    let bam_dir = cli.out_dir.join("bam");
    let logs_dir = cli.out_dir.join("logs");
    fs::create_dir_all(&bam_dir)
        .with_context(|| format!("failed to create BAM output directory {:?}", bam_dir))?;
    fs::create_dir_all(&logs_dir)
        .with_context(|| format!("failed to create log output directory {:?}", logs_dir))?;

    let mut bam_records: Vec<SampleManifestRecord> = Vec::with_capacity(records.len());
    for record in &records {
        ensure_path_exists(&record.reads, "FASTQ")?;
        let sample_stem = crate::tools::bam_to_bed::sanitize_sample_name(&record.sample);
        let bam_path = bam_dir.join(format!("{sample_stem}.sorted.bam"));
        let config = MapConfig {
            reference_fasta: &cli.reference_fasta,
            logs_dir: &logs_dir,
            sample_stem: &sample_stem,
            threads,
            min_mapq: cli.min_mapq,
            minimap2_args: &cli.minimap2_args,
        };
        map_fastq_to_bam(&config, &record.reads, &bam_path)?;
        let bam_path = bam_path
            .canonicalize()
            .with_context(|| format!("failed to canonicalize BAM path {:?}", bam_path))?;

        bam_records.push(SampleManifestRecord {
            sample: record.sample.clone(),
            reads: bam_path,
            group: record.group.clone(),
        });
    }

    let bam_manifest = cli.out_dir.join("samples.bam.tsv");
    write_manifest(&bam_manifest, &bam_records)?;
    let bam_manifest = bam_manifest
        .canonicalize()
        .with_context(|| format!("failed to canonicalize BAM manifest {:?}", bam_manifest))?;

    let bed_manifest =
        crate::tools::bam_to_bed::convert_records_to_bed_manifest(&bam_records, &cli.out_dir)?;
    let bed_manifest = bed_manifest
        .canonicalize()
        .with_context(|| format!("failed to canonicalize BED manifest {:?}", bed_manifest))?;

    println!("bam_manifest={}", bam_manifest.display());
    println!("bed_manifest={}", bed_manifest.display());

    Ok(bed_manifest)
}

fn run_map_from_args<I>(args: I) -> Result<()>
where
    I: IntoIterator<Item = OsString>,
{
    let cli = MapCli::parse_from(args);
    run_map(&cli).map(|_| ())
}

fn run_bam_to_bed_from_args<I>(args: I) -> Result<()>
where
    I: IntoIterator<Item = OsString>,
{
    let cli = BamToBedCli::parse_from(args);
    match (cli.input_bam.as_deref(), cli.manifest.as_deref()) {
        (Some(input_bam), None) => {
            let out_bed = cli.out_bed.as_deref().expect("required by clap");
            crate::tools::bam_to_bed::convert_single_bam_to_bed(input_bam, out_bed)?;
            println!("out_bed={}", out_bed.display());
            Ok(())
        }
        (None, Some(manifest)) => {
            let out_dir = cli.out_dir.as_deref().expect("required by clap");
            let bed_manifest =
                crate::tools::bam_to_bed::convert_manifest_to_bed_manifest(manifest, out_dir)?;
            println!("bed_manifest={}", bed_manifest.display());
            Ok(())
        }
        _ => anyhow::bail!("pass either --in-bam/--out-bed or --manifest/--out-dir"),
    }
}

fn run_gff_to_bed_from_args<I>(args: I) -> Result<()>
where
    I: IntoIterator<Item = OsString>,
{
    let cli = GffToBedCli::parse_from(args);
    crate::tools::gff_to_bed::convert_gff_to_bed(&cli.annotation_gff, &cli.out_bed)
}

fn push_optional_arg(args: &mut Vec<OsString>, flag: &str, value: Option<impl ToString>) {
    if let Some(value) = value {
        args.push(OsString::from(flag));
        args.push(OsString::from(value.to_string()));
    }
}

fn push_optional_path(args: &mut Vec<OsString>, flag: &str, value: Option<&Path>) {
    if let Some(value) = value {
        args.push(OsString::from(flag));
        args.push(value.as_os_str().to_owned());
    }
}

fn run_full_pipeline(cli: &RunCli) -> Result<()> {
    let map_cli = MapCli {
        manifest: cli.manifest.clone(),
        reference_fasta: cli.reference_fasta.clone(),
        out_dir: cli.out_dir.clone(),
        threads: cli.threads,
        min_mapq: cli.min_mapq,
        minimap2_args: cli.minimap2_args.clone(),
    };
    let bed_manifest = run_map(&map_cli)?;

    let annotation_bed = match (&cli.annotation_bed, &cli.annotation_gff) {
        (Some(path), None) => Some(path.clone()),
        (None, Some(gff)) => {
            let out_bed = cli.out_dir.join("annotation.bed");
            crate::tools::gff_to_bed::convert_gff_to_bed(gff, &out_bed)?;
            Some(
                out_bed
                    .canonicalize()
                    .with_context(|| format!("failed to canonicalize BED path {:?}", out_bed))?,
            )
        }
        (None, None) => None,
        (Some(_), Some(_)) => unreachable!("clap enforces conflict"),
    };

    let mut cluster_args = vec![
        OsString::from("trackclustertu cluster"),
        OsString::from("--manifest"),
        bed_manifest.as_os_str().to_owned(),
        OsString::from("--format"),
        OsString::from("bed6"),
        OsString::from("--out-dir"),
        cli.out_dir.as_os_str().to_owned(),
    ];
    push_optional_path(
        &mut cluster_args,
        "--annotation-bed",
        annotation_bed.as_deref(),
    );
    push_optional_arg(&mut cluster_args, "--threads", cli.threads);
    push_optional_arg(
        &mut cluster_args,
        "--score1-threshold",
        Some(cli.score1_threshold),
    );
    push_optional_arg(
        &mut cluster_args,
        "--score2-threshold",
        Some(cli.score2_threshold),
    );
    if cli.skip_score2_attachment {
        cluster_args.push(OsString::from("--skip-score2-attachment"));
    }
    push_optional_arg(&mut cluster_args, "--min-read-len", cli.min_read_len);
    push_optional_arg(&mut cluster_args, "--min-tu-count", cli.min_tu_count);
    if cli.timings {
        cluster_args.push(OsString::from("--timings"));
    }

    crate::tools::cluster_pipeline::run_cluster_from_args(cluster_args)
}

fn run_full_from_args<I>(args: I) -> Result<()>
where
    I: IntoIterator<Item = OsString>,
{
    let cli = RunCli::parse_from(args);
    run_full_pipeline(&cli)
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
            println!("{}", env!("CARGO_PKG_VERSION"));
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
