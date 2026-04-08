use std::fs;
use std::num::NonZeroUsize;
#[cfg(unix)]
use std::os::unix::fs::PermissionsExt;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use noodles::sam::alignment::io::Write as _;
use noodles::{
    bam,
    core::Position,
    sam::{
        self,
        alignment::{
            record::{
                cigar::{op::Kind as CigarKind, Op as CigarOp},
                Flags, MappingQuality,
            },
            record_buf::{Cigar, QualityScores, Sequence},
            RecordBuf,
        },
        header::record::value::{map::ReferenceSequence, Map},
    },
};

fn unique_tmp_dir(prefix: &str) -> std::path::PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock")
        .as_nanos();
    std::env::temp_dir().join(format!("{prefix}_{nanos}"))
}

fn build_header() -> sam::Header {
    sam::Header::builder()
        .add_reference_sequence(
            "chr1".to_owned(),
            Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()),
        )
        .build()
}

fn mapped_record(name: &str, start: u32, cigar: Cigar, mapq: u8, reverse: bool) -> RecordBuf {
    let mut flags = Flags::empty();
    if reverse {
        flags |= Flags::REVERSE_COMPLEMENTED;
    }

    let len = cigar.read_length();
    RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from((start + 1) as usize).unwrap())
        .set_mapping_quality(MappingQuality::try_from(mapq).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(vec![b'A'; len]))
        .set_quality_scores(QualityScores::from(vec![30; len]))
        .build()
}

fn match_cigar(len: usize) -> Cigar {
    vec![CigarOp::new(CigarKind::Match, len)]
        .into_iter()
        .collect()
}

fn make_executable(path: &std::path::Path) {
    #[cfg(unix)]
    {
        let mut perms = fs::metadata(path).unwrap().permissions();
        perms.set_mode(0o755);
        fs::set_permissions(path, perms).unwrap();
    }
}

#[test]
fn trackclustertu_run_executes_full_pipeline_with_fake_mapper_tools() {
    let tmp = unique_tmp_dir("trackclustertu_run_test");
    let tools_dir = tmp.join("bin");
    let out_dir = tmp.join("out");
    fs::create_dir_all(&tools_dir).unwrap();
    fs::create_dir_all(&out_dir).unwrap();

    let fixture_bam = tmp.join("fixture.bam");
    let header = build_header();
    let mut writer = bam::io::Writer::new(fs::File::create(&fixture_bam).unwrap());
    writer.write_header(&header).unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record("r1", 100, match_cigar(100), 60, false),
        )
        .unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record("r2", 101, match_cigar(100), 55, false),
        )
        .unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record("r3", 300, match_cigar(100), 50, false),
        )
        .unwrap();
    writer.try_finish().unwrap();

    let minimap2 = tools_dir.join("minimap2");
    fs::write(&minimap2, "#!/bin/sh\nprintf 'fake-sam\\n'\n").unwrap();
    make_executable(&minimap2);

    let samtools = tools_dir.join("samtools");
    fs::write(
        &samtools,
        r#"#!/bin/sh
cmd="$1"
shift
case "$cmd" in
  view)
    cat
    ;;
  sort)
    out=""
    while [ $# -gt 0 ]; do
      case "$1" in
        -o)
          out="$2"
          shift 2
          ;;
        -@)
          shift 2
          ;;
        -)
          shift
          ;;
        *)
          shift
          ;;
      esac
    done
    cat >/dev/null
    cp "$FAKE_BAM_PATH" "$out"
    ;;
  index)
    bam="$1"
    : > "${bam}.bai"
    ;;
  *)
    echo "unexpected samtools command: $cmd" >&2
    exit 1
    ;;
esac
"#,
    )
    .unwrap();
    make_executable(&samtools);

    let reference_fasta = tmp.join("ref.fa");
    fs::write(&reference_fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();

    let fastq = tmp.join("sampleA.fastq");
    fs::write(&fastq, "@r1\nACGTACGT\n+\nIIIIIIII\n").unwrap();

    let manifest = tmp.join("samples.fastq.tsv");
    fs::write(
        &manifest,
        format!(
            "sample\tgroup\treads\nsampleA\tcontrol\t{}\n",
            fastq.display()
        ),
    )
    .unwrap();

    let annotation_gff = tmp.join("genes.gff3");
    fs::write(
        &annotation_gff,
        concat!(
            "##gff-version 3\n",
            "chr1\tTest\tgene\t91\t150\t.\t+\t.\tID=id1;Name=geneA\n",
            "chr1\tTest\tgene\t161\t220\t.\t+\t.\tID=id2;Name=geneB\n",
            "chr1\tTest\tgene\t291\t310\t.\t+\t.\tID=id3;Name=geneC\n",
        ),
    )
    .unwrap();

    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let original_path = std::env::var_os("PATH").unwrap_or_default();
    let mut path_parts = vec![tools_dir.clone()];
    path_parts.extend(std::env::split_paths(&original_path));
    let joined_path = std::env::join_paths(path_parts).unwrap();

    let output = Command::new(exe)
        .env("PATH", joined_path)
        .env("FAKE_BAM_PATH", &fixture_bam)
        .args([
            "run",
            "--manifest",
            manifest.to_str().unwrap(),
            "--reference-fasta",
            reference_fasta.to_str().unwrap(),
            "--annotation-gff",
            annotation_gff.to_str().unwrap(),
            "--out-dir",
            out_dir.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let bam_manifest = fs::read_to_string(out_dir.join("samples.bam.tsv")).unwrap();
    assert!(bam_manifest.contains("sampleA"));
    assert!(bam_manifest.contains("sorted.bam"));

    assert_eq!(
        fs::read_to_string(out_dir.join("samples.bed.tsv")).unwrap(),
        format!(
            "sample\tgroup\treads\nsampleA\tcontrol\t{}\n",
            out_dir
                .join("bed/sampleA.bed")
                .canonicalize()
                .unwrap()
                .display()
        )
    );

    assert_eq!(
        fs::read_to_string(out_dir.join("tus.bed")).unwrap(),
        concat!(
            "chr1\t100\t200\tTU000001\t0\t+\n",
            "chr1\t300\t400\tTU000002\t0\t+\n",
        )
    );
    assert_eq!(
        fs::read_to_string(out_dir.join("tu_count.csv")).unwrap(),
        concat!("tu_id,count\n", "TU000001,2\n", "TU000002,1\n",)
    );
    assert_eq!(
        fs::read_to_string(out_dir.join("gene_count.csv")).unwrap(),
        concat!("gene_id,count\n", "geneA,2\n", "geneB,2\n", "geneC,1\n",)
    );

    let _ = fs::remove_dir_all(&tmp);
}

#[test]
fn trackclustertu_run_can_skip_score2_attachment() {
    let tmp = unique_tmp_dir("trackclustertu_run_skip_score2");
    let tools_dir = tmp.join("bin");
    let out_dir = tmp.join("out");
    fs::create_dir_all(&tools_dir).unwrap();
    fs::create_dir_all(&out_dir).unwrap();

    let fixture_bam = tmp.join("fixture.bam");
    let header = build_header();
    let mut writer = bam::io::Writer::new(fs::File::create(&fixture_bam).unwrap());
    writer.write_header(&header).unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record("r1", 100, match_cigar(100), 60, false),
        )
        .unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record("r2", 120, match_cigar(60), 55, false),
        )
        .unwrap();
    writer.try_finish().unwrap();

    let minimap2 = tools_dir.join("minimap2");
    fs::write(&minimap2, "#!/bin/sh\nprintf 'fake-sam\\n'\n").unwrap();
    make_executable(&minimap2);

    let samtools = tools_dir.join("samtools");
    fs::write(
        &samtools,
        r#"#!/bin/sh
cmd="$1"
shift
case "$cmd" in
  view)
    cat
    ;;
  sort)
    out=""
    while [ $# -gt 0 ]; do
      case "$1" in
        -o)
          out="$2"
          shift 2
          ;;
        -@)
          shift 2
          ;;
        -)
          shift
          ;;
        *)
          shift
          ;;
      esac
    done
    cat >/dev/null
    cp "$FAKE_BAM_PATH" "$out"
    ;;
  index)
    bam="$1"
    : > "${bam}.bai"
    ;;
  *)
    echo "unexpected samtools command: $cmd" >&2
    exit 1
    ;;
esac
"#,
    )
    .unwrap();
    make_executable(&samtools);

    let reference_fasta = tmp.join("ref.fa");
    fs::write(&reference_fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();

    let fastq = tmp.join("sampleA.fastq");
    fs::write(&fastq, "@r1\nACGTACGT\n+\nIIIIIIII\n").unwrap();

    let manifest = tmp.join("samples.fastq.tsv");
    fs::write(
        &manifest,
        format!(
            "sample\tgroup\treads\nsampleA\tcontrol\t{}\n",
            fastq.display()
        ),
    )
    .unwrap();

    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let original_path = std::env::var_os("PATH").unwrap_or_default();
    let mut path_parts = vec![tools_dir.clone()];
    path_parts.extend(std::env::split_paths(&original_path));
    let joined_path = std::env::join_paths(path_parts).unwrap();

    let output = Command::new(exe)
        .env("PATH", joined_path)
        .env("FAKE_BAM_PATH", &fixture_bam)
        .args([
            "run",
            "--manifest",
            manifest.to_str().unwrap(),
            "--reference-fasta",
            reference_fasta.to_str().unwrap(),
            "--out-dir",
            out_dir.to_str().unwrap(),
            "--skip-score2-attachment",
        ])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    assert_eq!(
        fs::read_to_string(out_dir.join("tus.bed")).unwrap(),
        concat!(
            "chr1\t100\t200\tTU000001\t0\t+\n",
            "chr1\t120\t180\tTU000002\t0\t+\n",
        )
    );
    assert_eq!(
        fs::read_to_string(out_dir.join("tu_count.csv")).unwrap(),
        concat!("tu_id,count\n", "TU000001,1\n", "TU000002,1\n",)
    );

    let _ = fs::remove_dir_all(&tmp);
}

#[test]
fn trackclustertu_run_forwards_three_prime_tolerance_to_cluster() {
    let tmp = unique_tmp_dir("trackclustertu_run_three_prime_tolerance");
    let tools_dir = tmp.join("bin");
    let out_dir = tmp.join("out");
    fs::create_dir_all(&tools_dir).unwrap();
    fs::create_dir_all(&out_dir).unwrap();

    let fixture_bam = tmp.join("fixture.bam");
    let header = build_header();
    let mut writer = bam::io::Writer::new(fs::File::create(&fixture_bam).unwrap());
    writer.write_header(&header).unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record("parent", 100, match_cigar(105), 60, false),
        )
        .unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record("child", 112, match_cigar(98), 55, false),
        )
        .unwrap();
    writer.try_finish().unwrap();

    let minimap2 = tools_dir.join("minimap2");
    fs::write(&minimap2, "#!/bin/sh\nprintf 'fake-sam\\n'\n").unwrap();
    make_executable(&minimap2);

    let samtools = tools_dir.join("samtools");
    fs::write(
        &samtools,
        r#"#!/bin/sh
cmd="$1"
shift
case "$cmd" in
  view)
    cat
    ;;
  sort)
    out=""
    while [ $# -gt 0 ]; do
      case "$1" in
        -o)
          out="$2"
          shift 2
          ;;
        -@)
          shift 2
          ;;
        -)
          shift
          ;;
        *)
          shift
          ;;
      esac
    done
    cat >/dev/null
    cp "$FAKE_BAM_PATH" "$out"
    ;;
  index)
    bam="$1"
    : > "${bam}.bai"
    ;;
  *)
    echo "unexpected samtools command: $cmd" >&2
    exit 1
    ;;
esac
"#,
    )
    .unwrap();
    make_executable(&samtools);

    let reference_fasta = tmp.join("ref.fa");
    fs::write(&reference_fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();

    let fastq = tmp.join("sampleA.fastq");
    fs::write(&fastq, "@r1\nACGTACGT\n+\nIIIIIIII\n").unwrap();

    let manifest = tmp.join("samples.fastq.tsv");
    fs::write(
        &manifest,
        format!(
            "sample\tgroup\treads\nsampleA\tcontrol\t{}\n",
            fastq.display()
        ),
    )
    .unwrap();

    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let original_path = std::env::var_os("PATH").unwrap_or_default();
    let mut path_parts = vec![tools_dir.clone()];
    path_parts.extend(std::env::split_paths(&original_path));
    let joined_path = std::env::join_paths(path_parts).unwrap();

    let output = Command::new(exe)
        .env("PATH", joined_path)
        .env("FAKE_BAM_PATH", &fixture_bam)
        .args([
            "run",
            "--manifest",
            manifest.to_str().unwrap(),
            "--reference-fasta",
            reference_fasta.to_str().unwrap(),
            "--out-dir",
            out_dir.to_str().unwrap(),
            "--three-prime-tolerance-bp",
            "0",
        ])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    assert_eq!(
        fs::read_to_string(out_dir.join("tus.bed")).unwrap(),
        concat!(
            "chr1\t100\t205\tTU000001\t0\t+\n",
            "chr1\t112\t210\tTU000002\t0\t+\n",
        )
    );
    assert_eq!(
        fs::read_to_string(out_dir.join("tu_count.csv")).unwrap(),
        concat!("tu_id,count\n", "TU000001,1\n", "TU000002,1\n",)
    );

    let _ = fs::remove_dir_all(&tmp);
}

#[test]
fn trackclustertu_run_forwards_five_prime_delta_to_cluster() {
    let tmp = unique_tmp_dir("trackclustertu_run_five_prime_delta");
    let tools_dir = tmp.join("bin");
    let out_dir = tmp.join("out");
    fs::create_dir_all(&tools_dir).unwrap();
    fs::create_dir_all(&out_dir).unwrap();

    let fixture_bam = tmp.join("fixture.bam");
    let header = build_header();
    let mut writer = bam::io::Writer::new(fs::File::create(&fixture_bam).unwrap());
    writer.write_header(&header).unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record("parent", 100, match_cigar(105), 60, false),
        )
        .unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record("child", 150, match_cigar(60), 55, false),
        )
        .unwrap();
    writer.try_finish().unwrap();

    let minimap2 = tools_dir.join("minimap2");
    fs::write(&minimap2, "#!/bin/sh\nprintf 'fake-sam\\n'\n").unwrap();
    make_executable(&minimap2);

    let samtools = tools_dir.join("samtools");
    fs::write(
        &samtools,
        r#"#!/bin/sh
cmd="$1"
shift
case "$cmd" in
  view)
    cat
    ;;
  sort)
    out=""
    while [ $# -gt 0 ]; do
      case "$1" in
        -o)
          out="$2"
          shift 2
          ;;
        -@)
          shift 2
          ;;
        -)
          shift
          ;;
        *)
          shift
          ;;
      esac
    done
    cat >/dev/null
    cp "$FAKE_BAM_PATH" "$out"
    ;;
  index)
    bam="$1"
    : > "${bam}.bai"
    ;;
  *)
    echo "unexpected samtools command: $cmd" >&2
    exit 1
    ;;
esac
"#,
    )
    .unwrap();
    make_executable(&samtools);

    let reference_fasta = tmp.join("ref.fa");
    fs::write(&reference_fasta, ">chr1\nACGTACGTACGTACGT\n").unwrap();

    let fastq = tmp.join("sampleA.fastq");
    fs::write(&fastq, "@r1\nACGTACGT\n+\nIIIIIIII\n").unwrap();

    let manifest = tmp.join("samples.fastq.tsv");
    fs::write(
        &manifest,
        format!(
            "sample\tgroup\treads\nsampleA\tcontrol\t{}\n",
            fastq.display()
        ),
    )
    .unwrap();

    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let original_path = std::env::var_os("PATH").unwrap_or_default();
    let mut path_parts = vec![tools_dir.clone()];
    path_parts.extend(std::env::split_paths(&original_path));
    let joined_path = std::env::join_paths(path_parts).unwrap();

    let output = Command::new(exe)
        .env("PATH", joined_path)
        .env("FAKE_BAM_PATH", &fixture_bam)
        .args([
            "run",
            "--manifest",
            manifest.to_str().unwrap(),
            "--reference-fasta",
            reference_fasta.to_str().unwrap(),
            "--out-dir",
            out_dir.to_str().unwrap(),
            "--max-5p-delta",
            "50",
        ])
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    assert_eq!(
        fs::read_to_string(out_dir.join("tus.bed")).unwrap(),
        "chr1\t100\t205\tTU000001\t0\t+\n"
    );
    assert_eq!(
        fs::read_to_string(out_dir.join("tu_count.csv")).unwrap(),
        concat!("tu_id,count\n", "TU000001,2\n",)
    );

    let _ = fs::remove_dir_all(&tmp);
}
