use std::fs;
use std::num::NonZeroUsize;
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
        .add_reference_sequence(
            "chr2".to_owned(),
            Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()),
        )
        .build()
}

fn mapped_record(
    name: &str,
    reference_sequence_id: usize,
    start: u32,
    cigar: Cigar,
    mapq: u8,
    reverse: bool,
    secondary: bool,
) -> RecordBuf {
    let mut flags = Flags::empty();
    if reverse {
        flags |= Flags::REVERSE_COMPLEMENTED;
    }
    if secondary {
        flags |= Flags::SECONDARY;
    }

    let len = cigar.read_length();
    RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_reference_sequence_id(reference_sequence_id)
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

fn skipped_cigar(left: usize, skip: usize, right: usize) -> Cigar {
    vec![
        CigarOp::new(CigarKind::Match, left),
        CigarOp::new(CigarKind::Skip, skip),
        CigarOp::new(CigarKind::Match, right),
    ]
    .into_iter()
    .collect()
}

fn unmapped_record(name: &str, len: usize) -> RecordBuf {
    RecordBuf::builder()
        .set_name(name)
        .set_flags(Flags::UNMAPPED)
        .set_sequence(Sequence::from(vec![b'A'; len]))
        .set_quality_scores(QualityScores::from(vec![30; len]))
        .build()
}

#[test]
fn trackclustertu_bam_to_bed_subcommand_converts_primary_alignments_to_bed6() {
    let tmp = unique_tmp_dir("trackclustertu_bam_to_bed_test");
    fs::create_dir_all(&tmp).unwrap();

    let bam_path = tmp.join("reads.bam");
    let bed_path = tmp.join("reads.bed");
    let header = build_header();

    let mut writer = bam::io::Writer::new(fs::File::create(&bam_path).unwrap());
    writer.write_header(&header).unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record("r_plus", 0, 99, match_cigar(10), 42, false, false),
        )
        .unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record("r_minus", 0, 149, match_cigar(6), 60, true, false),
        )
        .unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record("r_secondary", 0, 199, match_cigar(8), 50, false, true),
        )
        .unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record("r_chr2", 1, 49, match_cigar(5), 31, false, false),
        )
        .unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record(
                "r_spliced",
                0,
                299,
                skipped_cigar(10, 100, 10),
                55,
                false,
                false,
            ),
        )
        .unwrap();
    writer
        .write_alignment_record(&header, &unmapped_record("r_unmapped", 7))
        .unwrap();
    writer.try_finish().unwrap();

    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let output = Command::new(exe)
        .args([
            "bam-to-bed",
            "--in-bam",
            bam_path.to_str().unwrap(),
            "--out-bed",
            bed_path.to_str().unwrap(),
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
        fs::read_to_string(&bed_path).unwrap(),
        concat!(
            "chr1\t99\t109\tr_plus\t42\t+\n",
            "chr1\t149\t155\tr_minus\t60\t-\n",
            "chr2\t49\t54\tr_chr2\t31\t+\n",
        )
    );

    let _ = fs::remove_dir_all(&tmp);
}

#[test]
fn trackclustertu_bam_to_bed_manifest_rejects_sanitized_name_collisions() {
    let tmp = unique_tmp_dir("trackclustertu_bam_to_bed_collision_test");
    fs::create_dir_all(&tmp).unwrap();

    let manifest_path = tmp.join("samples.tsv");
    fs::write(
        &manifest_path,
        concat!(
            "sample\treads\n",
            "sample a\tmissing1.bam\n",
            "sample/a\tmissing2.bam\n",
        ),
    )
    .unwrap();

    let out_dir = tmp.join("out");
    let exe = env!("CARGO_BIN_EXE_trackclustertu");
    let output = Command::new(exe)
        .args([
            "bam-to-bed",
            "--manifest",
            manifest_path.to_str().unwrap(),
            "--out-dir",
            out_dir.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(
        !output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("filename sanitization") && stderr.contains("sample a"),
        "stderr did not explain sample-name collision:\n{stderr}",
    );

    let _ = fs::remove_dir_all(&tmp);
}
