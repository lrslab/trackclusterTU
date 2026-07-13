use std::fs;
use std::io::{BufRead, BufReader, Read, Write};
use std::num::NonZeroUsize;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use noodles::sam::alignment::io::Write as _;
use noodles::{
    bam, bgzf,
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
        header::record::value::{
            map::{self, header::tag as header_tag, ReferenceSequence},
            Map,
        },
    },
};
use trackclustertu::bam::{
    bam_to_bed6_with_evidence, BamConversionConfig, BamConversionError, BamConversionSummary,
    LibraryProfile,
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

fn build_coordinate_header() -> sam::Header {
    let mut header = build_header();
    let mut header_record = Map::<map::Header>::default();
    header_record.other_fields_mut().insert(
        header_tag::SORT_ORDER,
        map::header::sort_order::COORDINATE.into(),
    );
    *header.header_mut() = Some(header_record);
    header
}

fn corrupt_record_name_length(path: &std::path::Path, record_index: usize) {
    fn read_u32_le(src: &[u8], offset: usize) -> usize {
        usize::try_from(u32::from_le_bytes(
            src[offset..offset + 4].try_into().unwrap(),
        ))
        .unwrap()
    }

    let mut raw = Vec::new();
    bgzf::io::Reader::new(fs::File::open(path).unwrap())
        .read_to_end(&mut raw)
        .unwrap();
    assert_eq!(&raw[..4], b"BAM\x01");

    let mut offset = 4;
    let header_len = read_u32_le(&raw, offset);
    offset += 4 + header_len;
    let reference_count = read_u32_le(&raw, offset);
    offset += 4;
    for _ in 0..reference_count {
        let name_len = read_u32_le(&raw, offset);
        offset += 4 + name_len + 4;
    }

    for current_index in 0..=record_index {
        let block_size = read_u32_le(&raw, offset);
        let body_start = offset + 4;
        assert!(body_start + block_size <= raw.len());
        if current_index == record_index {
            // l_read_name is byte 8 of the BAM alignment payload. Zero is invalid, but the
            // surrounding length-prefixed record block remains intact and recoverable.
            raw[body_start + 8] = 0;
            break;
        }
        offset = body_start + block_size;
    }

    let mut writer = bgzf::io::Writer::new(fs::File::create(path).unwrap());
    writer.write_all(&raw).unwrap();
    writer.try_finish().unwrap();
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

fn soft_clipped_cigar(leading: usize, matched: usize, trailing: usize) -> Cigar {
    vec![
        CigarOp::new(CigarKind::HardClip, 3),
        CigarOp::new(CigarKind::SoftClip, leading),
        CigarOp::new(CigarKind::Match, matched),
        CigarOp::new(CigarKind::SoftClip, trailing),
        CigarOp::new(CigarKind::HardClip, 4),
    ]
    .into_iter()
    .collect()
}

fn mapped_record_with_flags(
    name: Option<&str>,
    start: u32,
    cigar: Cigar,
    mapq: u8,
    flags: Flags,
) -> RecordBuf {
    let len = cigar.read_length();
    let mut builder = RecordBuf::builder();
    if let Some(name) = name {
        builder = builder.set_name(name);
    }
    builder
        .set_flags(flags)
        .set_reference_sequence_id(0)
        .set_alignment_start(Position::try_from((start + 1) as usize).unwrap())
        .set_mapping_quality(MappingQuality::try_from(mapq).unwrap())
        .set_cigar(cigar)
        .set_sequence(Sequence::from(vec![b'A'; len]))
        .set_quality_scores(QualityScores::from(vec![30; len]))
        .build()
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
fn bam_library_api_reports_typed_config_and_evidence_errors() {
    let tmp = unique_tmp_dir("trackclustertu_bam_typed_input_errors");
    fs::create_dir_all(&tmp).unwrap();
    let missing_bam = tmp.join("missing.bam");
    let out_bed = tmp.join("reads.bed");

    let error = bam_to_bed6_with_evidence(
        &missing_bam,
        &out_bed,
        None,
        None,
        &BamConversionConfig {
            min_boundary_support: 0,
            ..BamConversionConfig::default()
        },
    )
    .unwrap_err();
    assert!(matches!(
        error,
        BamConversionError::InvalidMinimumBoundarySupport { value: 0 }
    ));

    let evidence_path = tmp.join("input.evidence.tsv");
    fs::write(
        &evidence_path,
        "read_name\tpoly_a_evidence\nread-1\tmaybe\n",
    )
    .unwrap();
    let error = bam_to_bed6_with_evidence(
        &missing_bam,
        &out_bed,
        None,
        Some(&evidence_path),
        &BamConversionConfig::default(),
    )
    .unwrap_err();
    match error {
        BamConversionError::EvidenceInvalidValue {
            path,
            line,
            column,
            value,
        } => {
            assert_eq!(path, evidence_path);
            assert_eq!(line, 2);
            assert_eq!(column, "poly_a_evidence");
            assert_eq!(value, "maybe");
        }
        other => panic!("unexpected error: {other:?}"),
    }

    let missing_evidence = tmp.join("missing.evidence.tsv");
    let error = bam_to_bed6_with_evidence(
        &missing_bam,
        &out_bed,
        None,
        Some(&missing_evidence),
        &BamConversionConfig::default(),
    )
    .unwrap_err();
    match error {
        BamConversionError::EvidenceOpen { path, source } => {
            assert_eq!(path, missing_evidence);
            assert_eq!(source.kind(), std::io::ErrorKind::NotFound);
        }
        other => panic!("unexpected error: {other:?}"),
    }

    let _ = fs::remove_dir_all(tmp);
}

#[test]
fn bam_library_api_preserves_bam_open_and_output_context() {
    let tmp = unique_tmp_dir("trackclustertu_bam_typed_io_errors");
    fs::create_dir_all(&tmp).unwrap();
    let missing_bam = tmp.join("missing.bam");
    let out_bed = tmp.join("reads.bed");

    let error = bam_to_bed6_with_evidence(
        &missing_bam,
        &out_bed,
        None,
        None,
        &BamConversionConfig::default(),
    )
    .unwrap_err();
    match error {
        BamConversionError::BamOpen {
            path,
            context,
            source,
        } => {
            assert_eq!(path, missing_bam);
            assert_eq!(context, "inspecting sort order");
            assert_eq!(source.kind(), std::io::ErrorKind::NotFound);
        }
        other => panic!("unexpected error: {other:?}"),
    }

    let bam_path = tmp.join("reads.bam");
    let header = build_header();
    let mut writer = bam::io::Writer::new(fs::File::create(&bam_path).unwrap());
    writer.write_header(&header).unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record("r1", 0, 10, match_cigar(10), 42, false, false),
        )
        .unwrap();
    writer.try_finish().unwrap();
    drop(writer);

    let missing_parent_output = tmp.join("missing-parent/reads.bed");
    let error = bam_to_bed6_with_evidence(
        &bam_path,
        &missing_parent_output,
        None,
        None,
        &BamConversionConfig::default(),
    )
    .unwrap_err();
    match error {
        BamConversionError::OutputCreate {
            path,
            output,
            source,
        } => {
            assert_eq!(path, missing_parent_output);
            assert_eq!(output, "BED");
            assert_eq!(source.kind(), std::io::ErrorKind::NotFound);
        }
        other => panic!("unexpected error: {other:?}"),
    }

    let _ = fs::remove_dir_all(tmp);
}

#[test]
fn bam_library_api_distinguishes_header_and_record_decode_failures() {
    let tmp = unique_tmp_dir("trackclustertu_bam_typed_decode_errors");
    fs::create_dir_all(&tmp).unwrap();
    let out_bed = tmp.join("reads.bed");

    let invalid_header = tmp.join("invalid-header.bam");
    fs::write(&invalid_header, b"not-a-bgzf-bam").unwrap();
    let error = bam_to_bed6_with_evidence(
        &invalid_header,
        &out_bed,
        None,
        None,
        &BamConversionConfig::default(),
    )
    .unwrap_err();
    match error {
        BamConversionError::BamHeaderRead {
            path,
            context,
            source,
        } => {
            assert_eq!(path, invalid_header);
            assert_eq!(context, "inspecting sort order");
            assert_ne!(source.kind(), std::io::ErrorKind::NotFound);
        }
        other => panic!("unexpected error: {other:?}"),
    }

    let invalid_record = tmp.join("invalid-record.bam");
    let header = build_header();
    let mut writer = bam::io::Writer::new(fs::File::create(&invalid_record).unwrap());
    writer.write_header(&header).unwrap();
    // Force the complete header into its own BGZF block so truncating the next
    // block exercises record decoding rather than header decoding.
    std::io::Write::flush(writer.get_mut()).unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record("corrupt-me", 0, 10, match_cigar(10), 42, false, false),
        )
        .unwrap();
    std::io::Write::flush(writer.get_mut()).unwrap();
    writer.try_finish().unwrap();
    drop(writer);
    let mut bytes = fs::read(&invalid_record).unwrap();
    // Remove the canonical 28-byte BGZF EOF marker and truncate the separately
    // flushed record block itself.
    bytes.truncate(bytes.len() - 28 - 4);
    fs::write(&invalid_record, bytes).unwrap();

    let error = bam_to_bed6_with_evidence(
        &invalid_record,
        &out_bed,
        None,
        None,
        &BamConversionConfig::default(),
    )
    .unwrap_err();
    match error {
        BamConversionError::BamRecordRead {
            path,
            record_ordinal,
            context,
            source,
        } => {
            assert_eq!(path, invalid_record);
            // BGZF decompression can yield the first record before discovering
            // the truncated block footer on the next read attempt.
            assert_eq!(record_ordinal, 2);
            assert_eq!(
                context,
                "buffering a BAM without a coordinate-sort declaration"
            );
            assert_ne!(source.kind(), std::io::ErrorKind::NotFound);
        }
        other => panic!("unexpected error: {other:?}"),
    }

    let _ = fs::remove_dir_all(tmp);
}

#[test]
fn intact_malformed_record_is_skipped_between_valid_records() {
    let tmp = unique_tmp_dir("trackclustertu_bam_recoverable_record_error");
    fs::create_dir_all(&tmp).unwrap();
    let bam_path = tmp.join("reads.bam");
    let bed_path = tmp.join("reads.bed");
    let evidence_path = tmp.join("reads.evidence.tsv");
    let header = build_coordinate_header();

    let mut writer = bam::io::Writer::new(fs::File::create(&bam_path).unwrap());
    writer.write_header(&header).unwrap();
    for (name, start) in [("good-before", 10), ("bad-middle", 20), ("good-after", 30)] {
        writer
            .write_alignment_record(
                &header,
                &mapped_record(name, 0, start, match_cigar(5), 42, false, false),
            )
            .unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);
    corrupt_record_name_length(&bam_path, 1);

    let summary = bam_to_bed6_with_evidence(
        &bam_path,
        &bed_path,
        Some(&evidence_path),
        None,
        &BamConversionConfig::default(),
    )
    .unwrap();

    assert_eq!(
        summary,
        BamConversionSummary {
            total: 3,
            retained: 2,
            malformed_record: 1,
            ..BamConversionSummary::default()
        }
    );
    assert_eq!(
        fs::read_to_string(&bed_path).unwrap(),
        concat!(
            "chr1\t10\t15\tgood-before\t42\t+\n",
            "chr1\t30\t35\tgood-after\t42\t+\n",
        )
    );

    let evidence = fs::read_to_string(&evidence_path).unwrap();
    let rows: Vec<Vec<&str>> = evidence
        .lines()
        .skip(2)
        .map(|line| line.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), 3);
    assert_eq!(rows[0][1], "1");
    assert_eq!(rows[0][2], "good-before");
    assert_eq!(rows[1][1], "2");
    assert_eq!(rows[1][2], ".");
    assert_eq!(rows[1][12], "malformed_record");
    assert_eq!(rows[2][1], "3");
    assert_eq!(rows[2][2], "good-after");

    let _ = fs::remove_dir_all(tmp);
}

#[test]
fn two_pass_support_skips_malformed_record_once_and_keeps_later_reads() {
    let tmp = unique_tmp_dir("trackclustertu_bam_two_pass_recoverable_record_error");
    fs::create_dir_all(&tmp).unwrap();
    let bam_path = tmp.join("reads.bam");
    let bed_path = tmp.join("reads.bed");
    let evidence_path = tmp.join("reads.evidence.tsv");
    let header = build_coordinate_header();

    let mut writer = bam::io::Writer::new(fs::File::create(&bam_path).unwrap());
    writer.write_header(&header).unwrap();
    for (name, start) in [
        ("a1", 10),
        ("a2", 10),
        ("bad-middle", 20),
        ("b1", 30),
        ("b2", 30),
    ] {
        writer
            .write_alignment_record(
                &header,
                &mapped_record(name, 0, start, match_cigar(5), 42, false, false),
            )
            .unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);
    corrupt_record_name_length(&bam_path, 2);

    let summary = bam_to_bed6_with_evidence(
        &bam_path,
        &bed_path,
        Some(&evidence_path),
        None,
        &BamConversionConfig {
            min_boundary_support: 2,
            ..BamConversionConfig::default()
        },
    )
    .unwrap();

    assert_eq!(
        summary,
        BamConversionSummary {
            total: 5,
            retained: 4,
            malformed_record: 1,
            ..BamConversionSummary::default()
        }
    );
    assert_eq!(
        fs::read_to_string(&bed_path).unwrap(),
        concat!(
            "chr1\t10\t15\ta1\t42\t+\n",
            "chr1\t10\t15\ta2\t42\t+\n",
            "chr1\t30\t35\tb1\t42\t+\n",
            "chr1\t30\t35\tb2\t42\t+\n",
        )
    );

    let evidence = fs::read_to_string(&evidence_path).unwrap();
    let rows: Vec<Vec<&str>> = evidence
        .lines()
        .skip(2)
        .map(|line| line.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), 5);
    assert_eq!(rows[2][1], "3");
    assert_eq!(rows[2][2], ".");
    assert_eq!(rows[2][12], "malformed_record");
    assert_eq!(rows[3][2], "b1");
    assert_eq!(rows[3][12], "retained");
    assert_eq!(rows[4][2], "b2");
    assert_eq!(rows[4][12], "retained");

    let _ = fs::remove_dir_all(tmp);
}

#[test]
fn trackclustertu_bam_to_bed_subcommand_converts_primary_alignments_to_bed6() {
    let tmp = unique_tmp_dir("trackclustertu_bam_to_bed_test");
    fs::create_dir_all(&tmp).unwrap();

    let bam_path = tmp.join("reads.bam");
    let bed_path = tmp.join("nested/bed/reads.bed");
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
fn trackclustertu_bam_to_bed_rejects_output_aliasing_input_without_truncation() {
    let tmp = unique_tmp_dir("trackclustertu_bam_alias_test");
    fs::create_dir_all(&tmp).unwrap();
    let bam_path = tmp.join("reads.bam");
    let header = build_header();
    let mut writer = bam::io::Writer::new(fs::File::create(&bam_path).unwrap());
    writer.write_header(&header).unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record("r1", 0, 10, match_cigar(10), 42, false, false),
        )
        .unwrap();
    writer.try_finish().unwrap();
    drop(writer);
    let original = fs::read(&bam_path).unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_trackclustertu"))
        .args([
            "bam-to-bed",
            "--in-bam",
            bam_path.to_str().unwrap(),
            "--out-bed",
            bam_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("aliases input"),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(fs::read(&bam_path).unwrap(), original);
    let _ = fs::remove_dir_all(tmp);
}

#[test]
fn coordinate_sorted_bam_streams_and_validates_record_order() {
    let tmp = unique_tmp_dir("trackclustertu_coordinate_sorted_bam");
    fs::create_dir_all(&tmp).unwrap();
    let bam_path = tmp.join("reads.bam");
    let bed_path = tmp.join("reads.bed");
    let evidence_path = tmp.join("reads.evidence.tsv");
    let header = build_coordinate_header();
    let mut writer = bam::io::Writer::new(fs::File::create(&bam_path).unwrap());
    writer.write_header(&header).unwrap();
    for (name, start) in [("r1", 10), ("r2", 20), ("r3", 30)] {
        writer
            .write_alignment_record(
                &header,
                &mapped_record(name, 0, start, match_cigar(5), 42, false, false),
            )
            .unwrap();
    }
    writer.try_finish().unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_trackclustertu"))
        .args([
            "bam-to-bed",
            "--in-bam",
            bam_path.to_str().unwrap(),
            "--out-bed",
            bed_path.to_str().unwrap(),
            "--out-evidence",
            evidence_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(fs::read_to_string(&bed_path).unwrap().lines().count(), 3);
    let evidence = fs::read_to_string(&evidence_path).unwrap();
    assert_eq!(evidence.lines().count(), 5);
    assert_eq!(
        evidence
            .lines()
            .skip(2)
            .map(|line| line.split('\t').nth(2).unwrap())
            .collect::<Vec<_>>(),
        ["r1", "r2", "r3"]
    );
    let _ = fs::remove_dir_all(tmp);
}

#[test]
fn coordinate_sorted_boundary_support_uses_two_pass_evidence_and_single_counts() {
    let tmp = unique_tmp_dir("trackclustertu_coordinate_two_pass");
    fs::create_dir_all(&tmp).unwrap();
    let bam_path = tmp.join("reads.bam");
    let bed_path = tmp.join("reads.bed");
    let input_evidence_path = tmp.join("input.evidence.tsv");
    let output_evidence_path = tmp.join("output.evidence.tsv");
    let header = build_coordinate_header();
    let mut writer = bam::io::Writer::new(fs::File::create(&bam_path).unwrap());
    writer.write_header(&header).unwrap();
    for (name, start, mapq) in [
        ("supported_z", 99, 30),
        ("supported_a", 99, 31),
        ("singleton", 199, 32),
        ("low_mapq", 299, 5),
    ] {
        writer
            .write_alignment_record(
                &header,
                &mapped_record(name, 0, start, match_cigar(10), mapq, false, false),
            )
            .unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);
    fs::write(
        &input_evidence_path,
        concat!(
            "read_name\tfive_prime_adapter_evidence\tthree_prime_adapter_evidence\n",
            "supported_z\tpresent\tpresent\n",
            "supported_a\tpresent\tpresent\n",
            "singleton\tpresent\tpresent\n",
            "low_mapq\tpresent\tpresent\n",
        ),
    )
    .unwrap();

    let summary = bam_to_bed6_with_evidence(
        &bam_path,
        &bed_path,
        Some(&output_evidence_path),
        Some(&input_evidence_path),
        &BamConversionConfig {
            min_mapq: 20,
            require_full_length: true,
            min_boundary_support: 2,
            library_profile: LibraryProfile::DirectCdna,
        },
    )
    .unwrap();

    assert_eq!(
        summary,
        BamConversionSummary {
            total: 4,
            retained: 2,
            below_min_mapq: 1,
            below_min_boundary_support: 1,
            ..BamConversionSummary::default()
        }
    );
    assert_eq!(
        fs::read_to_string(&bed_path).unwrap(),
        concat!(
            "chr1\t99\t109\tsupported_z\t30\t+\n",
            "chr1\t99\t109\tsupported_a\t31\t+\n",
        )
    );
    let evidence = fs::read_to_string(&output_evidence_path).unwrap();
    let rows: Vec<Vec<&str>> = evidence
        .lines()
        .skip(2)
        .map(|line| line.split('\t').collect())
        .collect();
    assert_eq!(
        rows.iter().map(|row| row[2]).collect::<Vec<_>>(),
        ["supported_z", "supported_a", "singleton", "low_mapq"]
    );
    assert_eq!(rows[0][12], "retained");
    assert_eq!(rows[1][12], "retained");
    assert_eq!(rows[2][12], "below_min_boundary_support");
    assert_eq!(rows[3][12], "below_min_mapq");

    let bed_without_sidecar = tmp.join("reads.without-sidecar.bed");
    let summary_without_sidecar = bam_to_bed6_with_evidence(
        &bam_path,
        &bed_without_sidecar,
        None,
        Some(&input_evidence_path),
        &BamConversionConfig {
            min_mapq: 20,
            require_full_length: true,
            min_boundary_support: 2,
            library_profile: LibraryProfile::DirectCdna,
        },
    )
    .unwrap();
    assert_eq!(summary_without_sidecar, summary);
    assert_eq!(
        fs::read_to_string(bed_without_sidecar).unwrap(),
        fs::read_to_string(bed_path).unwrap()
    );
    let _ = fs::remove_dir_all(tmp);
}

#[test]
fn large_coordinate_sorted_evidence_and_support_smoke_test() {
    const BOUNDARY_COUNT: usize = 16;
    const RECORDS_PER_BOUNDARY: usize = 1024;
    const TOTAL_RECORDS: usize = BOUNDARY_COUNT * RECORDS_PER_BOUNDARY;

    let tmp = unique_tmp_dir("trackclustertu_large_coordinate_two_pass");
    fs::create_dir_all(&tmp).unwrap();
    let bam_path = tmp.join("reads.bam");
    let bed_path = tmp.join("reads.bed");
    let evidence_path = tmp.join("reads.evidence.tsv");
    let header = build_coordinate_header();
    let mut writer = bam::io::Writer::new(fs::File::create(&bam_path).unwrap());
    writer.write_header(&header).unwrap();
    for boundary_idx in 0..BOUNDARY_COUNT {
        let start = u32::try_from(boundary_idx * 20).unwrap();
        for read_idx in 0..RECORDS_PER_BOUNDARY {
            let name = format!("b{boundary_idx:02}_r{read_idx:04}");
            writer
                .write_alignment_record(
                    &header,
                    &mapped_record(&name, 0, start, match_cigar(10), 30, false, false),
                )
                .unwrap();
        }
    }
    writer.try_finish().unwrap();
    drop(writer);

    let summary = bam_to_bed6_with_evidence(
        &bam_path,
        &bed_path,
        Some(&evidence_path),
        None,
        &BamConversionConfig {
            min_boundary_support: 2,
            ..BamConversionConfig::default()
        },
    )
    .unwrap();
    assert_eq!(summary.total, TOTAL_RECORDS as u64);
    assert_eq!(summary.retained, TOTAL_RECORDS as u64);
    assert_eq!(summary.below_min_boundary_support, 0);
    assert_eq!(
        BufReader::new(fs::File::open(&bed_path).unwrap())
            .lines()
            .count(),
        TOTAL_RECORDS
    );
    assert_eq!(
        BufReader::new(fs::File::open(&evidence_path).unwrap())
            .lines()
            .count(),
        TOTAL_RECORDS + 2
    );
    let _ = fs::remove_dir_all(tmp);
}

#[test]
fn falsely_coordinate_sorted_bam_fails_without_replacing_prior_output() {
    let tmp = unique_tmp_dir("trackclustertu_invalid_coordinate_order");
    fs::create_dir_all(&tmp).unwrap();
    let bam_path = tmp.join("reads.bam");
    let bed_path = tmp.join("reads.bed");
    let evidence_path = tmp.join("reads.evidence.tsv");
    let header = build_coordinate_header();
    let mut writer = bam::io::Writer::new(fs::File::create(&bam_path).unwrap());
    writer.write_header(&header).unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record("later", 0, 100, match_cigar(5), 42, false, false),
        )
        .unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record("earlier", 0, 50, match_cigar(5), 42, false, false),
        )
        .unwrap();
    writer.try_finish().unwrap();
    drop(writer);

    let raw_bed_path = tmp.join("raw.bed");
    let error = bam_to_bed6_with_evidence(
        &bam_path,
        &raw_bed_path,
        None,
        None,
        &BamConversionConfig {
            min_boundary_support: 2,
            ..BamConversionConfig::default()
        },
    )
    .unwrap_err();
    match error {
        BamConversionError::OutOfOrder {
            path,
            record_ordinal,
            current_reference_id,
            current_position,
            previous_reference_id,
            previous_position,
        } => {
            assert_eq!(path, bam_path);
            assert_eq!(record_ordinal, 2);
            assert_eq!(current_reference_id, 0);
            assert_eq!(current_position, 51);
            assert_eq!(previous_reference_id, 0);
            assert_eq!(previous_position, 101);
        }
        other => panic!("unexpected error: {other:?}"),
    }
    assert!(!raw_bed_path.exists());

    fs::write(&bed_path, "previous-success\n").unwrap();
    fs::write(&evidence_path, "previous-evidence\n").unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_trackclustertu"))
        .args([
            "bam-to-bed",
            "--in-bam",
            bam_path.to_str().unwrap(),
            "--out-bed",
            bed_path.to_str().unwrap(),
            "--out-evidence",
            evidence_path.to_str().unwrap(),
            "--min-boundary-support",
            "2",
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("out of order"),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(fs::read_to_string(&bed_path).unwrap(), "previous-success\n");
    assert_eq!(
        fs::read_to_string(&evidence_path).unwrap(),
        "previous-evidence\n"
    );
    let _ = fs::remove_dir_all(tmp);
}

#[test]
fn trackclustertu_bam_to_bed_emits_versioned_evidence_and_deterministic_counters() {
    let tmp = unique_tmp_dir("trackclustertu_bam_evidence_test");
    fs::create_dir_all(&tmp).unwrap();
    let bam_path = tmp.join("reads.bam");
    let bed_path = tmp.join("reads.bed");
    let input_evidence_path = tmp.join("input.evidence.tsv");
    let output_evidence_path = tmp.join("output.evidence.tsv");
    let header = build_header();

    let mut writer = bam::io::Writer::new(fs::File::create(&bam_path).unwrap());
    writer.write_header(&header).unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record_with_flags(
                Some("full_plus"),
                99,
                soft_clipped_cigar(5, 10, 2),
                42,
                Flags::empty(),
            ),
        )
        .unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record_with_flags(
                Some("full_minus"),
                149,
                soft_clipped_cigar(5, 10, 2),
                43,
                Flags::REVERSE_COMPLEMENTED,
            ),
        )
        .unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record_with_flags(Some("low_mapq"), 199, match_cigar(10), 5, Flags::empty()),
        )
        .unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record_with_flags(
                Some("secondary"),
                249,
                match_cigar(10),
                50,
                Flags::SECONDARY,
            ),
        )
        .unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record_with_flags(
                Some("supplementary"),
                299,
                match_cigar(10),
                50,
                Flags::SUPPLEMENTARY,
            ),
        )
        .unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record_with_flags(
                Some("spliced"),
                349,
                skipped_cigar(5, 20, 5),
                50,
                Flags::empty(),
            ),
        )
        .unwrap();
    writer
        .write_alignment_record(&header, &unmapped_record("unmapped", 7))
        .unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record_with_flags(None, 399, match_cigar(10), 50, Flags::empty()),
        )
        .unwrap();
    writer.try_finish().unwrap();

    fs::write(
        &input_evidence_path,
        concat!(
            "read_name\tpoly_a_evidence\tfive_prime_adapter_evidence\t",
            "three_prime_adapter_evidence\tfull_length_evidence\tlibrary_preparation\n",
            "full_plus\tpresent\tpresent\tunknown\tunknown\tdirect-rna-kit\n",
            "full_minus\tyes\tyes\t.\t.\tdirect-rna-kit\n",
            "low_mapq\tpresent\tpresent\tunknown\tunknown\tdirect-rna-kit\n",
        ),
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_trackclustertu"))
        .args([
            "bam-to-bed",
            "--in-bam",
            bam_path.to_str().unwrap(),
            "--out-bed",
            bed_path.to_str().unwrap(),
            "--in-evidence",
            input_evidence_path.to_str().unwrap(),
            "--out-evidence",
            output_evidence_path.to_str().unwrap(),
            "--library-profile",
            "direct-rna",
            "--require-full-length",
            "--min-mapq",
            "20",
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
            "chr1\t99\t109\tfull_plus\t42\t+\n",
            "chr1\t149\t159\tfull_minus\t43\t-\n",
        )
    );
    let evidence = fs::read_to_string(&output_evidence_path).unwrap();
    let mut lines = evidence.lines();
    assert_eq!(
        lines.next(),
        Some("#schema_version=trackclustertu.bam-evidence.v1")
    );
    let header = lines.next().unwrap();
    assert!(header.starts_with("schema_version\trecord_ordinal\tread_name"));
    let rows: Vec<&str> = lines.collect();
    assert_eq!(rows.len(), 8);
    assert!(rows[0].contains(
        "trackclustertu.bam-evidence.v1\t1\tfull_plus\tchr1\t99\t109\t+\t42\t5\t2\t3H5S10M2S4H\t1\tretained"
    ));
    assert!(rows[1].contains("\tfull_minus\tchr1\t149\t159\t-\t43\t2\t5\t3H5S10M2S4H\t1\tretained"));
    assert!(rows[2].contains("\tlow_mapq\t") && rows[2].contains("\tbelow_min_mapq\t"));
    assert!(rows[3].contains("\tsecondary\t") && rows[3].contains("\tsecondary\t"));
    assert!(rows[4].contains("\tsupplementary\t") && rows[4].contains("\tsupplementary\t"));
    assert!(rows[5].contains("\tspliced\t") && rows[5].contains("\tspliced_or_skipped\t"));
    assert!(rows[6].contains("\tunmapped\t") && rows[6].contains("\tunmapped\t"));
    assert!(rows[7].contains("\t.\t") && rows[7].contains("\tmissing_name\t"));

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("falling back to in-memory record buffering and BED sorting"));
    assert!(
        stderr.trim().ends_with(concat!(
            "bam_to_bed_counts\ttotal=8\tretained=2\tunmapped=1\tsecondary=1",
            "\tsupplementary=1\tspliced_or_skipped=1\tsoft_clipped=2\tmissing_name=1",
            "\tmalformed_record=0\tbelow_min_mapq=1\tmissing_full_length_evidence=0",
            "\tbelow_min_boundary_support=0\tinvalid_alignment=0",
        )),
        "stderr:\n{stderr}"
    );
    let _ = fs::remove_dir_all(tmp);
}

#[test]
fn trackclustertu_bam_to_bed_filters_low_exact_boundary_support() {
    let tmp = unique_tmp_dir("trackclustertu_bam_support_test");
    fs::create_dir_all(&tmp).unwrap();
    let bam_path = tmp.join("reads.bam");
    let bed_path = tmp.join("reads.bed");
    let evidence_path = tmp.join("reads.evidence.tsv");
    let header = build_header();
    let mut writer = bam::io::Writer::new(fs::File::create(&bam_path).unwrap());
    writer.write_header(&header).unwrap();
    for name in ["supported_a", "supported_b"] {
        writer
            .write_alignment_record(
                &header,
                &mapped_record(name, 0, 99, match_cigar(10), 30, false, false),
            )
            .unwrap();
    }
    writer
        .write_alignment_record(
            &header,
            &mapped_record("singleton", 0, 199, match_cigar(10), 30, false, false),
        )
        .unwrap();
    writer.try_finish().unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_trackclustertu"))
        .args([
            "bam-to-bed",
            "--in-bam",
            bam_path.to_str().unwrap(),
            "--out-bed",
            bed_path.to_str().unwrap(),
            "--out-evidence",
            evidence_path.to_str().unwrap(),
            "--min-boundary-support",
            "2",
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        fs::read_to_string(&bed_path).unwrap(),
        concat!(
            "chr1\t99\t109\tsupported_a\t30\t+\n",
            "chr1\t99\t109\tsupported_b\t30\t+\n",
        )
    );
    let evidence = fs::read_to_string(&evidence_path).unwrap();
    assert!(evidence
        .lines()
        .any(|line| line.contains("\tsingleton\t")
            && line.contains("\tbelow_min_boundary_support\t")));
    assert!(String::from_utf8_lossy(&output.stderr).contains("below_min_boundary_support=1"));
    let _ = fs::remove_dir_all(tmp);
}

#[test]
fn trackclustertu_bam_manifest_accepts_optional_evidence_and_profile_columns() {
    let tmp = unique_tmp_dir("trackclustertu_bam_manifest_evidence_test");
    fs::create_dir_all(&tmp).unwrap();
    let bam_path = tmp.join("reads.bam");
    let evidence_input = tmp.join("sample.input-evidence.tsv");
    let manifest_path = tmp.join("samples.tsv");
    let out_dir = tmp.join("out");
    let header = build_header();
    let mut writer = bam::io::Writer::new(fs::File::create(&bam_path).unwrap());
    writer.write_header(&header).unwrap();
    writer
        .write_alignment_record(
            &header,
            &mapped_record("r1", 0, 99, match_cigar(10), 30, false, false),
        )
        .unwrap();
    writer.try_finish().unwrap();
    fs::write(
        &evidence_input,
        concat!(
            "read_name\tfive_prime_adapter_evidence\tthree_prime_adapter_evidence\n",
            "r1\tpresent\tpresent\n",
        ),
    )
    .unwrap();
    fs::write(
        &manifest_path,
        format!(
            "sample\treads\tevidence\tlibrary_profile\nS1\t{}\t{}\tdirect-cdna\n",
            bam_path.display(),
            evidence_input.display(),
        ),
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_trackclustertu"))
        .args([
            "bam-to-bed",
            "--manifest",
            manifest_path.to_str().unwrap(),
            "--out-dir",
            out_dir.to_str().unwrap(),
            "--emit-evidence",
            "--require-full-length",
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
        fs::read_to_string(out_dir.join("bed/S1.bed")).unwrap(),
        "chr1\t99\t109\tr1\t30\t+\n"
    );
    let manifest = fs::read_to_string(out_dir.join("samples.bed.tsv")).unwrap();
    assert!(manifest.starts_with("sample\tgroup\treads\tevidence\tlibrary_profile\n"));
    assert!(manifest.contains("\tdirect-cdna\n"));
    let evidence = fs::read_to_string(out_dir.join("bed/S1.evidence.tsv")).unwrap();
    assert!(evidence.contains("\tdirect-cdna\tfive_prime_adapter_and_three_prime_adapter\t"));

    let cluster_out = tmp.join("cluster-out");
    let cluster = Command::new(env!("CARGO_BIN_EXE_trackclustertu"))
        .args([
            "cluster",
            "--manifest",
            out_dir.join("samples.bed.tsv").to_str().unwrap(),
            "--format",
            "bed6",
            "--out-dir",
            cluster_out.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        cluster.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&cluster.stdout),
        String::from_utf8_lossy(&cluster.stderr)
    );

    let membership = fs::read_to_string(cluster_out.join("membership.tsv")).unwrap();
    let membership_row = membership
        .lines()
        .find(|line| line.starts_with("S1::r1\t"))
        .expect("cluster membership row for evidence-backed read");
    assert_eq!(membership_row.split('\t').next_back(), Some("true"));

    let count = fs::read_to_string(cluster_out.join("tu_count.csv")).unwrap();
    let count_fields: Vec<&str> = count.lines().nth(1).unwrap().split(',').collect();
    assert_eq!(count_fields[1..], ["1", "1", "1", "1", "1"]);

    let sample_long = fs::read_to_string(cluster_out.join("tu_sample_long.tsv")).unwrap();
    let sample_fields: Vec<&str> = sample_long.lines().nth(1).unwrap().split('\t').collect();
    assert_eq!(sample_fields[1..], ["S1", "1", "1", "1", "1"]);
    let _ = fs::remove_dir_all(tmp);
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
