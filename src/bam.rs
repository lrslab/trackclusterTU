use std::fs::File;
use std::path::Path;

use anyhow::{Context, Result};
use noodles::{
    bam as noodles_bam,
    sam::alignment::{record::cigar::op::Kind as CigarKind, RecordBuf},
};

use crate::io::bed::{self, Bed6Record};
use crate::model::{Coord, Strand};

pub fn bam_to_bed6(bam_path: &Path, out_bed_path: &Path) -> Result<()> {
    let mut reader = File::open(bam_path)
        .with_context(|| format!("failed to open BAM {:?}", bam_path))
        .map(noodles_bam::io::Reader::new)?;
    let header = reader.read_header()?;

    let mut bed_records = Vec::new();
    let mut record = RecordBuf::default();
    while reader.read_record_buf(&header, &mut record)? != 0 {
        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }
        if record
            .cigar()
            .as_ref()
            .iter()
            .any(|op| op.kind() == CigarKind::Skip)
        {
            continue;
        }

        let Some(reference_sequence_id) = record.reference_sequence_id() else {
            continue;
        };
        let Some((reference_name, _)) = header
            .reference_sequences()
            .get_index(reference_sequence_id)
        else {
            continue;
        };
        let Some(alignment_start) = record.alignment_start() else {
            continue;
        };
        let Some(alignment_end) = record.alignment_end() else {
            continue;
        };

        let start = Coord::new((usize::from(alignment_start) - 1) as u32);
        let end = Coord::new(usize::from(alignment_end) as u32);
        let strand = if flags.is_reverse_complemented() {
            Strand::Minus
        } else {
            Strand::Plus
        };
        let score = record.mapping_quality().map(u8::from).unwrap_or(0) as u32;

        bed_records.push(Bed6Record {
            chrom: reference_name.to_string(),
            start,
            end,
            name: record
                .name()
                .map(|name| name.to_string())
                .unwrap_or_else(|| ".".to_owned()),
            score,
            strand,
            extra_fields: Vec::new(),
        });
    }

    bed_records.sort_by(|left, right| {
        left.chrom
            .cmp(&right.chrom)
            .then_with(|| left.start.cmp(&right.start))
            .then_with(|| left.end.cmp(&right.end))
            .then_with(|| left.strand.cmp(&right.strand))
            .then_with(|| left.name.cmp(&right.name))
    });
    bed::write_bed6(out_bed_path, bed_records.iter())?;
    Ok(())
}
