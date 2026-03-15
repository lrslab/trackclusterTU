use std::cmp::Ordering;

use thiserror::Error;

use super::{Coord, Interval, Strand};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Transcript {
    pub chrom: String,
    pub strand: Strand,
    pub tx_start: Coord,
    pub tx_end: Coord,
    pub name: String,
    pub score: u32,
    pub thick_start: Coord,
    pub thick_end: Coord,
    pub item_rgb: String,
    pub exons: Vec<Interval>,
    pub extra_fields: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Bed12Attrs {
    pub score: u32,
    pub thick_start: Coord,
    pub thick_end: Coord,
    pub item_rgb: String,
    pub extra_fields: Vec<String>,
}

#[derive(Error, Debug)]
pub enum TranscriptError {
    #[error("invalid transcript span: start {start} > end {end}")]
    InvalidSpan { start: Coord, end: Coord },

    #[error("expected at least 1 exon")]
    EmptyExons,

    #[error("exon is outside transcript span: exon {exon:?}, transcript [{tx_start}, {tx_end})")]
    ExonOutsideSpan {
        exon: Interval,
        tx_start: Coord,
        tx_end: Coord,
    },
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct JunctionSignature {
    pub chrom: String,
    pub strand: Strand,
    pub introns: Vec<Interval>,
}

impl Transcript {
    pub fn new(
        chrom: String,
        strand: Strand,
        tx_start: Coord,
        tx_end: Coord,
        name: String,
        mut exons: Vec<Interval>,
        bed: Bed12Attrs,
    ) -> Result<Self, TranscriptError> {
        if tx_start > tx_end {
            return Err(TranscriptError::InvalidSpan {
                start: tx_start,
                end: tx_end,
            });
        }
        if exons.is_empty() {
            return Err(TranscriptError::EmptyExons);
        }

        exons.sort_by(|left, right| match left.start.cmp(&right.start) {
            Ordering::Equal => left.end.cmp(&right.end),
            ordering => ordering,
        });

        for exon in &exons {
            if exon.start < tx_start || exon.end > tx_end {
                return Err(TranscriptError::ExonOutsideSpan {
                    exon: *exon,
                    tx_start,
                    tx_end,
                });
            }
        }

        Ok(Self {
            chrom,
            strand,
            tx_start,
            tx_end,
            name,
            score: bed.score,
            thick_start: bed.thick_start,
            thick_end: bed.thick_end,
            item_rgb: bed.item_rgb,
            exons,
            extra_fields: bed.extra_fields,
        })
    }

    pub fn introns(&self) -> Vec<Interval> {
        let mut introns = Vec::new();
        for window in self.exons.windows(2) {
            let left = window[0];
            let right = window[1];
            if left.end < right.start {
                introns.push(Interval {
                    start: left.end,
                    end: right.start,
                });
            }
        }
        introns
    }

    pub fn junction_signature(&self) -> JunctionSignature {
        JunctionSignature {
            chrom: self.chrom.clone(),
            strand: self.strand,
            introns: self.introns(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn introns_from_exons() {
        let exons = vec![
            Interval::new(Coord::new(100), Coord::new(150)).unwrap(),
            Interval::new(Coord::new(170), Coord::new(200)).unwrap(),
        ];
        let transcript = Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
            Coord::new(100),
            Coord::new(200),
            "tx1".to_owned(),
            exons,
            Bed12Attrs {
                score: 0,
                thick_start: Coord::new(100),
                thick_end: Coord::new(200),
                item_rgb: "0".to_owned(),
                extra_fields: Vec::new(),
            },
        )
        .unwrap();

        assert_eq!(
            transcript.introns(),
            vec![Interval::new(Coord::new(150), Coord::new(170)).unwrap()]
        );
    }

    #[test]
    fn junction_signature_is_stable() {
        let exons_a = vec![
            Interval::new(Coord::new(10), Coord::new(20)).unwrap(),
            Interval::new(Coord::new(30), Coord::new(40)).unwrap(),
        ];
        let exons_b = vec![
            Interval::new(Coord::new(10), Coord::new(20)).unwrap(),
            Interval::new(Coord::new(30), Coord::new(40)).unwrap(),
        ];

        let a = Transcript::new(
            "chr1".to_owned(),
            Strand::Minus,
            Coord::new(10),
            Coord::new(40),
            "a".to_owned(),
            exons_a,
            Bed12Attrs {
                score: 0,
                thick_start: Coord::new(10),
                thick_end: Coord::new(40),
                item_rgb: "0".to_owned(),
                extra_fields: Vec::new(),
            },
        )
        .unwrap();

        let b = Transcript::new(
            "chr1".to_owned(),
            Strand::Minus,
            Coord::new(10),
            Coord::new(40),
            "b".to_owned(),
            exons_b,
            Bed12Attrs {
                score: 999,
                thick_start: Coord::new(10),
                thick_end: Coord::new(40),
                item_rgb: "0".to_owned(),
                extra_fields: vec!["extra".to_owned()],
            },
        )
        .unwrap();

        assert_eq!(a.junction_signature(), b.junction_signature());
    }
}
