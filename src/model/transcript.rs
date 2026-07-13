//! Checked BED12-compatible transcripts and junction signatures.

use std::cmp::Ordering;

use thiserror::Error;

use super::{Coord, Interval, Strand};

#[derive(Clone, Debug, PartialEq, Eq)]
/// A transcript with a valid span and sorted, positive, non-overlapping exons.
pub struct Transcript {
    /// Reference sequence name.
    chrom: String,
    /// Transcript strand.
    strand: Strand,
    /// Inclusive transcript start.
    tx_start: Coord,
    /// Exclusive transcript end.
    tx_end: Coord,
    /// Transcript identifier.
    name: String,
    /// BED score.
    score: u32,
    /// BED thick-region start.
    thick_start: Coord,
    /// BED thick-region end.
    thick_end: Coord,
    /// BED RGB field.
    item_rgb: String,
    /// Sorted, positive, non-overlapping exon intervals within the transcript span.
    exons: Vec<Interval>,
    /// Columns retained after the standard BED12 fields.
    extra_fields: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
/// Non-structural BED12 attributes supplied when constructing a [`Transcript`].
pub struct Bed12Attrs {
    /// BED score.
    pub score: u32,
    /// BED thick-region start.
    pub thick_start: Coord,
    /// BED thick-region end.
    pub thick_end: Coord,
    /// BED RGB field.
    pub item_rgb: String,
    /// Columns retained after the standard BED12 fields.
    pub extra_fields: Vec<String>,
}

#[derive(Error, Debug)]
/// Validation error returned by [`Transcript::new`].
pub enum TranscriptError {
    /// The transcript start is greater than its end.
    #[error("invalid transcript span: start {start} > end {end}")]
    InvalidSpan {
        /// Proposed transcript start.
        start: Coord,
        /// Proposed transcript end.
        end: Coord,
    },

    /// No exons were supplied.
    #[error("expected at least 1 exon")]
    EmptyExons,

    /// An exon has zero length.
    #[error("exon must have positive length: {exon:?}")]
    EmptyExon {
        /// Invalid empty exon.
        exon: Interval,
    },

    /// An exon extends beyond the declared transcript span.
    #[error("exon is outside transcript span: exon {exon:?}, transcript [{tx_start}, {tx_end})")]
    ExonOutsideSpan {
        /// Exon outside the span.
        exon: Interval,
        /// Declared transcript start.
        tx_start: Coord,
        /// Declared transcript end.
        tx_end: Coord,
    },

    /// Two sorted exons overlap.
    #[error("transcript exons overlap: {left:?} and {right:?}")]
    OverlappingExons {
        /// Earlier overlapping exon.
        left: Interval,
        /// Later overlapping exon.
        right: Interval,
    },
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
/// Reference, strand, and intron coordinates used for junction comparison.
pub struct JunctionSignature {
    /// Reference sequence name.
    pub chrom: String,
    /// Transcript strand.
    pub strand: Strand,
    /// Ordered intron intervals.
    pub introns: Vec<Interval>,
}

impl Transcript {
    /// Construct a transcript after sorting and validating its exon intervals.
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

        exons.sort_by(|left, right| match left.start().cmp(&right.start()) {
            Ordering::Equal => left.end().cmp(&right.end()),
            ordering => ordering,
        });

        for exon in &exons {
            if exon.is_empty() {
                return Err(TranscriptError::EmptyExon { exon: *exon });
            }
            if exon.start() < tx_start || exon.end() > tx_end {
                return Err(TranscriptError::ExonOutsideSpan {
                    exon: *exon,
                    tx_start,
                    tx_end,
                });
            }
        }

        for window in exons.windows(2) {
            let left = window[0];
            let right = window[1];
            if left.end() > right.start() {
                return Err(TranscriptError::OverlappingExons { left, right });
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

    /// Reference sequence name.
    pub fn chrom(&self) -> &str {
        &self.chrom
    }

    /// Transcript strand.
    pub const fn strand(&self) -> Strand {
        self.strand
    }

    /// Inclusive transcript start.
    pub const fn tx_start(&self) -> Coord {
        self.tx_start
    }

    /// Exclusive transcript end.
    pub const fn tx_end(&self) -> Coord {
        self.tx_end
    }

    /// Transcript identifier/name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// BED score.
    pub const fn score(&self) -> u32 {
        self.score
    }

    /// BED thick-region start.
    pub const fn thick_start(&self) -> Coord {
        self.thick_start
    }

    /// BED thick-region end.
    pub const fn thick_end(&self) -> Coord {
        self.thick_end
    }

    /// BED RGB field.
    pub fn item_rgb(&self) -> &str {
        &self.item_rgb
    }

    /// Validated, sorted exon intervals.
    pub fn exons(&self) -> &[Interval] {
        &self.exons
    }

    /// Extra BED columns retained after the standard BED12 fields.
    pub fn extra_fields(&self) -> &[String] {
        &self.extra_fields
    }

    /// Derive the non-empty introns between adjacent exons.
    pub fn introns(&self) -> Vec<Interval> {
        let mut introns = Vec::new();
        for window in self.exons.windows(2) {
            let left = window[0];
            let right = window[1];
            if left.end() < right.start() {
                introns.push(
                    Interval::new(left.end(), right.start())
                        .expect("validated exon order yields a valid intron"),
                );
            }
        }
        introns
    }

    /// Return the reference, strand, and intron signature used for junction grouping.
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
    fn getters_expose_constructor_normalized_state() {
        let transcript = Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
            Coord::new(10),
            Coord::new(50),
            "tx".to_owned(),
            vec![
                Interval::new(Coord::new(30), Coord::new(50)).unwrap(),
                Interval::new(Coord::new(10), Coord::new(20)).unwrap(),
            ],
            bed_attrs(10, 50),
        )
        .unwrap();

        assert_eq!(transcript.chrom(), "chr1");
        assert_eq!(transcript.strand(), Strand::Plus);
        assert_eq!(transcript.tx_start(), Coord::new(10));
        assert_eq!(transcript.tx_end(), Coord::new(50));
        assert_eq!(transcript.name(), "tx");
        assert_eq!(
            transcript.exons(),
            &[
                Interval::new(Coord::new(10), Coord::new(20)).unwrap(),
                Interval::new(Coord::new(30), Coord::new(50)).unwrap(),
            ]
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

    fn bed_attrs(start: u32, end: u32) -> Bed12Attrs {
        Bed12Attrs {
            score: 0,
            thick_start: Coord::new(start),
            thick_end: Coord::new(end),
            item_rgb: "0".to_owned(),
            extra_fields: Vec::new(),
        }
    }

    #[test]
    fn transcript_rejects_empty_exons() {
        let error = Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
            Coord::new(10),
            Coord::new(20),
            "tx".to_owned(),
            vec![Interval::new(Coord::new(15), Coord::new(15)).unwrap()],
            bed_attrs(10, 20),
        )
        .unwrap_err();

        assert!(matches!(error, TranscriptError::EmptyExon { .. }));
    }

    #[test]
    fn transcript_rejects_overlapping_exons_after_sorting() {
        let error = Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
            Coord::new(10),
            Coord::new(40),
            "tx".to_owned(),
            vec![
                Interval::new(Coord::new(20), Coord::new(35)).unwrap(),
                Interval::new(Coord::new(10), Coord::new(25)).unwrap(),
            ],
            bed_attrs(10, 40),
        )
        .unwrap_err();

        assert!(matches!(error, TranscriptError::OverlappingExons { .. }));
    }

    #[test]
    fn transcript_rejects_exons_outside_declared_span() {
        let error = Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
            Coord::new(10),
            Coord::new(20),
            "tx".to_owned(),
            vec![Interval::new(Coord::new(10), Coord::new(21)).unwrap()],
            bed_attrs(10, 20),
        )
        .unwrap_err();

        assert!(matches!(error, TranscriptError::ExonOutsideSpan { .. }));
    }
}
