use std::collections::HashMap;

use crate::model::Transcript;

use super::StrandMode;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct PartitionKey {
    pub chrom: String,
    pub strand: Option<crate::model::Strand>,
}

pub fn partition(
    transcripts: &[Transcript],
    strand_mode: StrandMode,
) -> HashMap<PartitionKey, Vec<usize>> {
    let mut parts: HashMap<PartitionKey, Vec<usize>> = HashMap::new();
    for (index, transcript) in transcripts.iter().enumerate() {
        let key = PartitionKey {
            chrom: transcript.chrom.clone(),
            strand: strand_mode.key_strand(transcript.strand),
        };
        parts.entry(key).or_default().push(index);
    }
    parts
}

#[cfg(test)]
mod tests {
    use crate::interval::sort::sort_by_coord;
    use crate::model::{Bed12Attrs, Coord, Interval, Strand, Transcript};

    use super::*;

    fn make_tx(chrom: &str, strand: Strand, start: u32, end: u32, name: &str) -> Transcript {
        Transcript::new(
            chrom.to_owned(),
            strand,
            Coord::new(start),
            Coord::new(end),
            name.to_owned(),
            vec![Interval::new(Coord::new(start), Coord::new(end)).unwrap()],
            Bed12Attrs {
                score: 0,
                thick_start: Coord::new(start),
                thick_end: Coord::new(end),
                item_rgb: "0".to_owned(),
                extra_fields: Vec::new(),
            },
        )
        .unwrap()
    }

    #[test]
    fn partition_splits_by_chrom_only_when_ignore_strand() {
        let mut transcripts = vec![
            make_tx("chr1", Strand::Plus, 10, 20, "a"),
            make_tx("chr1", Strand::Minus, 11, 21, "b"),
            make_tx("chr2", Strand::Plus, 0, 1, "c"),
        ];
        sort_by_coord(&mut transcripts);

        let parts = partition(&transcripts, StrandMode::Ignore);
        assert_eq!(parts.len(), 2);

        let chr1 = parts
            .get(&PartitionKey {
                chrom: "chr1".to_owned(),
                strand: None,
            })
            .unwrap();
        assert_eq!(chr1.len(), 2);
    }

    #[test]
    fn partition_splits_by_chrom_and_strand_when_match_strand() {
        let mut transcripts = vec![
            make_tx("chr1", Strand::Plus, 10, 20, "a"),
            make_tx("chr1", Strand::Minus, 11, 21, "b"),
        ];
        sort_by_coord(&mut transcripts);

        let parts = partition(&transcripts, StrandMode::Match);
        assert_eq!(parts.len(), 2);
    }
}
