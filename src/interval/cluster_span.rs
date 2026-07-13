//! Sweep-like grouping of overlapping transcript spans.

use std::collections::BTreeMap;

use crate::model::{Interval, Transcript};

use super::StrandMode;

/// One connected group of overlapping transcript spans.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RangeCluster {
    /// Reference sequence name.
    pub chrom: String,
    /// Strand key when strand matching is enabled.
    pub strand: Option<crate::model::Strand>,
    /// Half-open span covering every member.
    pub span: Interval,
    /// Indices into the input transcript slice.
    pub members: Vec<usize>,
}

fn span(transcript: &Transcript) -> Interval {
    Interval::new(transcript.tx_start(), transcript.tx_end())
        .expect("validated transcript span has ordered bounds")
}

/// Group coordinate-sorted transcripts into connected span-overlap clusters.
///
/// `records_sorted` must first be ordered with [`super::sort_by_coord`]. Strict
/// half-open overlap joins records transitively; merely touching spans begin a
/// new cluster. [`StrandMode::Ignore`] partitions only by reference sequence,
/// while [`StrandMode::Match`] also partitions by strand.
///
/// Each reference/strand partition is maintained independently, so interleaved
/// strands cannot prematurely close one another's active cluster. Returned
/// clusters follow the encounter order of their first member, and each
/// [`RangeCluster::members`] list contains indices into `records_sorted`.
pub fn cluster_by_span(
    records_sorted: &[Transcript],
    strand_mode: StrandMode,
) -> Vec<RangeCluster> {
    let mut partitions: BTreeMap<(String, Option<crate::model::Strand>), Vec<RangeCluster>> =
        BTreeMap::new();

    for (index, transcript) in records_sorted.iter().enumerate() {
        let tx_span = span(transcript);
        let key_strand = strand_mode.key_strand(transcript.strand());
        let clusters = partitions
            .entry((transcript.chrom().to_owned(), key_strand))
            .or_default();

        match clusters.last_mut() {
            None => {
                clusters.push(RangeCluster {
                    chrom: transcript.chrom().to_owned(),
                    strand: key_strand,
                    span: tx_span,
                    members: vec![index],
                });
            }
            Some(cluster) => {
                if tx_span.start() < cluster.span.end() {
                    if tx_span.end() > cluster.span.end() {
                        cluster.span = Interval::new(cluster.span.start(), tx_span.end())
                            .expect("cluster extension retains ordered bounds");
                    }
                    cluster.members.push(index);
                } else {
                    clusters.push(RangeCluster {
                        chrom: transcript.chrom().to_owned(),
                        strand: key_strand,
                        span: tx_span,
                        members: vec![index],
                    });
                }
            }
        }
    }

    let mut clusters: Vec<RangeCluster> = partitions.into_values().flatten().collect();
    // Restore the encounter order of cluster starts from `records_sorted` after
    // processing each (contig, strand) partition independently.
    clusters.sort_by_key(|cluster| cluster.members[0]);

    clusters
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
    fn cluster_by_span_groups_overlapping_records() {
        let mut records = vec![
            make_tx("chr1", Strand::Plus, 10, 20, "a"),
            make_tx("chr1", Strand::Plus, 15, 25, "b"),
            make_tx("chr1", Strand::Plus, 30, 40, "c"),
        ];
        sort_by_coord(&mut records);

        let clusters = cluster_by_span(&records, StrandMode::Ignore);
        assert_eq!(clusters.len(), 2);
        assert_eq!(clusters[0].members.len(), 2);
        assert_eq!(clusters[1].members.len(), 1);
    }

    #[test]
    fn cluster_by_span_respects_strand_mode() {
        let mut records = vec![
            make_tx("chr1", Strand::Plus, 10, 20, "a"),
            make_tx("chr1", Strand::Minus, 15, 25, "b"),
        ];
        sort_by_coord(&mut records);

        let clusters = cluster_by_span(&records, StrandMode::Ignore);
        assert_eq!(clusters.len(), 1);

        let clusters = cluster_by_span(&records, StrandMode::Match);
        assert_eq!(clusters.len(), 2);
    }

    #[test]
    fn cluster_by_span_keeps_interleaved_strand_partitions_open() {
        let mut records = vec![
            make_tx("chr1", Strand::Plus, 10, 100, "plus_outer"),
            make_tx("chr1", Strand::Minus, 20, 30, "minus"),
            make_tx("chr1", Strand::Plus, 40, 50, "plus_inner"),
        ];
        sort_by_coord(&mut records);

        let clusters = cluster_by_span(&records, StrandMode::Match);
        assert_eq!(clusters.len(), 2);
        assert_eq!(clusters[0].strand, Some(Strand::Plus));
        assert_eq!(clusters[0].members, vec![0, 2]);
        assert_eq!(clusters[1].strand, Some(Strand::Minus));
        assert_eq!(clusters[1].members, vec![1]);
    }
}
