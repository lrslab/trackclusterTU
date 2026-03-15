use crate::model::{Interval, Transcript};

use super::StrandMode;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RangeCluster {
    pub chrom: String,
    pub strand: Option<crate::model::Strand>,
    pub span: Interval,
    pub members: Vec<usize>,
}

fn span(transcript: &Transcript) -> Interval {
    Interval {
        start: transcript.tx_start,
        end: transcript.tx_end,
    }
}

pub fn cluster_by_span(
    records_sorted: &[Transcript],
    strand_mode: StrandMode,
) -> Vec<RangeCluster> {
    let mut clusters: Vec<RangeCluster> = Vec::new();

    let mut current: Option<RangeCluster> = None;
    for (index, transcript) in records_sorted.iter().enumerate() {
        let tx_span = span(transcript);
        let key_strand = strand_mode.key_strand(transcript.strand);

        match current.as_mut() {
            None => {
                current = Some(RangeCluster {
                    chrom: transcript.chrom.clone(),
                    strand: key_strand,
                    span: tx_span,
                    members: vec![index],
                });
            }
            Some(cluster) => {
                if transcript.chrom != cluster.chrom || key_strand != cluster.strand {
                    clusters.push(current.take().unwrap());
                    current = Some(RangeCluster {
                        chrom: transcript.chrom.clone(),
                        strand: key_strand,
                        span: tx_span,
                        members: vec![index],
                    });
                    continue;
                }

                if tx_span.start < cluster.span.end {
                    if tx_span.end > cluster.span.end {
                        cluster.span.end = tx_span.end;
                    }
                    cluster.members.push(index);
                } else {
                    clusters.push(current.take().unwrap());
                    current = Some(RangeCluster {
                        chrom: transcript.chrom.clone(),
                        strand: key_strand,
                        span: tx_span,
                        members: vec![index],
                    });
                }
            }
        }
    }

    if let Some(cluster) = current {
        clusters.push(cluster);
    }

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
}
