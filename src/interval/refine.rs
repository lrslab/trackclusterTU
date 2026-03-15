use crate::model::Transcript;

pub fn exonic_overlap_bp(a: &Transcript, b: &Transcript) -> u32 {
    let mut total: u32 = 0;
    let mut ai: usize = 0;
    let mut bi: usize = 0;

    while ai < a.exons.len() && bi < b.exons.len() {
        let a_exon = a.exons[ai];
        let b_exon = b.exons[bi];

        if a_exon.end <= b_exon.start {
            ai += 1;
            continue;
        }
        if b_exon.end <= a_exon.start {
            bi += 1;
            continue;
        }

        total += a_exon.overlap_len(b_exon);

        if a_exon.end <= b_exon.end {
            ai += 1;
        } else {
            bi += 1;
        }
    }

    total
}

pub fn junctions_equal(a: &Transcript, b: &Transcript) -> bool {
    a.junction_signature() == b.junction_signature()
}

pub fn junctions_subset(a: &Transcript, b: &Transcript) -> bool {
    let a_introns = a.introns();
    let b_introns = b.introns();

    let mut ai: usize = 0;
    let mut bi: usize = 0;
    while ai < a_introns.len() && bi < b_introns.len() {
        match a_introns[ai].start.cmp(&b_introns[bi].start) {
            std::cmp::Ordering::Equal => match a_introns[ai].end.cmp(&b_introns[bi].end) {
                std::cmp::Ordering::Equal => {
                    ai += 1;
                    bi += 1;
                }
                std::cmp::Ordering::Greater => {
                    bi += 1;
                }
                std::cmp::Ordering::Less => {
                    return false;
                }
            },
            std::cmp::Ordering::Greater => {
                bi += 1;
            }
            std::cmp::Ordering::Less => {
                return false;
            }
        }
    }
    ai == a_introns.len()
}

#[cfg(test)]
mod tests {
    use crate::model::{Bed12Attrs, Coord, Strand, Transcript};

    use super::*;

    fn make_tx(chrom: &str, strand: Strand, name: &str, exons: &[(u32, u32)]) -> Transcript {
        let tx_start = exons.iter().map(|(s, _)| *s).min().unwrap_or(0);
        let tx_end = exons.iter().map(|(_, e)| *e).max().unwrap_or(0);
        let exons = exons
            .iter()
            .map(|(s, e)| crate::model::Interval::new(Coord::new(*s), Coord::new(*e)).unwrap())
            .collect::<Vec<_>>();

        Transcript::new(
            chrom.to_owned(),
            strand,
            Coord::new(tx_start),
            Coord::new(tx_end),
            name.to_owned(),
            exons,
            Bed12Attrs {
                score: 0,
                thick_start: Coord::new(tx_start),
                thick_end: Coord::new(tx_end),
                item_rgb: "0".to_owned(),
                extra_fields: Vec::new(),
            },
        )
        .unwrap()
    }

    #[test]
    fn exonic_overlap_bp_is_zero_when_disjoint() {
        let a = make_tx("chr1", Strand::Plus, "a", &[(0, 10)]);
        let b = make_tx("chr1", Strand::Plus, "b", &[(20, 30)]);
        assert_eq!(exonic_overlap_bp(&a, &b), 0);
    }

    #[test]
    fn exonic_overlap_bp_counts_shared_bases() {
        let a = make_tx("chr1", Strand::Plus, "a", &[(0, 10), (20, 30)]);
        let b = make_tx("chr1", Strand::Plus, "b", &[(5, 25)]);
        assert_eq!(exonic_overlap_bp(&a, &b), 10);
    }

    #[test]
    fn junction_subset_and_equal() {
        let a = make_tx("chr1", Strand::Plus, "a", &[(0, 10), (20, 30), (40, 50)]);
        let b = make_tx("chr1", Strand::Plus, "b", &[(0, 10), (20, 30), (40, 50)]);
        assert!(junctions_equal(&a, &b));
        assert!(junctions_subset(&a, &b));

        let c = make_tx("chr1", Strand::Plus, "c", &[(0, 10), (20, 30)]);
        assert!(!junctions_equal(&c, &a));
        assert!(junctions_subset(&c, &a));
        assert!(!junctions_subset(&a, &c));
    }
}
