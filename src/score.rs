use crate::interval::{intersection_len, total_len, union_len};
use crate::model::{Interval, Transcript};

pub fn score1_interval(a: Interval, b: Interval) -> f64 {
    let overlap = a.overlap_len(b) as u64;
    let union = a.len() as u64 + b.len() as u64 - overlap;
    if union == 0 {
        0.0
    } else {
        overlap as f64 / union as f64
    }
}

pub fn score2_interval(a: Interval, b: Interval) -> f64 {
    let overlap = a.overlap_len(b) as u64;
    let max_len = (a.len() as u64).max(b.len() as u64);
    if max_len == 0 {
        0.0
    } else {
        overlap as f64 / max_len as f64
    }
}

pub fn score1_intervals(a: &[Interval], b: &[Interval]) -> f64 {
    let overlap = intersection_len(a, b);
    let union = union_len(a, b);
    if union == 0 {
        0.0
    } else {
        overlap as f64 / union as f64
    }
}

pub fn score2_intervals(a: &[Interval], b: &[Interval]) -> f64 {
    let overlap = intersection_len(a, b);
    let max_len = total_len(a).max(total_len(b));
    if max_len == 0 {
        0.0
    } else {
        overlap as f64 / max_len as f64
    }
}

pub fn score1_transcripts(a: &Transcript, b: &Transcript, intron_weight: f64) -> f64 {
    debug_assert!(intron_weight >= 0.0);

    let exon_overlap = intersection_len(&a.exons, &b.exons) as f64;
    let exon_union = union_len(&a.exons, &b.exons) as f64;

    if intron_weight == 0.0 {
        if exon_union == 0.0 {
            0.0
        } else {
            exon_overlap / exon_union
        }
    } else {
        let a_introns = a.introns();
        let b_introns = b.introns();
        let intron_overlap = intersection_len(&a_introns, &b_introns) as f64;
        let intron_union = union_len(&a_introns, &b_introns) as f64;

        let overlap = exon_overlap + intron_weight * intron_overlap;
        let union = exon_union + intron_weight * intron_union;
        if union == 0.0 {
            0.0
        } else {
            overlap / union
        }
    }
}

pub fn score2_transcripts(a: &Transcript, b: &Transcript, intron_weight: f64) -> f64 {
    debug_assert!(intron_weight >= 0.0);

    let exon_overlap = intersection_len(&a.exons, &b.exons) as f64;
    let exon_len_a = total_len(&a.exons) as f64;
    let exon_len_b = total_len(&b.exons) as f64;

    if intron_weight == 0.0 {
        let max_len = exon_len_a.max(exon_len_b);
        if max_len == 0.0 {
            0.0
        } else {
            exon_overlap / max_len
        }
    } else {
        let a_introns = a.introns();
        let b_introns = b.introns();

        let intron_overlap = intersection_len(&a_introns, &b_introns) as f64;
        let intron_len_a = total_len(&a_introns) as f64;
        let intron_len_b = total_len(&b_introns) as f64;

        let overlap = exon_overlap + intron_weight * intron_overlap;
        let len_a = exon_len_a + intron_weight * intron_len_a;
        let len_b = exon_len_b + intron_weight * intron_len_b;
        let max_len = len_a.max(len_b);
        if max_len == 0.0 {
            0.0
        } else {
            overlap / max_len
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::model::{Bed12Attrs, Coord, Strand};

    use super::*;

    fn interval(start: u32, end: u32) -> Interval {
        Interval::new(Coord::new(start), Coord::new(end)).unwrap()
    }

    #[test]
    fn score1_interval_is_one_for_identical() {
        let a = interval(0, 10);
        assert_eq!(score1_interval(a, a), 1.0);
    }

    #[test]
    fn score1_interval_is_zero_for_disjoint() {
        let a = interval(0, 10);
        let b = interval(10, 20);
        assert_eq!(score1_interval(a, b), 0.0);
        assert_eq!(score2_interval(a, b), 0.0);
    }

    #[test]
    fn score2_interval_penalizes_containment_by_length_ratio() {
        let a = interval(0, 10);
        let b = interval(2, 8);
        assert_eq!(score2_interval(a, b), 0.6);
        assert_eq!(score2_interval(b, a), 0.6);
        assert_eq!(score1_interval(a, b), 0.6);
    }

    #[test]
    fn scores_return_zero_when_denominator_is_zero() {
        let empty = interval(5, 5);
        assert_eq!(score1_interval(empty, empty), 0.0);
        assert_eq!(score2_interval(empty, empty), 0.0);
    }

    fn make_tx(name: &str, exons: &[(u32, u32)]) -> Transcript {
        let tx_start = exons.iter().map(|(s, _)| *s).min().unwrap_or(0);
        let tx_end = exons.iter().map(|(_, e)| *e).max().unwrap_or(0);
        let exons = exons
            .iter()
            .map(|(s, e)| interval(*s, *e))
            .collect::<Vec<_>>();

        Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
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
    fn transcript_scores_are_symmetric_without_introns() {
        let a = make_tx("a", &[(0, 10)]);
        let b = make_tx("b", &[(5, 20)]);
        assert_eq!(
            score1_transcripts(&a, &b, 0.0),
            score1_transcripts(&b, &a, 0.0)
        );
        assert_eq!(
            score2_transcripts(&a, &b, 0.0),
            score2_transcripts(&b, &a, 0.0)
        );
    }

    #[test]
    fn score2_transcripts_penalize_containment_by_length_ratio() {
        let a = make_tx("a", &[(0, 10)]);
        let b = make_tx("b", &[(2, 8)]);

        assert_eq!(score2_transcripts(&a, &b, 0.0), 0.6);
        assert_eq!(score2_transcripts(&b, &a, 0.0), 0.6);
    }
}
