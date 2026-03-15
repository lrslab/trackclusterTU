use crate::model::Transcript;

fn span_len(tx: &Transcript) -> u32 {
    tx.tx_end.get() - tx.tx_start.get()
}

pub fn better_representative(candidate: &Transcript, current: &Transcript) -> bool {
    let candidate_len = span_len(candidate);
    let current_len = span_len(current);

    if candidate_len != current_len {
        return candidate_len > current_len;
    }

    candidate.name < current.name
}

#[cfg(test)]
mod tests {
    use crate::model::{Bed12Attrs, Coord, Interval, Strand, Transcript};

    use super::*;

    fn make_tx(start: u32, end: u32, name: &str) -> Transcript {
        Transcript::new(
            "chr1".to_owned(),
            Strand::Plus,
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
    fn prefers_longer_transcript() {
        let current = make_tx(0, 10, "b");
        let candidate = make_tx(0, 20, "a");
        assert!(better_representative(&candidate, &current));
    }

    #[test]
    fn tie_breaks_by_name() {
        let current = make_tx(0, 10, "b");
        let candidate = make_tx(0, 10, "a");
        assert!(better_representative(&candidate, &current));
    }
}
