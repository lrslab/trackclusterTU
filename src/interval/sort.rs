use crate::model::{Strand, Transcript};

fn strand_rank(strand: Strand) -> u8 {
    match strand {
        Strand::Plus => 0,
        Strand::Minus => 1,
        Strand::Unknown => 2,
    }
}

pub fn sort_by_coord(transcripts: &mut [Transcript]) {
    transcripts.sort_by(|left, right| {
        left.chrom
            .cmp(&right.chrom)
            .then_with(|| left.tx_start.cmp(&right.tx_start))
            .then_with(|| left.tx_end.cmp(&right.tx_end))
            .then_with(|| strand_rank(left.strand).cmp(&strand_rank(right.strand)))
    });
}

#[cfg(test)]
mod tests {
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
    fn sort_by_coord_is_stable() {
        let a = make_tx("chr1", Strand::Plus, 10, 20, "a");
        let b = make_tx("chr1", Strand::Plus, 10, 20, "b");

        let mut transcripts = vec![a.clone(), b.clone()];
        sort_by_coord(&mut transcripts);
        assert_eq!(transcripts[0].name, "a");
        assert_eq!(transcripts[1].name, "b");
    }
}
