use std::path::Path;

pub(crate) fn convert_gff_to_bed(annotation_gff: &Path, out_bed: &Path) -> anyhow::Result<usize> {
    let genes = crate::io::gff::read_gff3_genes(annotation_gff)?;
    crate::io::gff::write_gene_bed6(out_bed, &genes)?;
    Ok(genes.len())
}
