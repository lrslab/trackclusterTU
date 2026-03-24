use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use thiserror::Error;

use crate::model::{interval::IntervalError, strand::StrandParseError, Coord, Interval, Strand};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeneRecord {
    pub contig: String,
    pub strand: Strand,
    pub interval: Interval,
    pub id: String,
    pub feature_kind: String,
}

#[derive(Error, Debug)]
pub enum GffError {
    #[error("I/O error reading {path:?}: {source}")]
    IoRead {
        path: PathBuf,
        source: std::io::Error,
    },

    #[error("{path:?}:{line}: {source}")]
    Parse {
        path: PathBuf,
        line: usize,
        source: GffParseError,
    },
}

#[derive(Error, Debug)]
pub enum GffParseError {
    #[error("invalid integer for {field}: {value:?}")]
    InvalidInt { field: &'static str, value: String },

    #[error("expected 1-based GFF start > 0")]
    ZeroStart,

    #[error(transparent)]
    Strand(#[from] StrandParseError),

    #[error(transparent)]
    Interval(#[from] IntervalError),
}

pub fn read_gff3_genes<P: AsRef<Path>>(path: P) -> Result<Vec<GeneRecord>, GffError> {
    let path = path.as_ref().to_path_buf();
    let file = File::open(&path).map_err(|source| GffError::IoRead {
        path: path.clone(),
        source,
    })?;
    let reader = BufReader::new(file);
    let mut genes = Vec::new();

    for (line_idx, line_result) in reader.lines().enumerate() {
        let line_number = line_idx + 1;
        let line = line_result.map_err(|source| GffError::IoRead {
            path: path.clone(),
            source,
        })?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 9 || fields[2] != "gene" {
            continue;
        }

        let start_1based = parse_u32("start", fields[3]).map_err(|source| GffError::Parse {
            path: path.clone(),
            line: line_number,
            source,
        })?;
        let end_1based = parse_u32("end", fields[4]).map_err(|source| GffError::Parse {
            path: path.clone(),
            line: line_number,
            source,
        })?;
        if start_1based == 0 {
            return Err(GffError::Parse {
                path: path.clone(),
                line: line_number,
                source: GffParseError::ZeroStart,
            });
        }

        let strand = Strand::try_from(fields[6]).map_err(|source| GffError::Parse {
            path: path.clone(),
            line: line_number,
            source: source.into(),
        })?;
        let interval = Interval::new(Coord::new(start_1based - 1), Coord::new(end_1based))
            .map_err(|source| GffError::Parse {
                path: path.clone(),
                line: line_number,
                source: source.into(),
            })?;
        let attributes = parse_gff_attributes(fields[8]);
        let gene_id = attributes
            .get("Name")
            .or_else(|| attributes.get("gene"))
            .or_else(|| attributes.get("locus_tag"))
            .or_else(|| attributes.get("ID"))
            .cloned()
            .unwrap_or_else(|| format!("gene_{}", genes.len() + 1));

        genes.push(GeneRecord {
            contig: fields[0].to_owned(),
            strand,
            interval,
            id: gene_id,
            feature_kind: gene_bed_kind(attributes.get("gene_biotype")),
        });
    }

    Ok(genes)
}

pub fn write_gene_bed6<P: AsRef<Path>>(
    path: P,
    genes: &[GeneRecord],
) -> Result<(), std::io::Error> {
    let path = path.as_ref();
    let file = File::create(path)?;
    let mut writer = std::io::BufWriter::new(file);
    write_gene_bed6_to_writer(&mut writer, genes)
}

fn parse_u32(field: &'static str, value: &str) -> Result<u32, GffParseError> {
    value.parse::<u32>().map_err(|_| GffParseError::InvalidInt {
        field,
        value: value.to_owned(),
    })
}

fn parse_gff_attributes(raw: &str) -> HashMap<String, String> {
    raw.split(';')
        .filter_map(|field| field.split_once('='))
        .map(|(key, value)| (key.to_owned(), value.to_owned()))
        .collect()
}

fn gene_bed_kind(gene_biotype: Option<&String>) -> String {
    match gene_biotype.map(String::as_str) {
        Some("protein_coding") => "mRNA".to_owned(),
        Some("ncRNA") => "ncRNA".to_owned(),
        Some("tRNA") => "tRNA".to_owned(),
        Some("rRNA") => "rRNA".to_owned(),
        Some(other) => other.to_owned(),
        None => "gene".to_owned(),
    }
}

fn write_gene_bed6_to_writer<W: Write>(
    writer: &mut W,
    genes: &[GeneRecord],
) -> Result<(), std::io::Error> {
    for gene in genes {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t0\t{}",
            gene.contig,
            gene.interval.start.get(),
            gene.interval.end.get(),
            gene.id,
            gene.strand.as_char()
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gff_attributes_prefer_name() {
        let attrs = parse_gff_attributes("ID=gene-1;gene=thrL;Name=thrL;locus_tag=b0001");
        assert_eq!(attrs.get("Name").map(String::as_str), Some("thrL"));
        assert_eq!(attrs.get("locus_tag").map(String::as_str), Some("b0001"));
    }

    #[test]
    fn read_gff3_genes_converts_gene_features_to_bed_coords() {
        let dir = std::env::temp_dir().join(format!(
            "trackclustertu_gff_unit_{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("genes.gff3");
        std::fs::write(
            &path,
            concat!(
                "##gff-version 3\n",
                "chr1\tRefSeq\tgene\t11\t50\t.\t+\t.\tID=id1;gene=geneA;gene_biotype=protein_coding\n",
                "chr1\tRefSeq\tCDS\t11\t50\t.\t+\t0\tParent=id1\n",
                "chr1\tRefSeq\tgene\t61\t100\t.\t-\t.\tID=id2;locus_tag=b0002\n",
            ),
        )
        .unwrap();

        let genes = read_gff3_genes(&path).unwrap();
        assert_eq!(genes.len(), 2);
        assert_eq!(genes[0].contig, "chr1");
        assert_eq!(genes[0].interval.start.get(), 10);
        assert_eq!(genes[0].interval.end.get(), 50);
        assert_eq!(genes[0].id, "geneA");
        assert_eq!(genes[0].feature_kind, "mRNA");
        assert_eq!(genes[1].id, "b0002");
        assert_eq!(genes[1].feature_kind, "gene");

        let _ = std::fs::remove_dir_all(&dir);
    }
}
