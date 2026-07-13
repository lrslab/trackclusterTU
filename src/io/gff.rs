//! GFF3 gene reader and BED6 gene-track writer.
//!
//! GFF3's 1-based inclusive coordinates are converted to the crate's 0-based,
//! half-open [`crate::model::Interval`] representation. Percent escapes in attribute keys and
//! values are decoded without applying HTML form rules (so `+` remains `+`).

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use thiserror::Error;

use crate::model::{Coord, Interval, IntervalError, Strand, StrandParseError};

/// A gene feature extracted from a GFF3 file.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeneRecord {
    /// Reference sequence name from GFF3 column 1.
    pub contig: String,
    /// Gene strand from GFF3 column 7.
    pub strand: Strand,
    /// Zero-based, half-open gene interval.
    pub interval: Interval,
    /// Preferred gene identifier.
    ///
    /// Attributes are considered in the order `Name`, `gene`, `locus_tag`, and
    /// `ID`; a deterministic `gene_N` fallback is used when none is present.
    pub id: String,
    /// BED-facing feature category derived from `gene_biotype`.
    pub feature_kind: String,
}

/// Error produced while opening or parsing a GFF3 file.
#[derive(Error, Debug)]
pub enum GffError {
    /// The input file could not be opened or read.
    #[error("I/O error reading {path:?}: {source}")]
    IoRead {
        /// Path being read.
        path: PathBuf,
        /// Underlying filesystem error.
        source: std::io::Error,
    },

    /// A data row was not valid GFF3.
    #[error("{path:?}:{line}: {source}")]
    Parse {
        /// Path containing the invalid row.
        path: PathBuf,
        /// One-based source line number.
        line: usize,
        /// Row-level parse error.
        source: GffParseError,
    },
}

/// Error produced while decoding one GFF3 data row or attribute value.
#[derive(Error, Debug)]
pub enum GffParseError {
    /// A row does not contain the nine tab-separated GFF3 columns.
    #[error("expected exactly 9 tab-separated columns, got {got}")]
    WrongColumnCount {
        /// Number of columns observed.
        got: usize,
    },

    /// A coordinate field is not a valid nonnegative integer.
    #[error("invalid integer for {field}: {value:?}")]
    InvalidInt {
        /// GFF3 field name.
        field: &'static str,
        /// Invalid source text.
        value: String,
    },

    /// The 1-based inclusive start coordinate is zero.
    #[error("expected 1-based GFF start > 0")]
    ZeroStart,

    /// An attribute contains an incomplete or non-hexadecimal percent escape.
    #[error("invalid percent encoding in GFF3 attribute {value:?} at byte {offset}")]
    InvalidPercentEncoding {
        /// Full encoded attribute key or value.
        value: String,
        /// Zero-based byte offset of the invalid `%` escape.
        offset: usize,
    },

    /// Percent-decoded attribute bytes are not valid UTF-8.
    #[error("percent-decoded GFF3 attribute is not valid UTF-8: {value:?}")]
    InvalidAttributeUtf8 {
        /// Original encoded attribute key or value.
        value: String,
    },

    /// The strand column is invalid.
    #[error(transparent)]
    Strand(
        /// Underlying strand parse failure.
        #[from]
        StrandParseError,
    ),

    /// The converted coordinates do not form a valid interval.
    #[error(transparent)]
    Interval(
        /// Underlying interval validation failure.
        #[from]
        IntervalError,
    ),
}

/// Read `gene` features from a GFF3 file.
///
/// Blank/comment lines are ignored. Every data row must contain exactly nine
/// tab-separated columns; non-`gene` rows are otherwise skipped. Gene
/// coordinates are converted from 1-based inclusive to 0-based half-open.
///
/// # Errors
///
/// Returns [`GffError::IoRead`] for filesystem failures and [`GffError::Parse`]
/// with path and line context for malformed rows.
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
        if fields.len() != 9 {
            return Err(GffError::Parse {
                path: path.clone(),
                line: line_number,
                source: GffParseError::WrongColumnCount { got: fields.len() },
            });
        }
        if fields[2] != "gene" {
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
        let attributes = parse_gff_attributes(fields[8]).map_err(|source| GffError::Parse {
            path: path.clone(),
            line: line_number,
            source,
        })?;
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

/// Create a BED6 file from GFF-derived gene records.
///
/// The BED name is [`GeneRecord::id`], the score is zero, and the stored
/// zero-based half-open coordinates are written without conversion.
///
/// # Errors
///
/// Returns an I/O error if the output cannot be created, written, or flushed.
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

fn parse_gff_attributes(raw: &str) -> Result<HashMap<String, String>, GffParseError> {
    raw.split(';')
        .filter_map(|field| field.split_once('='))
        .map(|(key, value)| Ok((percent_decode(key)?, percent_decode(value)?)))
        .collect()
}

fn percent_decode(value: &str) -> Result<String, GffParseError> {
    let input = value.as_bytes();
    let mut decoded = Vec::with_capacity(input.len());
    let mut index = 0;

    while index < input.len() {
        if input[index] != b'%' {
            decoded.push(input[index]);
            index += 1;
            continue;
        }

        let Some(high) = input.get(index + 1).and_then(|byte| hex_value(*byte)) else {
            return Err(GffParseError::InvalidPercentEncoding {
                value: value.to_owned(),
                offset: index,
            });
        };
        let Some(low) = input.get(index + 2).and_then(|byte| hex_value(*byte)) else {
            return Err(GffParseError::InvalidPercentEncoding {
                value: value.to_owned(),
                offset: index,
            });
        };

        decoded.push((high << 4) | low);
        index += 3;
    }

    String::from_utf8(decoded).map_err(|_| GffParseError::InvalidAttributeUtf8 {
        value: value.to_owned(),
    })
}

fn hex_value(byte: u8) -> Option<u8> {
    match byte {
        b'0'..=b'9' => Some(byte - b'0'),
        b'a'..=b'f' => Some(byte - b'a' + 10),
        b'A'..=b'F' => Some(byte - b'A' + 10),
        _ => None,
    }
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
            gene.interval.start().get(),
            gene.interval.end().get(),
            gene.id,
            gene.strand.as_char()
        )?;
    }
    writer.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    struct FlushFailingWriter {
        bytes: Vec<u8>,
    }

    impl Write for FlushFailingWriter {
        fn write(&mut self, buffer: &[u8]) -> std::io::Result<usize> {
            self.bytes.extend_from_slice(buffer);
            Ok(buffer.len())
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Err(std::io::Error::new(
                std::io::ErrorKind::BrokenPipe,
                "deterministic flush failure",
            ))
        }
    }

    #[test]
    fn write_gene_bed6_surfaces_flush_errors() {
        let gene = GeneRecord {
            contig: "chr1".to_owned(),
            strand: Strand::Plus,
            interval: Interval::new(Coord::new(10), Coord::new(20)).unwrap(),
            id: "gene-1".to_owned(),
            feature_kind: "gene".to_owned(),
        };
        let mut writer = std::io::BufWriter::new(FlushFailingWriter { bytes: Vec::new() });

        let error = write_gene_bed6_to_writer(&mut writer, &[gene]).unwrap_err();

        assert_eq!(error.kind(), std::io::ErrorKind::BrokenPipe);
        assert_eq!(error.to_string(), "deterministic flush failure");
        assert_eq!(writer.get_ref().bytes, b"chr1\t10\t20\tgene-1\t0\t+\n");
    }

    #[test]
    fn gff_attributes_prefer_name() {
        let attrs = parse_gff_attributes("ID=gene-1;gene=thrL;Name=thrL;locus_tag=b0001").unwrap();
        assert_eq!(attrs.get("Name").map(String::as_str), Some("thrL"));
        assert_eq!(attrs.get("locus_tag").map(String::as_str), Some("b0001"));
    }

    #[test]
    fn gff_attributes_are_percent_decoded_without_form_url_rules() {
        let attrs = parse_gff_attributes("ID=gene%2D1;Name=alpha%20beta%3Bgamma%25+delta").unwrap();
        assert_eq!(attrs.get("ID").map(String::as_str), Some("gene-1"));
        assert_eq!(
            attrs.get("Name").map(String::as_str),
            Some("alpha beta;gamma%+delta")
        );
    }

    #[test]
    fn read_gff3_genes_rejects_malformed_non_gene_row() {
        let dir = std::env::temp_dir().join(format!(
            "trackclustertu_gff_malformed_unit_{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("genes.gff3");
        std::fs::write(&path, "chr1\tRefSeq\tCDS\t11\t50\t.\t+\t0\n").unwrap();

        let error = read_gff3_genes(&path).unwrap_err();
        assert!(matches!(
            error,
            GffError::Parse {
                line: 1,
                source: GffParseError::WrongColumnCount { got: 8 },
                ..
            }
        ));

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn read_gff3_genes_reports_invalid_attribute_encoding_with_line() {
        let dir = std::env::temp_dir().join(format!(
            "trackclustertu_gff_attribute_unit_{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("genes.gff3");
        std::fs::write(
            &path,
            "##gff-version 3\nchr1\tRefSeq\tgene\t11\t50\t.\t+\t.\tID=bad%ZZ\n",
        )
        .unwrap();

        let error = read_gff3_genes(&path).unwrap_err();
        assert!(matches!(
            error,
            GffError::Parse {
                line: 2,
                source: GffParseError::InvalidPercentEncoding { .. },
                ..
            }
        ));

        let _ = std::fs::remove_dir_all(&dir);
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
        assert_eq!(genes[0].interval.start().get(), 10);
        assert_eq!(genes[0].interval.end().get(), 50);
        assert_eq!(genes[0].id, "geneA");
        assert_eq!(genes[0].feature_kind, "mRNA");
        assert_eq!(genes[1].id, "b0002");
        assert_eq!(genes[1].feature_kind, "gene");

        let _ = std::fs::remove_dir_all(&dir);
    }
}
