use std::fs::File;
use std::io::{BufRead, BufReader, Lines, Write};
use std::path::{Path, PathBuf};

use thiserror::Error;

use crate::model::{
    interval::IntervalError, Bed12Attrs, Coord, Interval, Strand, Transcript, TranscriptError,
};

#[derive(Error, Debug)]
pub enum BedError {
    #[error("I/O error reading {path:?}: {source}")]
    IoRead {
        path: PathBuf,
        source: std::io::Error,
    },

    #[error("{path:?}:{line}: {source}")]
    Parse {
        path: PathBuf,
        line: usize,
        source: BedParseError,
    },

    #[error("I/O error writing {path:?}: {source}")]
    IoWrite {
        path: PathBuf,
        source: std::io::Error,
    },
}

#[derive(Error, Debug)]
pub enum BedParseError {
    #[error("expected at least 6 columns, got {got}")]
    TooFewColumnsBed6 { got: usize },

    #[error("expected at least 12 columns, got {got}")]
    TooFewColumns { got: usize },

    #[error("invalid integer for {field}: {value:?}")]
    InvalidInt { field: &'static str, value: String },

    #[error("invalid blockCount: {value:?}")]
    InvalidBlockCount { value: String },

    #[error("blockCount {block_count} does not match {field_name} length {list_len}")]
    BlockListLengthMismatch {
        block_count: usize,
        field_name: &'static str,
        list_len: usize,
    },

    #[error("block start+size overflows u32")]
    BlockOverflow,

    #[error(transparent)]
    Strand(#[from] crate::model::strand::StrandParseError),

    #[error(transparent)]
    Interval(#[from] IntervalError),

    #[error(transparent)]
    Transcript(#[from] TranscriptError),
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Bed6Record {
    pub chrom: String,
    pub start: Coord,
    pub end: Coord,
    pub name: String,
    pub score: u32,
    pub strand: Strand,
    pub extra_fields: Vec<String>,
}

pub struct Bed6Reader<R: BufRead> {
    path: PathBuf,
    line_number: usize,
    lines: Lines<R>,
}

pub fn read_bed6<P: AsRef<Path>>(path: P) -> Result<Bed6Reader<BufReader<File>>, BedError> {
    let path = path.as_ref().to_path_buf();
    let file = File::open(&path).map_err(|source| BedError::IoRead {
        path: path.clone(),
        source,
    })?;
    let reader = BufReader::new(file);
    Ok(Bed6Reader {
        path,
        line_number: 0,
        lines: reader.lines(),
    })
}

impl<R: BufRead> Iterator for Bed6Reader<R> {
    type Item = Result<Bed6Record, BedError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.line_number += 1;
            let line = match self.lines.next()? {
                Ok(line) => line,
                Err(source) => {
                    return Some(Err(BedError::IoRead {
                        path: self.path.clone(),
                        source,
                    }))
                }
            };

            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let record = match parse_bed6_line(line) {
                Ok(record) => record,
                Err(source) => {
                    return Some(Err(BedError::Parse {
                        path: self.path.clone(),
                        line: self.line_number,
                        source,
                    }))
                }
            };

            return Some(Ok(record));
        }
    }
}

pub struct Bed12Reader<R: BufRead> {
    path: PathBuf,
    line_number: usize,
    lines: Lines<R>,
}

pub fn read_bed12<P: AsRef<Path>>(path: P) -> Result<Bed12Reader<BufReader<File>>, BedError> {
    let path = path.as_ref().to_path_buf();
    let file = File::open(&path).map_err(|source| BedError::IoRead {
        path: path.clone(),
        source,
    })?;
    let reader = BufReader::new(file);
    Ok(Bed12Reader {
        path,
        line_number: 0,
        lines: reader.lines(),
    })
}

impl<R: BufRead> Iterator for Bed12Reader<R> {
    type Item = Result<Transcript, BedError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.line_number += 1;
            let line = match self.lines.next()? {
                Ok(line) => line,
                Err(source) => {
                    return Some(Err(BedError::IoRead {
                        path: self.path.clone(),
                        source,
                    }))
                }
            };

            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let record = match parse_bed12_line(line) {
                Ok(record) => record,
                Err(source) => {
                    return Some(Err(BedError::Parse {
                        path: self.path.clone(),
                        line: self.line_number,
                        source,
                    }))
                }
            };

            return Some(Ok(record));
        }
    }
}

fn parse_u32(field: &'static str, value: &str) -> Result<u32, BedParseError> {
    value.parse::<u32>().map_err(|_| BedParseError::InvalidInt {
        field,
        value: value.to_owned(),
    })
}

fn parse_usize(field: &'static str, value: &str) -> Result<usize, BedParseError> {
    value
        .parse::<usize>()
        .map_err(|_| BedParseError::InvalidInt {
            field,
            value: value.to_owned(),
        })
}

fn parse_comma_u32_list(value: &str) -> Result<Vec<u32>, BedParseError> {
    value
        .split(',')
        .filter(|token| !token.is_empty())
        .map(|token| parse_u32("blockList", token))
        .collect()
}

fn parse_bed6_line(line: &str) -> Result<Bed6Record, BedParseError> {
    let mut fields: Vec<&str> = line.split('\t').collect();
    if fields.len() == 1 {
        fields = line.split_whitespace().collect();
    }

    if fields.len() < 6 {
        return Err(BedParseError::TooFewColumnsBed6 { got: fields.len() });
    }

    let chrom = fields[0].to_owned();
    let start = Coord::new(parse_u32("chromStart", fields[1])?);
    let end = Coord::new(parse_u32("chromEnd", fields[2])?);
    Interval::new(start, end)?;

    let name = fields[3].to_owned();
    let score = parse_u32("score", fields[4]).unwrap_or(0);
    let strand = Strand::try_from(fields[5])?;

    let extra_fields = fields[6..]
        .iter()
        .map(|value| (*value).to_owned())
        .collect();

    Ok(Bed6Record {
        chrom,
        start,
        end,
        name,
        score,
        strand,
        extra_fields,
    })
}

fn parse_bed12_line(line: &str) -> Result<Transcript, BedParseError> {
    let mut fields: Vec<&str> = line.split('\t').collect();
    if fields.len() == 1 {
        fields = line.split_whitespace().collect();
    }

    if fields.len() < 12 {
        return Err(BedParseError::TooFewColumns { got: fields.len() });
    }

    let chrom = fields[0].to_owned();
    let tx_start = Coord::new(parse_u32("chromStart", fields[1])?);
    let tx_end = Coord::new(parse_u32("chromEnd", fields[2])?);
    let name = fields[3].to_owned();
    let score = parse_u32("score", fields[4]).unwrap_or(0);
    let strand = Strand::try_from(fields[5])?;
    let thick_start = Coord::new(parse_u32("thickStart", fields[6])?);
    let thick_end = Coord::new(parse_u32("thickEnd", fields[7])?);
    let item_rgb = fields[8].to_owned();
    let block_count =
        parse_usize("blockCount", fields[9]).map_err(|_| BedParseError::InvalidBlockCount {
            value: fields[9].to_owned(),
        })?;

    let block_sizes = parse_comma_u32_list(fields[10])?;
    if block_sizes.len() != block_count {
        return Err(BedParseError::BlockListLengthMismatch {
            block_count,
            field_name: "blockSizes",
            list_len: block_sizes.len(),
        });
    }

    let block_starts = parse_comma_u32_list(fields[11])?;
    if block_starts.len() != block_count {
        return Err(BedParseError::BlockListLengthMismatch {
            block_count,
            field_name: "blockStarts",
            list_len: block_starts.len(),
        });
    }

    let mut exons = Vec::with_capacity(block_count);
    for i in 0..block_count {
        let rel_start = block_starts[i];
        let block_size = block_sizes[i];

        let exon_start_u32 = tx_start
            .get()
            .checked_add(rel_start)
            .ok_or(BedParseError::BlockOverflow)?;
        let exon_end_u32 = exon_start_u32
            .checked_add(block_size)
            .ok_or(BedParseError::BlockOverflow)?;

        // Some legacy BEDs contain block sizes that extend past chromEnd; clamp to transcript end.
        let exon_end_u32 = exon_end_u32.min(tx_end.get());

        let exon = Interval::new(Coord::new(exon_start_u32), Coord::new(exon_end_u32))
            .map_err(|_| BedParseError::BlockOverflow)?;
        exons.push(exon);
    }

    let extra_fields = fields[12..]
        .iter()
        .map(|value| (*value).to_owned())
        .collect();

    Ok(Transcript::new(
        chrom,
        strand,
        tx_start,
        tx_end,
        name,
        exons,
        Bed12Attrs {
            score,
            thick_start,
            thick_end,
            item_rgb,
            extra_fields,
        },
    )?)
}

pub fn write_bed12<'a, P, I>(path: P, transcripts: I) -> Result<(), BedError>
where
    P: AsRef<Path>,
    I: IntoIterator<Item = &'a Transcript>,
{
    let path = path.as_ref().to_path_buf();
    let file = File::create(&path).map_err(|source| BedError::IoWrite {
        path: path.clone(),
        source,
    })?;
    let mut writer = std::io::BufWriter::new(file);
    write_bed12_to_writer(&mut writer, transcripts)
        .map_err(|source| BedError::IoWrite { path, source })?;
    Ok(())
}

pub fn write_bed12_to_writer<'a, W, I>(writer: &mut W, transcripts: I) -> Result<(), std::io::Error>
where
    W: Write,
    I: IntoIterator<Item = &'a Transcript>,
{
    for transcript in transcripts {
        let block_count = transcript.exons.len();
        let mut block_sizes = String::new();
        let mut block_starts = String::new();

        for exon in &transcript.exons {
            let size = exon.len();
            let rel_start = exon.start.get().saturating_sub(transcript.tx_start.get());
            block_sizes.push_str(&format!("{size},"));
            block_starts.push_str(&format!("{rel_start},"));
        }

        write!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            transcript.chrom,
            transcript.tx_start.get(),
            transcript.tx_end.get(),
            transcript.name,
            transcript.score,
            transcript.strand.as_char(),
            transcript.thick_start.get(),
            transcript.thick_end.get(),
            transcript.item_rgb,
            block_count,
            block_sizes,
            block_starts
        )?;

        for extra in &transcript.extra_fields {
            write!(writer, "\t{extra}")?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

pub fn write_bed6<'a, P, I>(path: P, records: I) -> Result<(), BedError>
where
    P: AsRef<Path>,
    I: IntoIterator<Item = &'a Bed6Record>,
{
    let path = path.as_ref().to_path_buf();
    let file = File::create(&path).map_err(|source| BedError::IoWrite {
        path: path.clone(),
        source,
    })?;
    let mut writer = std::io::BufWriter::new(file);
    write_bed6_to_writer(&mut writer, records)
        .map_err(|source| BedError::IoWrite { path, source })?;
    Ok(())
}

pub fn write_bed6_to_writer<'a, W, I>(writer: &mut W, records: I) -> Result<(), std::io::Error>
where
    W: Write,
    I: IntoIterator<Item = &'a Bed6Record>,
{
    for record in records {
        write!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}",
            record.chrom,
            record.start.get(),
            record.end.get(),
            record.name,
            record.score,
            record.strand.as_char()
        )?;

        for extra in &record.extra_fields {
            write!(writer, "\t{extra}")?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bed12_roundtrip_normalized() {
        let input = "chr1\t100\t200\ttx1\t0\t+\t100\t200\t0\t2\t50,30,\t0,70,\n";
        let transcript = parse_bed12_line(input.trim()).unwrap();

        let mut buffer = Vec::new();
        write_bed12_to_writer(&mut buffer, [&transcript]).unwrap();

        let output = String::from_utf8(buffer).unwrap();
        let reparsed = parse_bed12_line(output.trim()).unwrap();
        assert_eq!(transcript, reparsed);
    }

    #[test]
    fn bed6_roundtrip_normalized() {
        let input = "chr1\t0\t10\tread1\t42\t+\textra1\textra2\n";
        let record = parse_bed6_line(input.trim()).unwrap();

        let mut buffer = Vec::new();
        write_bed6_to_writer(&mut buffer, [&record]).unwrap();

        let output = String::from_utf8(buffer).unwrap();
        let reparsed = parse_bed6_line(output.trim()).unwrap();
        assert_eq!(record, reparsed);
    }
}
