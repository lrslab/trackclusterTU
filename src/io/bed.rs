//! Streaming BED6/BED12 readers and writers.
//!
//! Coordinates are represented internally and on disk as 0-based, half-open
//! intervals. Readers skip blank lines and lines whose first non-whitespace
//! character is `#`.

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use thiserror::Error;

use crate::model::{
    Bed12Attrs, Coord, Interval, IntervalError, Strand, Transcript, TranscriptError,
};

/// Error produced while opening, parsing, or writing a BED file.
#[derive(Error, Debug)]
pub enum BedError {
    /// An input file could not be opened or read.
    #[error("I/O error reading {path:?}: {source}")]
    IoRead {
        /// Path being read.
        path: PathBuf,
        /// Underlying filesystem error.
        source: std::io::Error,
    },

    /// A non-comment input row was not valid BED6 or BED12.
    #[error("{path:?}:{line}: {source}")]
    Parse {
        /// Path containing the invalid row.
        path: PathBuf,
        /// One-based source line number.
        line: usize,
        /// Row-level parse error.
        source: BedParseError,
    },

    /// An output file could not be created or written.
    #[error("I/O error writing {path:?}: {source}")]
    IoWrite {
        /// Path being written.
        path: PathBuf,
        /// Underlying filesystem error.
        source: std::io::Error,
    },
}

/// Error produced while decoding one BED data row.
#[derive(Error, Debug)]
pub enum BedParseError {
    /// One physical BED row is not valid UTF-8.
    #[error("invalid UTF-8: {message}")]
    InvalidUtf8 {
        /// UTF-8 decoding error reported by the buffered line reader.
        message: String,
    },

    /// A BED6 row contains fewer than the six required columns.
    #[error("expected at least 6 columns, got {got}")]
    TooFewColumnsBed6 {
        /// Number of columns observed.
        got: usize,
    },

    /// A BED12 row contains fewer than the twelve required columns.
    #[error("expected at least 12 columns, got {got}")]
    TooFewColumns {
        /// Number of columns observed.
        got: usize,
    },

    /// A numeric field contains a value that cannot be parsed as a nonnegative integer.
    #[error("invalid integer for {field}: {value:?}")]
    InvalidInt {
        /// BED field name.
        field: &'static str,
        /// Invalid source text.
        value: String,
    },

    /// The BED12 block count is not a valid nonnegative integer.
    #[error("invalid blockCount: {value:?}")]
    InvalidBlockCount {
        /// Invalid `blockCount` text.
        value: String,
    },

    /// `blockCount` differs from the number of parsed sizes or starts.
    #[error("blockCount {block_count} does not match {field_name} length {list_len}")]
    BlockListLengthMismatch {
        /// Declared `blockCount`.
        block_count: usize,
        /// Mismatched list field (`blockSizes` or `blockStarts`).
        field_name: &'static str,
        /// Number of entries actually parsed.
        list_len: usize,
    },

    /// Adding a block start or size would overflow the 32-bit coordinate domain.
    #[error("block start+size overflows u32")]
    BlockOverflow,

    /// The strand column is invalid.
    #[error(transparent)]
    Strand(
        /// Underlying strand parse failure.
        #[from]
        crate::model::StrandParseError,
    ),

    /// The row describes an invalid half-open interval.
    #[error(transparent)]
    Interval(
        /// Underlying interval validation failure.
        #[from]
        IntervalError,
    ),

    /// The BED12 fields violate a transcript invariant.
    #[error(transparent)]
    Transcript(
        /// Underlying transcript validation failure.
        #[from]
        TranscriptError,
    ),
}

/// One BED6 row plus any trailing, application-specific columns.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Bed6Record {
    /// Reference sequence name (BED `chrom`).
    pub chrom: String,
    /// Zero-based inclusive start coordinate.
    pub start: Coord,
    /// Zero-based exclusive end coordinate.
    pub end: Coord,
    /// Record identifier.
    pub name: String,
    /// BED score, parsed as a nonnegative integer.
    pub score: u32,
    /// Parsed BED strand.
    pub strand: Strand,
    /// Columns after the standard six, preserved in source order.
    pub extra_fields: Vec<String>,
}

/// Streaming iterator over validated BED6 rows.
///
/// Each item retains path and one-based line context in its [`BedError`]. Blank
/// and comment lines are skipped.
pub struct Bed6Reader<R: BufRead> {
    inner: BedLineReader<R>,
}

impl<R: BufRead> Bed6Reader<R> {
    /// Return the one-based physical line number of the most recently yielded row.
    pub const fn line_number(&self) -> usize {
        self.inner.line_number
    }
}

struct BedLineReader<R: BufRead> {
    path: PathBuf,
    reader: R,
    buffer: Vec<u8>,
    line_number: usize,
}

impl<R: BufRead> BedLineReader<R> {
    fn new(path: PathBuf, reader: R) -> Self {
        Self {
            path,
            reader,
            buffer: Vec::new(),
            line_number: 0,
        }
    }

    fn next_data_line(&mut self) -> Option<Result<String, BedError>> {
        loop {
            self.buffer.clear();
            let bytes_read = match self.reader.read_until(b'\n', &mut self.buffer) {
                Ok(0) => return None,
                Ok(bytes_read) => bytes_read,
                Err(source) => {
                    return Some(Err(BedError::IoRead {
                        path: self.path.clone(),
                        source,
                    }))
                }
            };
            debug_assert!(bytes_read > 0);
            self.line_number += 1;

            let bytes = strip_line_ending(&self.buffer);
            let first_non_ascii_whitespace =
                bytes.iter().position(|byte| !byte.is_ascii_whitespace());
            if first_non_ascii_whitespace.is_none()
                || first_non_ascii_whitespace.is_some_and(|index| bytes[index] == b'#')
            {
                continue;
            }

            let decoded = match std::str::from_utf8(bytes) {
                Ok(decoded) => decoded,
                Err(source) => {
                    return Some(Err(BedError::Parse {
                        path: self.path.clone(),
                        line: self.line_number,
                        source: BedParseError::InvalidUtf8 {
                            message: source.to_string(),
                        },
                    }))
                }
            };
            if decoded.trim().is_empty() || decoded.trim_start().starts_with('#') {
                continue;
            }
            return Some(Ok(decoded.trim_end().to_owned()));
        }
    }
}

fn strip_line_ending(mut line: &[u8]) -> &[u8] {
    if line.last() == Some(&b'\n') {
        line = &line[..line.len() - 1];
    }
    if line.last() == Some(&b'\r') {
        line = &line[..line.len() - 1];
    }
    line
}

/// Open a BED6 file and return a streaming reader.
///
/// The file is opened immediately; individual rows are parsed as the returned
/// iterator advances. Rows may contain additional columns, available through
/// [`Bed6Record::extra_fields`].
///
/// # Errors
///
/// Returns [`BedError::IoRead`] if `path` cannot be opened. Row-level I/O and
/// parse failures are returned by the iterator.
pub fn read_bed6<P: AsRef<Path>>(path: P) -> Result<Bed6Reader<BufReader<File>>, BedError> {
    let path = path.as_ref().to_path_buf();
    let file = File::open(&path).map_err(|source| BedError::IoRead {
        path: path.clone(),
        source,
    })?;
    let reader = BufReader::new(file);
    Ok(Bed6Reader {
        inner: BedLineReader::new(path, reader),
    })
}

impl<R: BufRead> Iterator for Bed6Reader<R> {
    type Item = Result<Bed6Record, BedError>;

    fn next(&mut self) -> Option<Self::Item> {
        let line = match self.inner.next_data_line()? {
            Ok(line) => line,
            Err(error) => return Some(Err(error)),
        };
        Some(parse_bed6_line(&line).map_err(|source| BedError::Parse {
            path: self.inner.path.clone(),
            line: self.inner.line_number,
            source,
        }))
    }
}

/// Streaming iterator over validated BED12 transcript rows.
///
/// Every row is converted into a checked [`Transcript`]. Blank and comment
/// lines are skipped.
pub struct Bed12Reader<R: BufRead> {
    inner: BedLineReader<R>,
}

impl<R: BufRead> Bed12Reader<R> {
    /// Return the one-based physical line number of the most recently yielded row.
    pub const fn line_number(&self) -> usize {
        self.inner.line_number
    }
}

/// Open a BED12 file and return a streaming transcript reader.
///
/// # Errors
///
/// Returns [`BedError::IoRead`] if `path` cannot be opened. Row-level I/O and
/// parse failures are returned by the iterator.
pub fn read_bed12<P: AsRef<Path>>(path: P) -> Result<Bed12Reader<BufReader<File>>, BedError> {
    let path = path.as_ref().to_path_buf();
    let file = File::open(&path).map_err(|source| BedError::IoRead {
        path: path.clone(),
        source,
    })?;
    let reader = BufReader::new(file);
    Ok(Bed12Reader {
        inner: BedLineReader::new(path, reader),
    })
}

impl<R: BufRead> Iterator for Bed12Reader<R> {
    type Item = Result<Transcript, BedError>;

    fn next(&mut self) -> Option<Self::Item> {
        let line = match self.inner.next_data_line()? {
            Ok(line) => line,
            Err(error) => return Some(Err(error)),
        };
        Some(parse_bed12_line(&line).map_err(|source| BedError::Parse {
            path: self.inner.path.clone(),
            line: self.inner.line_number,
            source,
        }))
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
    // BED scores are required integers. Invalid values fail the record instead
    // of being silently rewritten to zero.
    let score = parse_u32("score", fields[4])?;
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
    let score = parse_u32("score", fields[4])?;
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

/// Create `path` and write transcripts as normalized BED12 rows.
///
/// Block sizes and relative starts are derived from each transcript's exon
/// intervals; trailing extra fields are retained.
///
/// # Errors
///
/// Returns [`BedError::IoWrite`] if the file cannot be created, written, or
/// flushed.
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
    write_bed12_to_writer(&mut writer, transcripts).map_err(|source| BedError::IoWrite {
        path: path.clone(),
        source,
    })?;
    writer
        .flush()
        .map_err(|source| BedError::IoWrite { path, source })?;
    Ok(())
}

/// Write transcripts as normalized BED12 rows to an existing byte writer.
///
/// # Errors
///
/// Returns the first I/O error reported by `writer`.
pub fn write_bed12_to_writer<'a, W, I>(writer: &mut W, transcripts: I) -> Result<(), std::io::Error>
where
    W: Write,
    I: IntoIterator<Item = &'a Transcript>,
{
    for transcript in transcripts {
        let block_count = transcript.exons().len();
        let mut block_sizes = String::new();
        let mut block_starts = String::new();

        for exon in transcript.exons() {
            let size = exon.len();
            let rel_start = exon
                .start()
                .get()
                .saturating_sub(transcript.tx_start().get());
            block_sizes.push_str(&format!("{size},"));
            block_starts.push_str(&format!("{rel_start},"));
        }

        write!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            transcript.chrom(),
            transcript.tx_start().get(),
            transcript.tx_end().get(),
            transcript.name(),
            transcript.score(),
            transcript.strand().as_char(),
            transcript.thick_start().get(),
            transcript.thick_end().get(),
            transcript.item_rgb(),
            block_count,
            block_sizes,
            block_starts
        )?;

        for extra in transcript.extra_fields() {
            write!(writer, "\t{extra}")?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

/// Create `path` and write BED6 records, including trailing extra fields.
///
/// # Errors
///
/// Returns [`BedError::IoWrite`] if the file cannot be created, written, or
/// flushed.
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
    write_bed6_to_writer(&mut writer, records).map_err(|source| BedError::IoWrite {
        path: path.clone(),
        source,
    })?;
    writer
        .flush()
        .map_err(|source| BedError::IoWrite { path, source })?;
    Ok(())
}

/// Write BED6 records, including trailing fields, to an existing byte writer.
///
/// This is useful for streaming into an already-open file or an in-memory
/// buffer.
///
/// # Examples
///
/// ```
/// use trackclustertu::io::bed::{write_bed6_to_writer, Bed6Record};
/// use trackclustertu::model::{Coord, Strand};
///
/// let record = Bed6Record {
///     chrom: "chr1".into(),
///     start: Coord::new(10),
///     end: Coord::new(20),
///     name: "read-1".into(),
///     score: 42,
///     strand: Strand::Plus,
///     extra_fields: vec!["full_length".into()],
/// };
/// let mut output = Vec::new();
/// write_bed6_to_writer(&mut output, [&record])?;
/// assert_eq!(
///     String::from_utf8(output)?,
///     "chr1\t10\t20\tread-1\t42\t+\tfull_length\n",
/// );
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
///
/// # Errors
///
/// Returns the first I/O error reported by `writer`.
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

    #[test]
    fn bed6_invalid_score_is_not_silently_normalized() {
        let error = parse_bed6_line("chr1\t0\t10\tread1\tnot-a-score\t+").unwrap_err();
        assert!(matches!(
            error,
            BedParseError::InvalidInt { field: "score", .. }
        ));
    }

    #[test]
    fn bed12_invalid_score_is_not_silently_normalized() {
        let error =
            parse_bed12_line("chr1\t0\t10\ttx1\tnot-a-score\t+\t0\t10\t0\t1\t10,\t0,").unwrap_err();
        assert!(matches!(
            error,
            BedParseError::InvalidInt { field: "score", .. }
        ));
    }

    #[test]
    fn bed_reader_reports_invalid_score_line_number() {
        let input = b"# header\nchr1\t0\t10\tread1\tnot-a-score\t+\n";
        let mut reader = Bed6Reader {
            inner: BedLineReader::new(PathBuf::from("reads.bed"), BufReader::new(&input[..])),
        };

        let error = reader.next().unwrap().unwrap_err();
        assert!(matches!(
            error,
            BedError::Parse {
                line: 2,
                source: BedParseError::InvalidInt { field: "score", .. },
                ..
            }
        ));
        assert_eq!(reader.line_number(), 2);
    }

    #[test]
    fn bed6_reader_treats_invalid_utf8_as_a_recoverable_row_error() {
        let mut input = b"chr1\t0\t10\tgood1\t0\t+\nchr1\t20\t30\t".to_vec();
        input.push(0xff);
        input.extend_from_slice(b"\t0\t+\nchr1\t40\t50\tgood2\t0\t+\n");
        let mut reader = Bed6Reader {
            inner: BedLineReader::new(PathBuf::from("reads.bed"), BufReader::new(input.as_slice())),
        };

        assert_eq!(reader.next().unwrap().unwrap().name, "good1");
        assert_eq!(reader.line_number(), 1);
        assert!(matches!(
            reader.next().unwrap().unwrap_err(),
            BedError::Parse {
                line: 2,
                source: BedParseError::InvalidUtf8 { .. },
                ..
            }
        ));
        assert_eq!(reader.line_number(), 2);
        assert_eq!(reader.next().unwrap().unwrap().name, "good2");
        assert_eq!(reader.line_number(), 3);
        assert!(reader.next().is_none());
    }

    #[test]
    fn bed6_reader_skips_invalid_utf8_comment_before_decoding() {
        let mut input = b" \t# metadata ".to_vec();
        input.push(0xff);
        input.extend_from_slice(b"\nchr1\t40\t50\tgood\t0\t+\n");
        let mut reader = Bed6Reader {
            inner: BedLineReader::new(PathBuf::from("reads.bed"), BufReader::new(input.as_slice())),
        };

        assert_eq!(reader.next().unwrap().unwrap().name, "good");
        assert_eq!(reader.line_number(), 2);
        assert!(reader.next().is_none());
    }

    #[test]
    fn bed12_reader_treats_invalid_utf8_as_a_recoverable_row_error() {
        let mut input = b"chr1\t0\t10\tgood1\t0\t+\t0\t10\t0\t1\t10,\t0,\nchr1\t20\t30\t".to_vec();
        input.push(0xff);
        input.extend_from_slice(
            b"\t0\t+\t20\t30\t0\t1\t10,\t0,\nchr1\t40\t50\tgood2\t0\t+\t40\t50\t0\t1\t10,\t0,\n",
        );
        let mut reader = Bed12Reader {
            inner: BedLineReader::new(
                PathBuf::from("reads.bed12"),
                BufReader::new(input.as_slice()),
            ),
        };

        assert_eq!(reader.next().unwrap().unwrap().name(), "good1");
        assert_eq!(reader.line_number(), 1);
        assert!(matches!(
            reader.next().unwrap().unwrap_err(),
            BedError::Parse {
                line: 2,
                source: BedParseError::InvalidUtf8 { .. },
                ..
            }
        ));
        assert_eq!(reader.next().unwrap().unwrap().name(), "good2");
        assert_eq!(reader.line_number(), 3);
        assert!(reader.next().is_none());
    }

    #[test]
    fn bed12_reader_skips_invalid_utf8_comment_before_decoding() {
        let mut input = b"# metadata ".to_vec();
        input.push(0xff);
        input.extend_from_slice(b"\nchr1\t40\t50\tgood\t0\t+\t40\t50\t0\t1\t10,\t0,\n");
        let mut reader = Bed12Reader {
            inner: BedLineReader::new(
                PathBuf::from("reads.bed12"),
                BufReader::new(input.as_slice()),
            ),
        };

        assert_eq!(reader.next().unwrap().unwrap().name(), "good");
        assert_eq!(reader.line_number(), 2);
        assert!(reader.next().is_none());
    }

    #[test]
    fn bed12_rejects_empty_block() {
        let error = parse_bed12_line("chr1\t0\t10\ttx1\t0\t+\t0\t10\t0\t1\t0,\t0,").unwrap_err();
        assert!(matches!(
            error,
            BedParseError::Transcript(TranscriptError::EmptyExon { .. })
        ));
    }

    #[test]
    fn bed12_rejects_overlapping_blocks() {
        let error =
            parse_bed12_line("chr1\t0\t20\ttx1\t0\t+\t0\t20\t0\t2\t10,10,\t0,5,").unwrap_err();
        assert!(matches!(
            error,
            BedParseError::Transcript(TranscriptError::OverlappingExons { .. })
        ));
    }

    #[test]
    fn bed12_rejects_block_past_transcript_end() {
        let error = parse_bed12_line("chr1\t0\t10\ttx1\t0\t+\t0\t10\t0\t1\t6,\t5,").unwrap_err();
        assert!(matches!(
            error,
            BedParseError::Transcript(TranscriptError::ExonOutsideSpan { .. })
        ));
    }
}
