//! Crate-private standards-compliant CSV/TSV records.
//!
//! BED and GFF intentionally do not use this layer: they are strict genomic
//! interchange formats with their own escaping and column rules.
//! An unquoted physical line whose first non-whitespace byte is `#` is reserved
//! for schema/comment metadata and is not returned as data. The writer quotes
//! fields containing `#`, so application data remains round-trippable.

use std::cell::RefCell;
use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::rc::Rc;

use csv::{ReaderBuilder, StringRecord, Terminator, WriterBuilder};
use thiserror::Error;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Delimiter {
    Comma,
    Tab,
}

impl Delimiter {
    const fn byte(self) -> u8 {
        match self {
            Self::Comma => b',',
            Self::Tab => b'\t',
        }
    }
}

#[derive(Debug, Error)]
pub(crate) enum DelimitedError {
    #[error("failed to open delimited table {path:?}: {source}")]
    Open {
        path: PathBuf,
        source: std::io::Error,
    },

    #[error("failed to read delimited table {path:?}{location}: {source}", location = format_location(*record, *line))]
    Read {
        path: PathBuf,
        record: Option<u64>,
        line: Option<u64>,
        source: csv::Error,
    },

    #[error("failed to create delimited table {path:?}: {source}")]
    Create {
        path: PathBuf,
        source: std::io::Error,
    },

    #[error("failed to write metadata for delimited table {path:?}: {source}")]
    MetadataWrite {
        path: PathBuf,
        source: std::io::Error,
    },

    #[error("failed to write record {record} to delimited table {path:?}: {source}")]
    Write {
        path: PathBuf,
        record: u64,
        source: csv::Error,
    },

    #[error("failed to flush delimited table {path:?}: {source}")]
    Flush {
        path: PathBuf,
        source: std::io::Error,
    },
}

fn format_location(record: Option<u64>, line: Option<u64>) -> String {
    match (record, line) {
        (Some(record), Some(line)) => format!(" at record {record}, line {line}"),
        (Some(record), None) => format!(" at record {record}"),
        (None, Some(line)) => format!(" at line {line}"),
        (None, None) => String::new(),
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct DelimitedRecord {
    fields: StringRecord,
    line_number: u64,
}

impl DelimitedRecord {
    pub(crate) fn fields(&self) -> &StringRecord {
        &self.fields
    }

    pub(crate) fn line_number(&self) -> u64 {
        self.line_number
    }
}

pub(crate) struct DelimitedReader<R: Read> {
    path: PathBuf,
    inner: csv::Reader<CommentNormalizingReader<R>>,
    line_infos: Rc<RefCell<VecDeque<LineInfo>>>,
}

impl DelimitedReader<File> {
    pub(crate) fn open(path: &Path, delimiter: Delimiter) -> Result<Self, DelimitedError> {
        let file = File::open(path).map_err(|source| DelimitedError::Open {
            path: path.to_path_buf(),
            source,
        })?;
        Ok(Self::from_reader(path, delimiter, file))
    }
}

impl<R: Read> DelimitedReader<R> {
    pub(crate) fn from_reader(path: &Path, delimiter: Delimiter, reader: R) -> Self {
        let line_infos = Rc::new(RefCell::new(VecDeque::new()));
        let reader =
            CommentNormalizingReader::new(reader, delimiter.byte(), Rc::clone(&line_infos));
        let inner = ReaderBuilder::new()
            .delimiter(delimiter.byte())
            .has_headers(false)
            .flexible(true)
            .comment(Some(b'#'))
            .from_reader(reader);
        Self {
            path: path.to_path_buf(),
            inner,
            line_infos,
        }
    }

    pub(crate) fn records(&mut self) -> DelimitedRecords<'_, R> {
        DelimitedRecords {
            path: &self.path,
            inner: self.inner.records(),
            line_infos: &self.line_infos,
        }
    }

    pub(crate) fn path(&self) -> &Path {
        &self.path
    }
}

struct CommentNormalizingReader<R: Read> {
    inner: BufReader<R>,
    line: Vec<u8>,
    offset: usize,
    delimiter: u8,
    in_quotes: bool,
    at_field_start: bool,
    byte_offset: u64,
    physical_line: u64,
    line_infos: Rc<RefCell<VecDeque<LineInfo>>>,
}

#[derive(Clone, Copy, Debug)]
struct LineInfo {
    start: u64,
    physical_line: u64,
}

impl<R: Read> CommentNormalizingReader<R> {
    fn new(inner: R, delimiter: u8, line_infos: Rc<RefCell<VecDeque<LineInfo>>>) -> Self {
        Self {
            inner: BufReader::new(inner),
            line: Vec::new(),
            offset: 0,
            delimiter,
            in_quotes: false,
            at_field_start: true,
            byte_offset: 0,
            physical_line: 1,
            line_infos,
        }
    }

    fn fill_line(&mut self) -> std::io::Result<bool> {
        self.line.clear();
        self.offset = 0;
        if self.inner.read_until(b'\n', &mut self.line)? == 0 {
            return Ok(false);
        }

        let first_non_whitespace = self
            .line
            .iter()
            .position(|byte| !byte.is_ascii_whitespace());
        let is_comment =
            !self.in_quotes && first_non_whitespace.is_some_and(|index| self.line[index] == b'#');
        if is_comment {
            if let Some(index) = first_non_whitespace.filter(|index| *index > 0) {
                self.line[0] = b'#';
                self.line[index] = b' ';
            }
            self.at_field_start = true;
            self.finish_physical_line();
            return Ok(true);
        }

        if !self.in_quotes && first_non_whitespace.is_some() {
            self.line_infos.borrow_mut().push_back(LineInfo {
                start: self.byte_offset,
                physical_line: self.physical_line,
            });
        }

        let mut index = 0;
        while index < self.line.len() {
            let byte = self.line[index];
            if self.in_quotes {
                if byte == b'"' {
                    if self.line.get(index + 1) == Some(&b'"') {
                        index += 2;
                        continue;
                    }
                    self.in_quotes = false;
                }
            } else if byte == b'"' && self.at_field_start {
                self.in_quotes = true;
                self.at_field_start = false;
            } else {
                self.at_field_start = byte == self.delimiter || byte == b'\n';
            }
            index += 1;
        }
        self.finish_physical_line();
        Ok(true)
    }

    fn finish_physical_line(&mut self) {
        self.byte_offset += self.line.len() as u64;
        self.physical_line += 1;
    }
}

impl<R: Read> Read for CommentNormalizingReader<R> {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        if buffer.is_empty() {
            return Ok(0);
        }
        if self.offset == self.line.len() && !self.fill_line()? {
            return Ok(0);
        }
        let available = &self.line[self.offset..];
        let copied = available.len().min(buffer.len());
        buffer[..copied].copy_from_slice(&available[..copied]);
        self.offset += copied;
        Ok(copied)
    }
}

pub(crate) struct DelimitedRecords<'a, R: Read> {
    path: &'a Path,
    inner: csv::StringRecordsIter<'a, CommentNormalizingReader<R>>,
    line_infos: &'a RefCell<VecDeque<LineInfo>>,
}

fn take_physical_line(line_infos: &RefCell<VecDeque<LineInfo>>, position: &csv::Position) -> u64 {
    let mut line_infos = line_infos.borrow_mut();
    let matching_index = line_infos
        .iter()
        .rposition(|line| line.start <= position.byte())
        .or_else(|| {
            line_infos
                .iter()
                .position(|line| line.start >= position.byte())
        });
    let Some(matching_index) = matching_index else {
        return position.line();
    };
    let physical_line = line_infos[matching_index].physical_line;
    line_infos.drain(..=matching_index);
    physical_line
}

impl<R: Read> Iterator for DelimitedRecords<'_, R> {
    type Item = Result<DelimitedRecord, DelimitedError>;

    fn next(&mut self) -> Option<Self::Item> {
        let result = self.inner.next()?;
        match result {
            Ok(fields) => {
                let position = fields.position();
                Some(Ok(DelimitedRecord {
                    line_number: position
                        .map_or(0, |position| take_physical_line(self.line_infos, position)),
                    fields,
                }))
            }
            Err(source) => {
                let position = source.position();
                Some(Err(DelimitedError::Read {
                    path: self.path.to_path_buf(),
                    record: position.map(|position| position.record() + 1),
                    line: position.map(|position| take_physical_line(self.line_infos, position)),
                    source,
                }))
            }
        }
    }
}

pub(crate) struct DelimitedWriter<W: Write> {
    path: PathBuf,
    inner: csv::Writer<W>,
    next_record: u64,
}

impl DelimitedWriter<BufWriter<File>> {
    pub(crate) fn create(
        path: &Path,
        delimiter: Delimiter,
        metadata_lines: &[String],
    ) -> Result<Self, DelimitedError> {
        Self::create_with_flexible_rows(path, delimiter, metadata_lines, false)
    }

    pub(crate) fn create_flexible(
        path: &Path,
        delimiter: Delimiter,
        metadata_lines: &[String],
    ) -> Result<Self, DelimitedError> {
        Self::create_with_flexible_rows(path, delimiter, metadata_lines, true)
    }

    fn create_with_flexible_rows(
        path: &Path,
        delimiter: Delimiter,
        metadata_lines: &[String],
        flexible: bool,
    ) -> Result<Self, DelimitedError> {
        let file = File::create(path).map_err(|source| DelimitedError::Create {
            path: path.to_path_buf(),
            source,
        })?;
        let mut writer = BufWriter::new(file);
        for line in metadata_lines {
            debug_assert!(line.starts_with('#'));
            writer
                .write_all(line.as_bytes())
                .and_then(|()| writer.write_all(b"\n"))
                .map_err(|source| DelimitedError::MetadataWrite {
                    path: path.to_path_buf(),
                    source,
                })?;
        }
        Ok(if flexible {
            Self::from_writer_with_flexible_rows(path, delimiter, writer, true)
        } else {
            Self::from_writer(path, delimiter, writer)
        })
    }
}

impl<W: Write> DelimitedWriter<W> {
    pub(crate) fn from_writer(path: &Path, delimiter: Delimiter, writer: W) -> Self {
        Self::from_writer_with_flexible_rows(path, delimiter, writer, false)
    }

    fn from_writer_with_flexible_rows(
        path: &Path,
        delimiter: Delimiter,
        writer: W,
        flexible: bool,
    ) -> Self {
        let inner = WriterBuilder::new()
            .delimiter(delimiter.byte())
            .has_headers(false)
            .flexible(flexible)
            .terminator(Terminator::Any(b'\n'))
            .comment(Some(b'#'))
            .from_writer(writer);
        Self {
            path: path.to_path_buf(),
            inner,
            next_record: 1,
        }
    }

    pub(crate) fn write_record<I, T>(&mut self, record: I) -> Result<(), DelimitedError>
    where
        I: IntoIterator<Item = T>,
        T: AsRef<[u8]>,
    {
        let record_number = self.next_record;
        self.inner
            .write_record(record)
            .map_err(|source| DelimitedError::Write {
                path: self.path.clone(),
                record: record_number,
                source,
            })?;
        self.next_record += 1;
        Ok(())
    }

    pub(crate) fn flush(&mut self) -> Result<(), DelimitedError> {
        self.inner.flush().map_err(|source| DelimitedError::Flush {
            path: self.path.clone(),
            source,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csv_and_tsv_round_trip_adversarial_fields_and_trailing_empty() {
        for delimiter in [Delimiter::Comma, Delimiter::Tab] {
            let path = Path::new("memory.table");
            let mut bytes = Vec::new();
            {
                let mut writer = DelimitedWriter::from_writer(path, delimiter, &mut bytes);
                writer
                    .write_record([
                        "#leading-comment-byte",
                        "plain",
                        "comma,value",
                        "tab\tvalue",
                        "quote\"value",
                        "line1\r\nline2",
                        "",
                    ])
                    .unwrap();
                writer.flush().unwrap();
            }

            let mut reader = DelimitedReader::from_reader(path, delimiter, bytes.as_slice());
            let records = reader.records().collect::<Result<Vec<_>, _>>().unwrap();
            assert_eq!(records.len(), 1);
            assert_eq!(
                records[0].fields().iter().collect::<Vec<_>>(),
                vec![
                    "#leading-comment-byte",
                    "plain",
                    "comma,value",
                    "tab\tvalue",
                    "quote\"value",
                    "line1\r\nline2",
                    "",
                ]
            );
        }
    }

    #[test]
    fn comments_are_ignored_without_losing_physical_line_context() {
        let input = b"#schema=v1\n  # indented comment\n\"#data\"\tvalue\n\"multi\n  # quoted data\"\tlast\n";
        let path = Path::new("comments.tsv");
        let mut reader = DelimitedReader::from_reader(path, Delimiter::Tab, &input[..]);
        let records = reader.records().collect::<Result<Vec<_>, _>>().unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].line_number(), 3);
        assert_eq!(records[0].fields().get(0), Some("#data"));
        assert_eq!(records[1].line_number(), 4);
        assert_eq!(records[1].fields().get(0), Some("multi\n  # quoted data"));
    }

    #[test]
    fn indented_comments_are_filtered_before_utf8_decoding() {
        let input = b"  # invalid comment bytes: \xff\xfe\nok\tvalue\n";
        let path = Path::new("comment-utf8.tsv");
        let mut reader = DelimitedReader::from_reader(path, Delimiter::Tab, &input[..]);
        let records = reader.records().collect::<Result<Vec<_>, _>>().unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].line_number(), 2);
        assert_eq!(records[0].fields().get(0), Some("ok"));
    }

    #[test]
    fn comment_tracking_memory_is_bounded_by_one_physical_line() {
        let input = (0..20_000)
            .map(|index| format!("  # comment {index}\n"))
            .chain(std::iter::once("value\tlast\n".to_owned()))
            .collect::<String>();
        let path = Path::new("many-comments.tsv");
        let mut reader = DelimitedReader::from_reader(path, Delimiter::Tab, input.as_bytes());
        let records = reader.records().collect::<Result<Vec<_>, _>>().unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].line_number(), 20_001);
        assert!(reader.inner.get_ref().line.capacity() < 1_024);
        assert!(reader.line_infos.borrow().is_empty());
    }

    #[test]
    fn data_record_tracking_is_pruned_as_records_are_consumed() {
        let input = (0..20_000)
            .map(|index| format!("value-{index}\tlast\n"))
            .collect::<String>();
        let path = Path::new("many-records.tsv");
        let mut reader = DelimitedReader::from_reader(path, Delimiter::Tab, input.as_bytes());
        let line_infos = Rc::clone(&reader.line_infos);
        let mut records = reader.records();
        for expected_line in 1..=20_000 {
            let record = records.next().unwrap().unwrap();
            assert_eq!(record.line_number(), expected_line);
            assert!(line_infos.borrow().len() <= 1);
        }
        assert!(records.next().is_none());
        assert!(line_infos.borrow().is_empty());
    }

    #[test]
    fn byte_positions_inside_a_bad_record_map_to_that_record_start() {
        let line_infos = RefCell::new(VecDeque::from([
            LineInfo {
                start: 10,
                physical_line: 3,
            },
            LineInfo {
                start: 30,
                physical_line: 4,
            },
        ]));
        let mut position = csv::Position::new();
        position.set_byte(20).set_line(99);
        assert_eq!(take_physical_line(&line_infos, &position), 3);
        assert_eq!(line_infos.borrow().front().unwrap().physical_line, 4);
    }

    #[test]
    fn malformed_records_report_path_record_and_physical_line() {
        let input = b"#schema=v1\nok\tvalue\nbad\t\xff\n";
        let path = Path::new("malformed.tsv");
        let mut reader = DelimitedReader::from_reader(path, Delimiter::Tab, &input[..]);
        let mut records = reader.records();

        assert_eq!(records.next().unwrap().unwrap().line_number(), 2);
        let error = records.next().unwrap().unwrap_err();
        match &error {
            DelimitedError::Read {
                path: error_path,
                record: Some(record),
                line: Some(line),
                ..
            } => {
                assert_eq!(error_path, path);
                assert_eq!(*record, 2);
                assert_eq!(*line, 3);
            }
            other => panic!("unexpected error: {other:?}"),
        }
        assert!(error.to_string().contains("malformed.tsv"));
    }
}
