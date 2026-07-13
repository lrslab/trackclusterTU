//! Genomic strand values and strict text parsing.

use thiserror::Error;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
/// Genomic alignment or feature strand.
pub enum Strand {
    /// Forward (`+`) strand.
    Plus,
    /// Reverse (`-`) strand.
    Minus,
    /// Unspecified (`.`) strand.
    Unknown,
}

#[derive(Error, Debug)]
/// Error returned when text is not `+`, `-`, or `.`.
pub enum StrandParseError {
    /// Invalid strand text.
    #[error("invalid strand {value:?}")]
    Invalid {
        /// Original invalid value.
        value: String,
    },
}

impl Strand {
    /// Return the BED/GFF character representing the strand.
    pub fn as_char(self) -> char {
        match self {
            Self::Plus => '+',
            Self::Minus => '-',
            Self::Unknown => '.',
        }
    }
}

impl TryFrom<&str> for Strand {
    type Error = StrandParseError;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "+" => Ok(Self::Plus),
            "-" => Ok(Self::Minus),
            "." => Ok(Self::Unknown),
            _ => Err(StrandParseError::Invalid {
                value: value.to_owned(),
            }),
        }
    }
}
