use thiserror::Error;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Strand {
    Plus,
    Minus,
    Unknown,
}

#[derive(Error, Debug)]
pub enum StrandParseError {
    #[error("invalid strand {value:?}")]
    Invalid { value: String },
}

impl Strand {
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
