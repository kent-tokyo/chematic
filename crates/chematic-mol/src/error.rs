//! Error types for MOL/SDF parsing.

/// Errors that can occur while parsing a MOL V2000 or SDF file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MolParseError {
    /// The header block (lines 1–3) could not be parsed.
    InvalidHeader { line: usize, detail: String },
    /// The counts line (line 4) is missing or malformed.
    InvalidCountLine { line: usize, detail: String },
    /// An atom-block line could not be parsed.
    InvalidAtomLine { line: usize, detail: String },
    /// A bond-block line could not be parsed.
    InvalidBondLine { line: usize, detail: String },
    /// The element symbol on an atom line is not recognised.
    UnknownElement { symbol: String, line: usize },
    /// The input ended before the molecule was complete.
    UnexpectedEnd,
    /// A V3000 (Extended Ctab) structural error.
    V3000ParseError { line: usize, msg: String },
    /// An IO error occurred while reading the input stream.
    Io(String),
}

impl std::fmt::Display for MolParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidHeader { line, detail } => {
                write!(f, "invalid header at line {line}: {detail}")
            }
            Self::InvalidCountLine { line, detail } => {
                write!(f, "invalid counts line at line {line}: {detail}")
            }
            Self::InvalidAtomLine { line, detail } => {
                write!(f, "invalid atom line at line {line}: {detail}")
            }
            Self::InvalidBondLine { line, detail } => {
                write!(f, "invalid bond line at line {line}: {detail}")
            }
            Self::UnknownElement { symbol, line } => {
                write!(f, "unknown element symbol '{symbol}' at line {line}")
            }
            Self::UnexpectedEnd => {
                write!(f, "unexpected end of input")
            }
            Self::V3000ParseError { line, msg } => {
                write!(f, "V3000 parse error at line {line}: {msg}")
            }
            Self::Io(msg) => write!(f, "IO error: {msg}"),
        }
    }
}

impl std::error::Error for MolParseError {}
