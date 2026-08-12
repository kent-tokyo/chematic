//! Error types for SMILES parsing and writing.

use core::fmt;

/// Errors produced during SMILES parsing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SmilesError {
    /// Input string ended unexpectedly.
    UnexpectedEnd { pos: usize },
    /// Unrecognised element symbol inside a bracket atom.
    UnknownElement { symbol: String, pos: usize },
    /// A ring-closure digit was opened but never closed (or vice-versa).
    UnmatchedRingClosure { ring_num: u8, pos: usize },
    /// Mismatched parentheses.
    MismatchedParentheses { pos: usize },
    /// A bracket atom `[...]` could not be parsed.
    InvalidBracketAtom { detail: String, pos: usize },
    /// Conflicting bond types at both ends of a ring closure.
    ConflictingRingBond { ring_num: u8, pos: usize },
    /// Empty SMILES string.
    EmptyInput,
    /// Branch nesting exceeded the safe recursion limit.
    NestingTooDeep { pos: usize },
    /// Trailing input after a structurally complete SMILES chain could not be
    /// parsed (e.g. an unrecognised atom symbol or stray punctuation).
    UnexpectedCharacter { pos: usize },
    /// A recognized OpenSMILES extended chirality class (`@TH`/`@AL`/`@TB`/`@OH`) or
    /// an out-of-range square-planar tag (e.g. `@SP4`) -- syntactically valid, but not
    /// implemented. `class` is the tag text without the leading `@` (e.g. `"TB4"`).
    /// Only `@SP1`/`@SP2`/`@SP3` are currently supported, alongside plain `@`/`@@`.
    UnsupportedChiralityClass { class: String, pos: usize },
}

impl fmt::Display for SmilesError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEnd { pos } => write!(f, "unexpected end of input at position {pos}"),
            Self::UnknownElement { symbol, pos } => {
                write!(f, "unknown element '{symbol}' at position {pos}")
            }
            Self::UnmatchedRingClosure { ring_num, pos } => {
                write!(f, "unmatched ring closure {ring_num} at position {pos}")
            }
            Self::MismatchedParentheses { pos } => {
                write!(f, "mismatched parenthesis at position {pos}")
            }
            Self::InvalidBracketAtom { detail, pos } => {
                write!(f, "invalid bracket atom at position {pos}: {detail}")
            }
            Self::ConflictingRingBond { ring_num, pos } => write!(
                f,
                "conflicting bond types for ring closure {ring_num} at position {pos}"
            ),
            Self::EmptyInput => write!(f, "SMILES input is empty"),
            Self::NestingTooDeep { pos } => write!(f, "branch nesting too deep at position {pos}"),
            Self::UnexpectedCharacter { pos } => {
                write!(f, "unexpected character at position {pos}")
            }
            Self::UnsupportedChiralityClass { class, pos } => write!(
                f,
                "unsupported chirality class '@{class}' at position {pos}: only tetrahedral (@/@@) and square-planar (@SP1/@SP2/@SP3) are currently supported"
            ),
        }
    }
}

impl std::error::Error for SmilesError {}
