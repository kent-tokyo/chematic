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
        }
    }
}

impl std::error::Error for SmilesError {}
