//! Error types for SMILES parsing and writing.

use core::fmt;

/// Errors produced during SMILES parsing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SmilesError {
    /// The input or parsed graph exceeded a configured resource limit.
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
    /// Input string ended unexpectedly.
    UnexpectedEnd { pos: usize },
    /// Unrecognised element symbol inside a bracket atom.
    UnknownElement { symbol: String, pos: usize },
    /// A ring-closure digit was opened but never closed (or vice-versa).
    UnmatchedRingClosure { ring_num: u8, pos: usize },
    /// An extended SMILES+ `%(n)` ring closure above the legacy `u8` range
    /// was opened but never closed (or vice-versa).
    UnmatchedExtendedRingClosure { ring_num: u32, pos: usize },
    /// Mismatched parentheses.
    MismatchedParentheses { pos: usize },
    /// A bracket atom `[...]` could not be parsed.
    InvalidBracketAtom { detail: String, pos: usize },
    /// Conflicting bond types at both ends of a ring closure.
    ConflictingRingBond { ring_num: u8, pos: usize },
    /// An extended SMILES+ `%(n)` ring closure above the legacy `u8` range
    /// declared conflicting bond types at its two ends.
    ConflictingExtendedRingBond { ring_num: u32, pos: usize },
    /// Empty SMILES string.
    EmptyInput,
    /// Branch nesting exceeded the safe recursion limit.
    NestingTooDeep { pos: usize },
    /// Trailing input after a structurally complete SMILES chain could not be
    /// parsed (e.g. an unrecognised atom symbol or stray punctuation).
    UnexpectedCharacter { pos: usize },
    /// An OpenSMILES extended chirality class with a permutation number
    /// outside its range (`@TH3`, `@AL0`, `@SP4`, `@TB21`, `@OH31`, ...), which
    /// RDKit also rejects. `class` is the tag text without the leading `@`
    /// (e.g. `"TB21"`).
    UnsupportedChiralityClass { class: String, pos: usize },
    /// A neutral oxygen or fluorine has more bonds (bond orders plus explicit
    /// hydrogens) than its element permits, e.g. the three-coordinate oxygen
    /// of `O=O1C=CC=C1` or `CO(C)C` (#769). Charged atoms (`[O+]`, `[o+]`)
    /// and the hypervalent states of other elements are not checked here.
    /// `atom` is the zero-based atom index in input order.
    InvalidValence {
        element: &'static str,
        atom: usize,
        valence: u32,
        max: u32,
    },
}

impl fmt::Display for SmilesError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(f, "SMILES {resource} exceeds limit {limit} (got {actual})"),
            Self::UnexpectedEnd { pos } => write!(f, "unexpected end of input at position {pos}"),
            Self::UnknownElement { symbol, pos } => {
                write!(f, "unknown element '{symbol}' at position {pos}")
            }
            Self::UnmatchedRingClosure { ring_num, pos } => {
                write!(f, "unmatched ring closure {ring_num} at position {pos}")
            }
            Self::UnmatchedExtendedRingClosure { ring_num, pos } => {
                write!(
                    f,
                    "unmatched extended ring closure {ring_num} at position {pos}"
                )
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
            Self::ConflictingExtendedRingBond { ring_num, pos } => write!(
                f,
                "conflicting bond types for extended ring closure {ring_num} at position {pos}"
            ),
            Self::EmptyInput => write!(f, "SMILES input is empty"),
            Self::NestingTooDeep { pos } => write!(f, "branch nesting too deep at position {pos}"),
            Self::UnexpectedCharacter { pos } => {
                write!(f, "unexpected character at position {pos}")
            }
            Self::UnsupportedChiralityClass { class, pos } => write!(
                f,
                "unsupported chirality class '@{class}' at position {pos}: valid are @TH1-2, @AL1-2, @SP1-3, @TB1-20 and @OH1-30 (or the class without a number)"
            ),
            Self::InvalidValence {
                element,
                atom,
                valence,
                max,
            } => write!(
                f,
                "explicit valence {valence} for neutral {element} atom #{atom} is greater than permitted ({max})"
            ),
        }
    }
}

impl std::error::Error for SmilesError {}
