//! Checked enhanced-atom-stereo reaction CXSMILES (RDKit #9732).
//!
//! Supports `a`, `&n` and `on` atom groups. Group numbers form one namespace
//! across reactants, agents and products; the same kind/number joins members
//! across participants. Other CX fields are refused. Legacy APIs are unchanged.

use crate::{Reaction, ReactionParseLimits, RxnError};
use chematic_core::{AtomIdx, Molecule, StereoGroup, StereoGroupKind};

mod parse;
mod write;

/// An invalid or unsupported enhanced-stereo reaction representation.
#[derive(Debug)]
pub enum ReactionCxError {
    Reaction(RxnError),
    Invalid(String),
    Unsupported(String),
}

impl std::fmt::Display for ReactionCxError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Reaction(error) => error.fmt(f),
            Self::Invalid(message) => write!(f, "reaction CXSMILES: {message}"),
            Self::Unsupported(message) => write!(f, "reaction CXSMILES unsupported: {message}"),
        }
    }
}
impl std::error::Error for ReactionCxError {}
impl From<RxnError> for ReactionCxError {
    fn from(error: RxnError) -> Self {
        Self::Reaction(error)
    }
}

/// Parse reaction SMILES with checked enhanced tetrahedral atom groups.
///
/// Global CX atom indices traverse reactants, agents, products and each
/// component in source order. Group members must have tetrahedral annotations.
/// Unknown CX fields, titles, overlaps and invalid references are errors.
pub fn parse_reaction_cxsmiles(input: &str) -> Result<Reaction, ReactionCxError> {
    parse_reaction_cxsmiles_with_limits(input, &ReactionParseLimits::default())
}

/// Parse with caller-selected reaction limits; the byte limit includes CX data.
pub fn parse_reaction_cxsmiles_with_limits(
    input: &str,
    limits: &ReactionParseLimits,
) -> Result<Reaction, ReactionCxError> {
    parse::parse(input, limits)
}

/// Write enhanced stereo after remapping every member through the actual
/// emitted atom order. `canonical` canonicalizes each molecule and sorts
/// components within each role, never across the two arrow boundaries.
///
/// This is not a globally unique key for enhanced-stereo reactions. Group
/// numbers are preserved, including groups shared across participants. Bond
/// groups and group members removed by molecular canonicalization are errors.
pub fn write_reaction_cxsmiles(
    reaction: &Reaction,
    canonical: bool,
) -> Result<String, ReactionCxError> {
    write::write(reaction, canonical)
}

type GroupKey = (u8, u32);
fn group_key(kind: &StereoGroupKind) -> Result<GroupKey, ReactionCxError> {
    match *kind {
        StereoGroupKind::Absolute => Ok((0, 0)),
        StereoGroupKind::And(n) if n > 0 => Ok((1, n)),
        StereoGroupKind::Or(n) if n > 0 => Ok((2, n)),
        _ => Err(ReactionCxError::Invalid(
            "AND/OR group number must be positive".into(),
        )),
    }
}
fn group_head((kind, number): GroupKey) -> String {
    match kind {
        0 => "a".into(),
        1 => format!("&{number}"),
        _ => format!("o{number}"),
    }
}
