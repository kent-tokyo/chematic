//! Explicit AUTO/CDX/CDXML molecular projection (RDKit issue #9604).
//! CDX is a bounded concrete-molecule profile; no query/reaction expansion.
mod binary;
mod properties;
mod xml;

use crate::{CdxmlError, CdxmlParseLimits, parse_cdxml_all_with_limits};
use chematic_core::Molecule;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChemDrawFormat {
    Auto,
    Cdx,
    Cdxml,
}
impl ChemDrawFormat {
    pub fn from_name(value: &str) -> Option<Self> {
        [Self::Auto, Self::Cdx, Self::Cdxml]
            .into_iter()
            .find(|f| value.eq_ignore_ascii_case(f.name()))
    }
    pub fn name(self) -> &'static str {
        match self {
            Self::Auto => "AUTO",
            Self::Cdx => "CDX",
            Self::Cdxml => "CDXML",
        }
    }
}
#[derive(Debug, Clone, Copy, Default)]
pub struct ChemDrawOptions {
    /// CDX only: omit known presentation properties/objects, reporting each tag.
    /// Unsupported chemistry remains an error, irrespective of this flag.
    pub allow_metadata_loss: bool,
}
#[derive(Debug, Clone, Copy)]
pub struct ChemDrawLimits {
    pub max_input_bytes: usize,
    pub max_items: usize,
    pub max_depth: usize,
    pub max_property_bytes: usize,
    pub max_output_bytes: usize,
    pub cdxml: CdxmlParseLimits,
}
impl Default for ChemDrawLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 16 * 1024 * 1024,
            max_items: 200_000,
            max_depth: 32,
            max_property_bytes: 1024 * 1024,
            max_output_bytes: 32 * 1024 * 1024,
            cdxml: CdxmlParseLimits::default(),
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ChemDrawError {
    Invalid {
        offset: usize,
        detail: &'static str,
    },
    Unsupported {
        offset: usize,
        tag: u16,
    },
    MetadataLoss {
        offset: usize,
        tag: u16,
    },
    Limit {
        resource: &'static str,
        limit: usize,
    },
    Cdxml(CdxmlError),
}
impl core::fmt::Display for ChemDrawError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid { offset, detail } => {
                write!(f, "invalid ChemDraw block at byte {offset}: {detail}")
            }
            Self::Unsupported { offset, tag } => write!(
                f,
                "unsupported CDX chemical/object tag 0x{tag:04x} at byte {offset}"
            ),
            Self::MetadataLoss { offset, tag } => write!(
                f,
                "CDX metadata loss requires opt-in: tag 0x{tag:04x} at byte {offset}"
            ),
            Self::Limit { resource, limit } => {
                write!(f, "ChemDraw {resource} exceeds limit {limit}")
            }
            Self::Cdxml(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for ChemDrawError {}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChemDrawOmission {
    pub offset: usize,
    pub tag: u16,
}
pub type ChemDrawMolecule = (Molecule, Vec<(f64, f64)>);
pub struct ChemDrawBlock {
    pub format: ChemDrawFormat,
    /// Exact source bytes; projected XML is not a source-preserving replacement.
    pub source: Vec<u8>,
    /// Molecular XML projection for CDX; exact UTF-8 source for CDXML.
    pub projected_cdxml: String,
    pub molecules: Vec<ChemDrawMolecule>,
    /// CDX omissions only. CDXML retains the existing graph-only reader policy.
    pub omissions: Vec<ChemDrawOmission>,
}
pub fn parse_chemdraw_block(
    input: &[u8],
    format: ChemDrawFormat,
    options: &ChemDrawOptions,
) -> Result<ChemDrawBlock, ChemDrawError> {
    parse_chemdraw_block_with_limits(input, format, options, &ChemDrawLimits::default())
}
pub fn parse_chemdraw_block_with_limits(
    input: &[u8],
    format: ChemDrawFormat,
    options: &ChemDrawOptions,
    limits: &ChemDrawLimits,
) -> Result<ChemDrawBlock, ChemDrawError> {
    if input.len() > limits.max_input_bytes {
        return Err(ChemDrawError::Limit {
            resource: "input bytes",
            limit: limits.max_input_bytes,
        });
    }
    let format = if format == ChemDrawFormat::Auto {
        if input.starts_with(b"VjCD0100") {
            ChemDrawFormat::Cdx
        } else {
            ChemDrawFormat::Cdxml
        }
    } else {
        format
    };
    let (xml, omissions) = if format == ChemDrawFormat::Cdx {
        binary::decode(input, options, limits)?
    } else {
        (
            std::str::from_utf8(input)
                .map_err(|_| ChemDrawError::Invalid {
                    offset: 0,
                    detail: "CDXML requires UTF-8 bytes",
                })?
                .to_owned(),
            Vec::new(),
        )
    };
    if xml.len() > limits.max_output_bytes {
        return Err(ChemDrawError::Limit {
            resource: "projected XML bytes",
            limit: limits.max_output_bytes,
        });
    }
    xml::validate(&xml, limits)?;
    let mut molecules =
        parse_cdxml_all_with_limits(&xml, &limits.cdxml).map_err(ChemDrawError::Cdxml)?;
    if format == ChemDrawFormat::Cdx {
        for (m, _) in &mut molecules {
            let aromatic: Vec<_> = m
                .bonds()
                .filter(|(_, b)| b.order == chematic_core::BondOrder::Aromatic)
                .flat_map(|(_, b)| [b.atom1, b.atom2])
                .collect();
            for a in aromatic {
                m.set_atom_aromatic(a, true);
            }
        }
    }
    // The legacy CDXML reader computes an E/Z CIP cache but does not populate
    // the direction metadata consumed by SMILES/MOL writers. Complete that
    // geometry-derived representation only in this new projection API.
    for (m, coords) in &mut molecules {
        chematic_perception::apply_ez_directions_from_2d(m, coords);
    }
    Ok(ChemDrawBlock {
        format,
        source: input.to_vec(),
        projected_cdxml: xml,
        molecules,
        omissions,
    })
}
