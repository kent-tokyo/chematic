//! Strict transport matching Python's kind and WASM's type coordinate tags.
use super::checked::error;
use super::{OrcaAtom, OrcaBlock, OrcaCoords, OrcaInput, OrcaWriteError, OrcaWriteLimits};

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct Document {
    comments: Vec<String>,
    keywords: Vec<String>,
    blocks: Vec<Block>,
    coords: Coordinates,
}
#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct Block {
    name: String,
    raw: String,
    has_end: bool,
}
#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct Coordinates {
    #[serde(alias = "type")]
    kind: String,
    charge: i32,
    multiplicity: u32,
    atoms: Vec<Atom>,
}
#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct Atom {
    element: String,
    x: f64,
    y: f64,
    z: f64,
    frozen: [bool; 3],
    extra: Option<String>,
}
impl Document {
    fn into_input(self, limits: &OrcaWriteLimits) -> Result<OrcaInput, OrcaWriteError> {
        if self.coords.kind != "xyz" || self.coords.atoms.len() > limits.max_atoms {
            return Err(error("bounded embedded xyz geometry required"));
        }
        let atoms = self
            .coords
            .atoms
            .into_iter()
            .map(|a| {
                let element = chematic_core::Element::from_symbol(&a.element)
                    .ok_or_else(|| error("unknown element"))?;
                Ok(OrcaAtom {
                    element,
                    x: a.x,
                    y: a.y,
                    z: a.z,
                    frozen: a.frozen,
                    extra: a.extra,
                })
            })
            .collect::<Result<Vec<_>, OrcaWriteError>>()?;
        Ok(OrcaInput {
            comments: self.comments,
            keywords: self.keywords,
            blocks: self
                .blocks
                .into_iter()
                .map(|b| OrcaBlock {
                    name: b.name,
                    raw: b.raw,
                    has_end: b.has_end,
                })
                .collect(),
            coords: Some(OrcaCoords::Xyz {
                charge: self.coords.charge,
                multiplicity: self.coords.multiplicity,
                atoms,
            }),
        })
    }
}

/// Strict JSON adapter for checked structured output. Required fields are never
/// defaulted; unknown/duplicate fields and bool numeric scalars are refused.
/// `kind` or `type` identifies the xyz variant, but both together are refused.
/// Optional atom `extra` may be omitted/null. Input bytes are checked before
/// JSON decoding; decoded allocation is bounded by bytes, not a process cap.
pub fn write_orca_input_json_checked_with_limits(
    input_json: &str,
    limits: &OrcaWriteLimits,
) -> Result<String, OrcaWriteError> {
    if input_json.len() > limits.max_input_bytes {
        return Err(error("JSON input byte limit exceeded"));
    }
    let input: Document = serde_json::from_str(input_json).map_err(error)?;
    super::write_orca_input_checked_with_limits(&input.into_input(limits)?, limits)
}

/// Checked JSON output with default transport/text/atom budgets.
pub fn write_orca_input_json_checked(input_json: &str) -> Result<String, OrcaWriteError> {
    write_orca_input_json_checked_with_limits(input_json, &OrcaWriteLimits::default())
}
