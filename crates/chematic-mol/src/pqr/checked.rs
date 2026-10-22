//! Checked atom-record output, independently based on the APBS PQR field layout.
use super::{PqrAtomRecord, infer_element};

/// Output budgets; atom count and emitted UTF-8 bytes are bounded independently.
#[derive(Debug, Clone, Copy)]
pub struct PqrWriteLimits {
    pub max_atoms: usize,
    pub max_output_bytes: usize,
}
impl Default for PqrWriteLimits {
    fn default() -> Self {
        Self {
            max_atoms: 100_000,
            max_output_bytes: 64 * 1024 * 1024,
        }
    }
}

/// Checked output failure. `atom` is zero based; None denotes a document limit.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PqrWriteError {
    pub atom: Option<usize>,
    pub detail: String,
}
impl core::fmt::Display for PqrWriteError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "PQR output {:?}: {}", self.atom, self.detail)
    }
}
impl std::error::Error for PqrWriteError {}

fn error(atom: Option<usize>, detail: &str) -> PqrWriteError {
    PqrWriteError {
        atom,
        detail: detail.into(),
    }
}
fn token(value: &str) -> bool {
    !value.is_empty() && value.len() <= 64 && value.bytes().all(|b| b.is_ascii_graphic())
}
fn validate(atom: &PqrAtomRecord, index: usize) -> Result<(), PqrWriteError> {
    if !matches!(atom.group_pdb.as_str(), "ATOM" | "HETATM") {
        return Err(error(Some(index), "record must be ATOM or HETATM"));
    }
    if !token(&atom.atom_name)
        || !token(&atom.res_name)
        || atom.chain_id.as_deref().is_some_and(|s| !token(s))
        || atom.icode.is_some_and(|c| !c.is_ascii_alphabetic())
    {
        return Err(error(
            Some(index),
            "names/chain must be single ASCII tokens of 1..=64 bytes; insertion code must be one ASCII letter",
        ));
    }
    if infer_element(&atom.group_pdb, &atom.res_name, &atom.atom_name) != Some(atom.element) {
        return Err(error(
            Some(index),
            "element cannot be retained by the PQR name-inference rule",
        ));
    }
    if ![atom.x, atom.y, atom.z, atom.charge, atom.radius]
        .iter()
        .all(|v| v.is_finite())
        || atom.radius < 0.0
    {
        return Err(error(
            Some(index),
            "coordinates/charge/radius must be finite; radius must be nonnegative",
        ));
    }
    Ok(())
}
fn scalar(value: f64, index: usize) -> Result<String, PqrWriteError> {
    let text = format!("{value:.16e}");
    if text.parse::<f64>().map(f64::to_bits).ok() != Some(value.to_bits()) {
        return Err(error(
            Some(index),
            "numeric field cannot retain its f64 value",
        ));
    }
    Ok(text)
}

/// Write supplied PQR atom fields without rounding, using default output budgets.
pub fn write_pqr_checked(atoms: &[PqrAtomRecord]) -> Result<String, PqrWriteError> {
    write_pqr_checked_with_limits(atoms, &PqrWriteLimits::default())
}

/// Retain ordered atom metadata and all five finite f64 fields bit for bit,
/// including signed zero/subnormals. PQR uses Angstrom coordinates/radii and
/// electron-unit partial charges; no values are computed or rescaled. No bonds,
/// charge model, duplicate-serial policy or chemistry correctness are inferred.
///
/// ASCII single-token names/chain and an optional alphabetic insertion code
/// survive reparse. The supplied element must match PQR's documented inference
/// heuristic, since the format has no element column. Empty input, incompatible
/// fields, negative radii and limits are errors without partial output. This is
/// atom-record output; original non-atom records/formatting are not represented.
pub fn write_pqr_checked_with_limits(
    atoms: &[PqrAtomRecord],
    limits: &PqrWriteLimits,
) -> Result<String, PqrWriteError> {
    if atoms.is_empty() || atoms.len() > limits.max_atoms {
        return Err(error(
            None,
            "atom count must be nonzero and within the output limit",
        ));
    }
    let mut output = String::new();
    for (index, atom) in atoms.iter().enumerate() {
        validate(atom, index)?;
        let mut line = format!(
            "{} {} {} {} ",
            atom.group_pdb, atom.serial, atom.atom_name, atom.res_name
        );
        if let Some(chain) = &atom.chain_id {
            line.push_str(chain);
            line.push(' ');
        }
        line.push_str(&atom.res_seq.to_string());
        if let Some(code) = atom.icode {
            line.push(code);
        }
        for value in [atom.x, atom.y, atom.z, atom.charge, atom.radius] {
            line.push(' ');
            line.push_str(&scalar(value, index)?);
        }
        line.push('\n');
        if output
            .len()
            .checked_add(line.len())
            .is_none_or(|size| size > limits.max_output_bytes)
        {
            return Err(error(None, "output byte limit exceeded"));
        }
        output
            .try_reserve(line.len())
            .map_err(|_| error(None, "output allocation failed"))?;
        output.push_str(&line);
    }
    Ok(output)
}
