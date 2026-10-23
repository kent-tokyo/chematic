//! Bounded structured ORCA output, retaining finite Cartesian f64 values.
use super::{OrcaCoords, OrcaInput, OrcaInputParseLimits, parse_orca_input_with_limits};
use std::fmt::Write;

/// Independent transport-input, emitted-text and atom budgets.
#[derive(Debug, Clone, Copy)]
pub struct OrcaWriteLimits {
    pub max_input_bytes: usize,
    pub max_output_bytes: usize,
    pub max_atoms: usize,
}
impl Default for OrcaWriteLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 * 1024 * 1024,
            max_output_bytes: 64 * 1024 * 1024,
            max_atoms: 100_000,
        }
    }
}

/// No partial output is returned on validation, preservation or budget failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OrcaWriteError(pub String);
impl core::fmt::Display for OrcaWriteError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "ORCA checked output: {}", self.0)
    }
}
impl std::error::Error for OrcaWriteError {}
pub(super) fn error(detail: impl ToString) -> OrcaWriteError {
    OrcaWriteError(detail.to_string())
}

struct Output {
    text: String,
    cap: usize,
}
impl Write for Output {
    fn write_str(&mut self, text: &str) -> std::fmt::Result {
        if self
            .text
            .len()
            .checked_add(text.len())
            .is_none_or(|n| n > self.cap)
        {
            return Err(std::fmt::Error);
        }
        self.text
            .try_reserve(text.len())
            .map_err(|_| std::fmt::Error)?;
        self.text.push_str(text);
        Ok(())
    }
}

fn render(input: &OrcaInput, out: &mut Output) -> std::fmt::Result {
    for comment in &input.comments {
        writeln!(out, "{comment}")?;
    }
    if !input.keywords.is_empty() {
        out.write_str("! ")?;
        for (index, keyword) in input.keywords.iter().enumerate() {
            if index != 0 {
                out.write_char(' ')?;
            }
            out.write_str(keyword)?;
        }
        out.write_char('\n')?;
    }
    for block in &input.blocks {
        write!(out, "%{}", block.name)?;
        if block.has_end {
            write!(out, "\n{}\nend\n", block.raw)?;
        } else {
            writeln!(out, " {}", block.raw)?;
        }
    }
    let Some(OrcaCoords::Xyz {
        charge,
        multiplicity,
        atoms,
    }) = &input.coords
    else {
        return Err(std::fmt::Error);
    };
    writeln!(out, "* xyz {charge} {multiplicity}")?;
    for atom in atoms {
        out.write_str(atom.element.symbol())?;
        for (value, frozen) in [atom.x, atom.y, atom.z].into_iter().zip(atom.frozen) {
            // 17 significant decimal digits retain every finite binary64 value.
            write!(out, " {value:.16e}")?;
            if frozen {
                out.write_char('$')?;
            }
        }
        if let Some(extra) = &atom.extra {
            write!(out, " {extra}")?;
        }
        out.write_char('\n')?;
    }
    out.write_str("*\n")
}

fn retained(before: &OrcaInput, after: &OrcaInput) -> bool {
    if before.comments != after.comments
        || before.keywords != after.keywords
        || before.blocks != after.blocks
    {
        return false;
    }
    let (
        Some(OrcaCoords::Xyz {
            charge: ac,
            multiplicity: am,
            atoms: a,
        }),
        Some(OrcaCoords::Xyz {
            charge: bc,
            multiplicity: bm,
            atoms: b,
        }),
    ) = (&before.coords, &after.coords)
    else {
        return false;
    };
    ac == bc
        && am == bm
        && a.len() == b.len()
        && a.iter().zip(b).all(|(a, b)| {
            a.element == b.element
                && a.frozen == b.frozen
                && a.extra == b.extra
                && [a.x, a.y, a.z].map(f64::to_bits) == [b.x, b.y, b.z].map(f64::to_bits)
        })
}

/// Write one embedded xyz input with no coordinate rounding or unit conversion.
pub fn write_orca_input_checked(input: &OrcaInput) -> Result<String, OrcaWriteError> {
    write_orca_input_checked_with_limits(input, &OrcaWriteLimits::default())
}

/// Preserve all represented comments, keywords, opaque blocks and ordered xyz
/// atom fields on reparse, including coordinate f64 bits and freeze flags.
/// Only the Cartesian Angs/Bohrs/default-unit profile is supported. Conflicting
/// units, %coords, nonfinite/empty geometry, spin zero, unrepresentable fields
/// and budgets are refused. Unknown job settings remain opaque, not validated.
/// This does not preserve source layout or restore fields omitted during parse.
pub fn write_orca_input_checked_with_limits(
    input: &OrcaInput,
    limits: &OrcaWriteLimits,
) -> Result<String, OrcaWriteError> {
    super::cartesian::coordinate_units(&input.keywords).map_err(error)?;
    if input
        .blocks
        .iter()
        .any(|b| b.name.eq_ignore_ascii_case("coords"))
    {
        return Err(error("%coords is outside the checked Cartesian profile"));
    }
    let Some(OrcaCoords::Xyz {
        atoms,
        multiplicity,
        ..
    }) = &input.coords
    else {
        return Err(error("one embedded xyz geometry is required"));
    };
    if atoms.is_empty() || atoms.len() > limits.max_atoms || *multiplicity == 0 {
        return Err(error(
            "nonempty bounded geometry and positive multiplicity required",
        ));
    }
    if atoms
        .iter()
        .any(|a| ![a.x, a.y, a.z].into_iter().all(f64::is_finite))
    {
        return Err(error("coordinates must be finite"));
    }
    let mut out = Output {
        text: String::new(),
        cap: limits.max_output_bytes,
    };
    render(input, &mut out).map_err(|_| error("output byte limit or allocation failed"))?;
    let parsed = parse_orca_input_with_limits(
        &out.text,
        &OrcaInputParseLimits {
            max_input_bytes: limits.max_output_bytes,
            max_atoms: limits.max_atoms,
            ..Default::default()
        },
    )
    .map_err(error)?;
    if !retained(input, &parsed) {
        return Err(error("represented fields would change on reparse"));
    }
    Ok(out.text)
}
