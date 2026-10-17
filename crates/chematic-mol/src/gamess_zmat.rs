//! GAMESS-US Cartesian geometry plus explicit simple internal coordinates.
//! Uses the public GAMESS input description; no GAMESS/Open Babel code is used.
use crate::gzmat::{ZMatrixError, ZMatrixRow, cartesian_to_zmatrix};
use chematic_core::Element;

pub const MAX_GAMESS_ZMAT_ATOMS: usize = 1024;
const MAX_OUTPUT_BYTES: usize = 256 * 1024;

fn error(atom: Option<usize>, detail: impl Into<String>) -> ZMatrixError {
    ZMatrixError {
        atom,
        detail: detail.into(),
    }
}

/// Geometry template options; no Hamiltonian, basis or SCF type is inferred.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GamessZmatWriteOptions {
    pub charge: i32,
    pub multiplicity: u32,
    /// Significant decimal digits for Cartesian coordinates, 1..=17.
    pub significant_digits: u8,
    pub title: String,
}
impl GamessZmatWriteOptions {
    pub fn new(charge: i32, multiplicity: u32) -> Self {
        Self {
            charge,
            multiplicity,
            significant_digits: 17,
            title: "chematic GAMESS internal coordinates".into(),
        }
    }
    fn validate(&self) -> Result<(), ZMatrixError> {
        if !(1..=17).contains(&self.significant_digits) || self.multiplicity == 0 {
            return Err(error(
                None,
                "digits must be 1..=17 and multiplicity must be positive",
            ));
        }
        if self.title.trim().is_empty()
            || self.title.len() > 80
            || !self.title.is_ascii()
            || self
                .title
                .chars()
                .any(|c| c.is_control() || matches!(c, '$' | '!'))
        {
            return Err(error(
                None,
                "title must be 1..80 printable ASCII bytes without $ or !",
            ));
        }
        Ok(())
    }
}

fn coordinate_lines(
    elements: &[Element],
    positions: &[[f64; 3]],
    digits: u8,
) -> Result<(Vec<String>, Vec<[f64; 3]>), ZMatrixError> {
    let mut lines = Vec::new();
    let mut rounded = Vec::new();
    for (i, (element, position)) in elements.iter().zip(positions).enumerate() {
        let values: Vec<_> = position
            .iter()
            .map(|v| {
                let value = if *v == 0.0 { 0.0 } else { *v };
                format!("{:.*E}", usize::from(digits) - 1, value)
            })
            .collect();
        let mut row = [0.; 3];
        for (slot, value) in row.iter_mut().zip(&values) {
            *slot = value
                .parse()
                .map_err(|_| error(Some(i), "invalid formatted coordinate"))?;
        }
        let line = format!(
            "{} {}.0 {} {} {}\n",
            element.symbol(),
            element.atomic_number(),
            values[0],
            values[1],
            values[2]
        );
        if line.len() - 1 > 80 {
            return Err(error(
                Some(i),
                "Cartesian card exceeds 80 columns; reduce significant_digits",
            ));
        }
        lines.push(line);
        rounded.push(row);
    }
    Ok((lines, rounded))
}

fn append_internals(output: &mut String, rows: &[ZMatrixRow], count: usize) {
    if count == 0 {
        return;
    }
    output.push_str(" $ZMAT IZMAT(1)=\n");
    let mut emitted = 0;
    for (i, row) in rows.iter().enumerate() {
        for kind in 1..=row.references.len() {
            output.push_str(&format!(" {kind},{}", i + 1));
            for reference in &row.references[..kind] {
                output.push_str(&format!(",{}", reference + 1));
            }
            emitted += 1;
            output.push_str(if emitted == count { "\n" } else { ",\n" });
        }
    }
    output.push_str(" $END\n");
}

/// Write a GAMESS-US geometry template with C1/COORD=UNIQUE Cartesian $DATA
/// and caller-defined $ZMAT IZMAT distances/angles/torsions. References are
/// zero-based distinct earlier atoms, as in `cartesian_to_zmatrix`.
///
/// NZVAR is 0 for one atom, 1 for two, and 3N-6 for non-linear larger inputs.
/// One-atom output omits $ZMAT. Linear/near-linear reference frames fail.
/// Absolute coordinates and atom order are retained in $DATA; no basis,
/// isotope, topology, electronic-state plausibility or job execution is supplied.
/// All cards fit 80 columns; 17 digits round-trip the original finite f64s.
pub fn write_gamess_zmat(
    elements: &[Element],
    positions_angstrom: &[[f64; 3]],
    references: &[Vec<usize>],
    options: &GamessZmatWriteOptions,
) -> Result<String, ZMatrixError> {
    options.validate()?;
    if elements.is_empty() || elements.len() > MAX_GAMESS_ZMAT_ATOMS {
        return Err(error(None, "GAMESS ZMAT atom count must be 1..=1024"));
    }
    let rows = cartesian_to_zmatrix(elements, positions_angstrom, references)?;
    let (lines, rounded) =
        coordinate_lines(elements, positions_angstrom, options.significant_digits)?;
    // Lower precision must not destroy a required internal-coordinate frame.
    cartesian_to_zmatrix(elements, &rounded, references)?;
    let count: usize = rows.iter().map(|row| row.references.len()).sum();
    let mut output = String::new();
    output
        .try_reserve(MAX_OUTPUT_BYTES)
        .map_err(|_| error(None, "GAMESS output allocation failed"))?;
    output.push_str(&format!(
        " $CONTRL COORD=UNIQUE UNITS=ANGS\n ICHARG={} MULT={} NZVAR={}\n $END\n $DATA\n{}\nC1\n",
        options.charge, options.multiplicity, count, options.title
    ));
    for line in lines {
        output.push_str(&line);
    }
    output.push_str(" $END\n");
    append_internals(&mut output, &rows, count);
    if output.len() > MAX_OUTPUT_BYTES {
        return Err(error(None, "GAMESS output byte limit exceeded"));
    }
    Ok(output)
}

/// Convert exactly one plain XYZ frame (all atoms retained) to the geometry
/// template. Input is bounded to 1 MiB, 1024 atoms, and 4096 bytes per line.
/// Multiple frames or trailing non-XYZ content are rejected.
pub fn xyz_to_gamess_zmat(
    xyz: &str,
    references: &[Vec<usize>],
    options: &GamessZmatWriteOptions,
) -> Result<String, ZMatrixError> {
    let limits = crate::XyzParseLimits {
        max_input_bytes: 1 << 20,
        max_atoms_per_frame: MAX_GAMESS_ZMAT_ATOMS,
        max_frames: 1,
        max_line_bytes: 4096,
    };
    let frames = crate::parse_xyz_all_with_limits(xyz, &limits)
        .map_err(|e| error(None, format!("XYZ input: {e}")))?;
    let frame = frames
        .first()
        .ok_or_else(|| error(None, "exactly one XYZ frame required"))?;
    let mut lines = xyz.lines();
    if lines.next().and_then(|s| s.trim().parse::<usize>().ok()) != Some(frame.atoms.len()) {
        return Err(error(None, "plain XYZ must begin with its atom count"));
    }
    if lines.next().is_some_and(|s| s.contains('=')) {
        return Err(error(None, "XYZ key/value comment metadata is unsupported"));
    }
    for (i, line) in lines.take(frame.atoms.len()).enumerate() {
        if line.split_whitespace().count() != 4 {
            return Err(error(Some(i), "plain XYZ requires exactly element x y z"));
        }
    }
    let elements: Vec<_> = frame.atoms.iter().map(|a| a.element).collect();
    let positions: Vec<_> = frame.atoms.iter().map(|a| [a.x, a.y, a.z]).collect();
    write_gamess_zmat(&elements, &positions, references, options)
}
