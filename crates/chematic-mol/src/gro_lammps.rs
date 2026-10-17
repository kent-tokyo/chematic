//! Explicit GRO nm/nm-per-ps to LAMMPS atomic geometry conversion.
//! Derived from the GROMACS GRO and LAMMPS units/triclinic-box manuals.
use crate::gro::{GroFrame, GroParseLimits, parse_gro_with_limits};
use crate::lammps_data::{
    LammpsAtom, LammpsAtomStyle, LammpsBox, LammpsData, LammpsMass, LammpsVelocity,
};

/// LAMMPS input-script units. Both use Angstroms and grams/mole.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GroLammpsUnits {
    Real,
    Metal,
}
impl GroLammpsUnits {
    /// nm/ps to Angstrom/fs (real) or Angstrom/ps (metal).
    pub fn velocity_scale(self) -> f64 {
        match self {
            Self::Real => 0.01,
            Self::Metal => 10.0,
        }
    }
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Real => "real",
            Self::Metal => "metal",
        }
    }
}
/// Invalid input, geometry profile, type assignment or numeric scaling.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GroLammpsError(pub String);
impl std::fmt::Display for GroLammpsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for GroLammpsError {}

/// Atomic coordinate document plus exact source and row identity.
/// Atom ID i+1 corresponds to source_frame.atoms[i], including wrapped/repeated
/// GRO IDs. No chemical topology or charges are inferred.
#[derive(Debug, Clone, PartialEq)]
pub struct GroLammpsConversion {
    pub data: LammpsData,
    pub units: GroLammpsUnits,
    pub source_frame: GroFrame,
}
fn scaled(value: f64, factor: f64) -> Result<f64, GroLammpsError> {
    let result = value * factor;
    if !result.is_finite() || (value != 0.0 && result == 0.0) {
        return Err(GroLammpsError(
            "unit conversion overflow or underflow".into(),
        ));
    }
    Ok(result)
}
fn triple(values: [f64; 3], factor: f64) -> Result<[f64; 3], GroLammpsError> {
    Ok([
        scaled(values[0], factor)?,
        scaled(values[1], factor)?,
        scaled(values[2], factor)?,
    ])
}
fn simulation_box(frame: &GroFrame) -> Result<LammpsBox, GroLammpsError> {
    let [a, b, c] = frame.box_vectors_nm;
    if a[1] != 0.0 || a[2] != 0.0 || b[2] != 0.0 {
        return Err(GroLammpsError(
            "only GRO restricted triclinic boxes are supported (v1y=v1z=v2z=0)".into(),
        ));
    }
    let hi = triple([a[0], b[1], c[2]], 10.0)?;
    let tilt = triple([b[0], c[0], c[1]], 10.0)?;
    let result = LammpsBox {
        lo: [0.0; 3],
        hi,
        tilt: (tilt != [0.0; 3]).then_some(tilt),
    };
    result.validate().map_err(GroLammpsError)?;
    Ok(result)
}

/// Convert one GRO frame to atomic data, with type IDs in source-row order and
/// dense masses (mass[i] belongs to type i+1). Units are mandatory.
/// Positive orthogonal/restricted-triclinic cells only; no rotation, shifting
/// or wrapping. Velocities must be on all rows or none. Limits: 16 MiB input,
/// 100,000 atoms, 1024-byte lines and 4096 types. Names/IDs remain in the sidecar.
pub fn gro_to_lammps_atomic(
    text: &str,
    atom_types: &[i64],
    masses: &[f64],
    units: GroLammpsUnits,
) -> Result<GroLammpsConversion, GroLammpsError> {
    if masses.is_empty()
        || masses.len() > 4096
        || masses.iter().any(|m| !m.is_finite() || *m <= 0.0)
    {
        return Err(GroLammpsError(
            "supply 1 to 4096 positive finite masses in grams/mole".into(),
        ));
    }
    let limits = GroParseLimits {
        max_input_bytes: 16 << 20,
        max_atoms: 100_000,
        max_line_len: 1024,
    };
    let source = parse_gro_with_limits(text, &limits).map_err(|e| GroLammpsError(e.to_string()))?;
    if source.atoms.len() != atom_types.len()
        || atom_types.iter().any(|&t| t < 1 || t > masses.len() as i64)
    {
        return Err(GroLammpsError(
            "atom_types must have one valid 1-based mass-table type per source row".into(),
        ));
    }
    let present = source
        .atoms
        .iter()
        .filter(|a| a.velocity_nm_per_ps.is_some())
        .count();
    if present != 0 && present != source.atoms.len() {
        return Err(GroLammpsError(
            "mixed present/missing GRO velocities cannot be converted".into(),
        ));
    }
    let mut data = LammpsData {
        counts: vec![
            ("atoms".into(), source.atoms.len() as i64),
            ("atom types".into(), masses.len() as i64),
        ],
        atom_style: LammpsAtomStyle::Atomic,
        simulation_box: simulation_box(&source)?,
        masses: masses
            .iter()
            .enumerate()
            .map(|(i, &mass)| LammpsMass {
                atom_type: i as i64 + 1,
                mass,
            })
            .collect(),
        atoms: Vec::with_capacity(source.atoms.len()),
        velocities: Vec::new(),
        bonds: Vec::new(),
        unparsed_sections: Vec::new(),
    };
    for (i, (atom, &atom_type)) in source.atoms.iter().zip(atom_types).enumerate() {
        let id = i as i64 + 1;
        let [x, y, z] = triple(atom.position_nm, 10.0)?;
        data.atoms.push(LammpsAtom {
            id,
            atom_type,
            molecule_id: None,
            charge: None,
            x,
            y,
            z,
            image: None,
        });
        if let Some(velocity) = atom.velocity_nm_per_ps {
            let [vx, vy, vz] = triple(velocity, units.velocity_scale())?;
            data.velocities.push(LammpsVelocity {
                atom_id: id,
                vx,
                vy,
                vz,
            });
        }
    }
    Ok(GroLammpsConversion {
        data,
        units,
        source_frame: source,
    })
}
