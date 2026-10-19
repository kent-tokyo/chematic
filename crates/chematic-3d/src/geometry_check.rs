//! Read-only diagnostics for caller-supplied Cartesian coordinates (Angstrom).
//!
//! This is a bounded geometric screen, not chemical sanitization, an energy
//! calculation, or a proof that a conformer is physically realizable. Reference
//! lengths reuse the embedding model; clashes use a caller-selected absolute
//! cutoff, including pairs in separate fragments. No coordinates are generated.

use chematic_core::{AtomIdx, BondIdx, BondOrder, Coords3D, Molecule, Point3};

use crate::dg_fft::ideal_bond_length;
use crate::stereo_constraints::{StereoVerification, verify_stereo};

pub const MAX_GEOMETRY_CHECK_ATOMS: usize = 2048;
/// Absolute coordinate-component bound, in Angstrom, for stable stereo arithmetic.
pub const MAX_GEOMETRY_CHECK_COORDINATE: f64 = 1_000_000.0;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeometryCheckConfig {
    /// Flag a bond when abs(actual-reference)/reference exceeds this value.
    /// Must be finite and in [0, 1); default 0.5.
    pub bond_relative_tolerance: f64,
    /// Flag unconnected atom pairs strictly closer than this many Angstrom.
    /// Must be finite and positive; default 1.2. Not a van der Waals test.
    pub nonbonded_min_distance: f64,
}

impl Default for GeometryCheckConfig {
    fn default() -> Self {
        Self {
            bond_relative_tolerance: 0.5,
            nonbonded_min_distance: 1.2,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GeometryCheckError {
    EmptyMolecule,
    AtomLimit,
    CoordinateCount { expected: usize, actual: usize },
    NonfiniteCoordinate { atom: AtomIdx },
    CoordinateRange { atom: AtomIdx },
    BondTolerance,
    ClashCutoff,
}

impl std::fmt::Display for GeometryCheckError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyMolecule => write!(f, "geometry check requires a nonempty molecule"),
            Self::AtomLimit => write!(
                f,
                "geometry check supports at most {MAX_GEOMETRY_CHECK_ATOMS} atoms"
            ),
            Self::CoordinateCount { expected, actual } => {
                write!(f, "expected {expected} coordinates, got {actual}")
            }
            Self::NonfiniteCoordinate { atom } => {
                write!(f, "non-finite coordinate at atom {}", atom.0)
            }
            Self::CoordinateRange { atom } => write!(
                f,
                "coordinate at atom {} exceeds +/-{MAX_GEOMETRY_CHECK_COORDINATE} Angstrom",
                atom.0
            ),
            Self::BondTolerance => {
                write!(f, "bond_relative_tolerance must be finite and in [0, 1)")
            }
            Self::ClashCutoff => write!(f, "nonbonded_min_distance must be finite and positive"),
        }
    }
}

impl std::error::Error for GeometryCheckError {}

fn distance(a: Point3, b: Point3) -> f64 {
    // Squaring tiny finite deltas can underflow to zero and falsely report a
    // clash for a still smaller caller cutoff. Scaled hypot retains the norm.
    libm::hypot(libm::hypot(a.x - b.x, a.y - b.y), a.z - b.z)
}

#[derive(Debug, Clone, PartialEq)]
pub struct BondGeometryViolation {
    pub bond: BondIdx,
    pub atoms: [AtomIdx; 2],
    pub distance: f64,
    pub reference_length: f64,
    pub relative_deviation: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct NonbondedClash {
    pub atoms: [AtomIdx; 2],
    pub distance: f64,
}

#[derive(Debug, Clone)]
pub struct GeometryCheckReport {
    /// Thresholds actually used, retained for reproducible batch decisions.
    pub config: GeometryCheckConfig,
    pub atom_count: usize,
    pub bonds_checked: usize,
    /// Unsupported bond orders or pseudoatom endpoints; never silently passed.
    pub skipped_bonds: Vec<BondIdx>,
    /// None when no bond had a supported reference length.
    pub max_relative_bond_deviation: Option<f64>,
    pub bond_violations: Vec<BondGeometryViolation>,
    /// All graph-unconnected pairs below the cutoff, including 1-3 pairs and
    /// separate fragments. Sorted by original atom indices.
    pub nonbonded_clashes: Vec<NonbondedClash>,
    /// Number of enhanced stereo groups whose group semantics were not checked.
    pub stereo_groups_not_evaluated: usize,
    /// Declared tetrahedral and directed E/Z checks, including unknown outcomes.
    /// Unspecified stereochemistry is not inferred. Other stereo classes and
    /// enhanced stereo-group semantics are outside this screen.
    pub stereo: StereoVerification,
}

/// Inspect coordinates in the original molecule's atom order without mutation.
///
/// Checks ordinary single/double/triple/aromatic bonds (including directional
/// single bonds). Dative, zero, quadruple and query bonds, and bonds touching
/// pseudoatoms, have no supported length reference and are reported as skipped.
/// Any graph bond still excludes its endpoints from the nonbonded clash screen.
/// Cost is O(N^2) with at most 2048 atoms; output can contain O(N^2) clashes.
/// Invalid shapes, non-finite values and values outside +/-1e6 A are errors;
/// finite but collapsed/stretched coordinates return diagnostics, not errors.
pub fn check_geometry(
    mol: &Molecule,
    coords: &Coords3D,
    config: &GeometryCheckConfig,
) -> Result<GeometryCheckReport, GeometryCheckError> {
    let n = mol.atom_count();
    if n == 0 {
        return Err(GeometryCheckError::EmptyMolecule);
    }
    if n > MAX_GEOMETRY_CHECK_ATOMS {
        return Err(GeometryCheckError::AtomLimit);
    }
    if coords.atom_count() != n {
        return Err(GeometryCheckError::CoordinateCount {
            expected: n,
            actual: coords.atom_count(),
        });
    }
    if !config.bond_relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&config.bond_relative_tolerance)
    {
        return Err(GeometryCheckError::BondTolerance);
    }
    if !config.nonbonded_min_distance.is_finite() || config.nonbonded_min_distance <= 0.0 {
        return Err(GeometryCheckError::ClashCutoff);
    }
    for (i, p) in coords.points.iter().enumerate() {
        let atom = AtomIdx(i as u32);
        if !p.is_finite() {
            return Err(GeometryCheckError::NonfiniteCoordinate { atom });
        }
        if [p.x, p.y, p.z]
            .iter()
            .any(|v| v.abs() > MAX_GEOMETRY_CHECK_COORDINATE)
        {
            return Err(GeometryCheckError::CoordinateRange { atom });
        }
    }
    let mut report = GeometryCheckReport {
        config: *config,
        atom_count: n,
        bonds_checked: 0,
        skipped_bonds: Vec::new(),
        max_relative_bond_deviation: None,
        bond_violations: Vec::new(),
        nonbonded_clashes: Vec::new(),
        stereo_groups_not_evaluated: mol.stereo_groups().len(),
        stereo: verify_stereo(mol, coords),
    };
    // Dense, bounded adjacency avoids a graph lookup per potential clash pair.
    let mut connected = vec![false; n * n];
    for (idx, bond) in mol.bonds() {
        let [a, b] = [bond.atom1, bond.atom2];
        connected[a.0 as usize * n + b.0 as usize] = true;
        connected[b.0 as usize * n + a.0 as usize] = true;
        let supported = matches!(
            bond.order,
            BondOrder::Single
                | BondOrder::Up
                | BondOrder::Down
                | BondOrder::Double
                | BondOrder::Triple
                | BondOrder::Aromatic
        ) && !mol.atom(a).wildcard
            && !mol.atom(b).wildcard;
        if !supported {
            report.skipped_bonds.push(idx);
            continue;
        }
        let reference_length = ideal_bond_length(mol, a, b);
        let distance = distance(coords.get(a), coords.get(b));
        let relative_deviation = (distance - reference_length).abs() / reference_length;
        report.bonds_checked += 1;
        report.max_relative_bond_deviation = Some(
            report
                .max_relative_bond_deviation
                .unwrap_or(0.0)
                .max(relative_deviation),
        );
        if relative_deviation > config.bond_relative_tolerance {
            report.bond_violations.push(BondGeometryViolation {
                bond: idx,
                atoms: [a, b],
                distance,
                reference_length,
                relative_deviation,
            });
        }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if connected[i * n + j] {
                continue;
            }
            let atoms = [AtomIdx(i as u32), AtomIdx(j as u32)];
            let distance = distance(coords.get(atoms[0]), coords.get(atoms[1]));
            if distance < config.nonbonded_min_distance {
                report
                    .nonbonded_clashes
                    .push(NonbondedClash { atoms, distance });
            }
        }
    }
    Ok(report)
}
