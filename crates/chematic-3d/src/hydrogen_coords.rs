//! Add explicit-hydrogen coordinates to an existing heavy-atom conformer.
//!
//! The placement rules follow RDKit's BSD-licensed
//! `MolOps::setTerminalAtomCoords()` (RDKit 2026.09.1): hydrogens are appended
//! in atom order, each new atom sees the neighbors that existed when it was
//! appended, and SP/SP2/SP3 geometry determines its initial direction.  This
//! is substantially cheaper than embedding all hydrogens with distance
//! geometry and provides a sound starting point for a full-atom force-field
//! refinement.

use chematic_core::{AtomIdx, Coords3D, Molecule, Point3};

use crate::minimize::{
    ForceFieldBridgeError, ForceFieldPolicy, MinimizeConfig, minimize_with_policy,
};

const ZERO_TOL_SQ: f64 = 1.0e-16;

/// Failure to construct an explicit-hydrogen conformer from heavy coordinates.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum HydrogenCoordinateError {
    CoordinateCount { expected: usize, got: usize },
    MoleculeView(String),
    InvalidAddedHydrogen { atom: usize },
}

impl core::fmt::Display for HydrogenCoordinateError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::CoordinateCount { expected, got } => {
                write!(f, "expected {expected} atom coordinates, got {got}")
            }
            Self::MoleculeView(reason) => write!(f, "cannot build RDKit molecule view: {reason}"),
            Self::InvalidAddedHydrogen { atom } => {
                write!(f, "added hydrogen {atom} is not terminal")
            }
        }
    }
}

impl std::error::Error for HydrogenCoordinateError {}

/// Full-atom MMFF94 refinement of a heavy-atom conformer.
#[derive(Debug, Clone)]
pub struct ExplicitHydrogenMmffResult {
    /// Refined coordinates for the original molecule's atoms.
    pub heavy_coords: Coords3D,
    /// Refined coordinates for the original atoms followed by added H atoms.
    pub all_atom_coords: Coords3D,
    pub energy: f64,
    pub converged: bool,
    pub iterations: usize,
    pub max_residual_force: f64,
}

/// Typed failure from coordinate materialization or MMFF94 refinement.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub enum ExplicitHydrogenMmffError {
    HydrogenCoordinates(HydrogenCoordinateError),
    ForceField(ForceFieldBridgeError),
}

impl core::fmt::Display for ExplicitHydrogenMmffError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::HydrogenCoordinates(error) => error.fmt(f),
            Self::ForceField(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for ExplicitHydrogenMmffError {}

impl From<HydrogenCoordinateError> for ExplicitHydrogenMmffError {
    fn from(error: HydrogenCoordinateError) -> Self {
        Self::HydrogenCoordinates(error)
    }
}

impl From<ForceFieldBridgeError> for ExplicitHydrogenMmffError {
    fn from(error: ForceFieldBridgeError) -> Self {
        Self::ForceField(error)
    }
}

fn norm_sq(v: Point3) -> f64 {
    v.dot(&v)
}

fn normalized(v: Point3) -> Option<Point3> {
    (norm_sq(v) >= ZERO_TOL_SQ).then(|| v.scale(1.0 / v.norm()))
}

/// RDKit `Point3D::getPerpendicular()`.
fn perpendicular(v: Point3) -> Point3 {
    let raw = if v.x != 0.0 {
        if v.y != 0.0 {
            Point3::new(v.y, -v.x, 0.0)
        } else if v.z != 0.0 {
            Point3::new(v.z, 0.0, -v.x)
        } else {
            Point3::new(0.0, 1.0, 0.0)
        }
    } else if v.y != 0.0 {
        if v.z != 0.0 {
            Point3::new(0.0, v.z, -v.y)
        } else {
            Point3::new(1.0, 0.0, 0.0)
        }
    } else {
        Point3::new(1.0, 0.0, 0.0)
    };
    raw.normalize()
}

fn rotate(v: Point3, axis: Point3, angle: f64) -> Point3 {
    let axis = axis.normalize();
    let (sin, cos) = angle.sin_cos();
    v.scale(cos)
        .add(&axis.cross(&v).scale(sin))
        .add(&axis.scale(axis.dot(&v) * (1.0 - cos)))
}

fn point_away(parent: Point3, neighbor: Point3) -> Option<Point3> {
    normalized(parent.sub(&neighbor))
}

fn bond_length(atomic_number: u32) -> f64 {
    chematic_chem::rdkit_single_bond_radius(1)
        + chematic_chem::rdkit_single_bond_radius(atomic_number as u8)
}

fn direction_from_three_neighbors(
    view: &chematic_smiles::RdkitMolView,
    coords: &Coords3D,
    parent: usize,
    neighbors: &[usize],
) -> Option<Point3> {
    let parent_pos = coords.points[parent];
    let v1 = point_away(parent_pos, coords.points[neighbors[0]])?;
    let v2 = point_away(parent_pos, coords.points[neighbors[1]])?;
    let v3 = point_away(parent_pos, coords.points[neighbors[2]])?;
    let mut direction = if v3.dot(&v1.cross(&v2)).abs() < 0.1 {
        let mut normal = v1.cross(&v2);
        if norm_sq(normal) < ZERO_TOL_SQ {
            normal = v1.cross(&v3).scale(-1.0);
        }
        if norm_sq(normal) < ZERO_TOL_SQ {
            normal = v2.cross(&v3);
        }
        if norm_sq(normal) < ZERO_TOL_SQ {
            return None;
        }
        if matches!(view.atoms[parent].chiral_tag, 1 | 2) {
            let volume = normal.sub(&v3).dot(&v1.sub(&v3).cross(&v2.sub(&v3)));
            if (view.atoms[parent].chiral_tag == 2 && volume < 0.0)
                || (view.atoms[parent].chiral_tag == 1 && volume > 0.0)
            {
                normal = normal.scale(-1.0);
            }
        }
        normal
    } else {
        v1.add(&v2).add(&v3)
    };
    direction = normalized(direction).unwrap_or(Point3::new(0.0, 0.0, 1.0));
    Some(direction)
}

fn initial_hydrogen_direction(
    view: &chematic_smiles::RdkitMolView,
    coords: &Coords3D,
    parent: usize,
    hydrogen: usize,
) -> Option<Point3> {
    let parent_pos = coords.points[parent];
    let others: Vec<usize> = view
        .neighbors(parent)
        .filter(|&neighbor| neighbor != hydrogen && neighbor < hydrogen)
        .collect();
    match others.len() + 1 {
        1 => Some(Point3::new(0.0, 0.0, 1.0)),
        2 => {
            let neighbor = others[0];
            let away = point_away(parent_pos, coords.points[neighbor])?;
            match view.atoms[parent].hybridization {
                4 => Some(rotate(
                    away,
                    perpendicular(away),
                    (180.0_f64 - 109.471).to_radians(),
                )),
                3 => {
                    let mut axis = perpendicular(away);
                    if view.degree(neighbor) > 1 {
                        let bond = view
                            .bond_between(parent, neighbor)
                            .expect("neighbors are bonded");
                        let bond = &view.bonds[bond];
                        if (bond.is_aromatic || bond.bond_type == 2 || bond.is_conjugated)
                            && let Some(next) =
                                view.neighbors(neighbor).find(|&atom| atom != parent)
                        {
                            let along = coords.points[next].sub(&coords.points[neighbor]);
                            let cross = along.cross(&away);
                            if norm_sq(cross) >= ZERO_TOL_SQ {
                                axis = cross.normalize();
                            }
                        }
                    }
                    Some(rotate(away, axis, 60.0_f64.to_radians()))
                }
                _ => Some(away),
            }
        }
        3 => {
            let v1 = point_away(parent_pos, coords.points[others[0]])?;
            let v2 = point_away(parent_pos, coords.points[others[1]])?;
            let mut direction = normalized(v1.add(&v2))?;
            if view.atoms[parent].hybridization == 4 {
                let axis = v1.cross(&v2).cross(&direction);
                if norm_sq(axis) >= ZERO_TOL_SQ {
                    direction = rotate(direction, axis, (109.471_f64 / 2.0).to_radians());
                }
            }
            Some(direction)
        }
        4 => direction_from_three_neighbors(view, coords, parent, &others[..3]),
        _ => Some(Point3::new(0.0, 0.0, 1.0)),
    }
}

/// Return `mol` with all implicit hydrogens made explicit and coordinates
/// appended for those hydrogens. Original atom indices and coordinates are
/// unchanged.
pub fn add_hydrogens_with_coords(
    mol: &Molecule,
    coords: &Coords3D,
) -> Result<(Molecule, Coords3D), HydrogenCoordinateError> {
    let original_atoms = mol.atom_count();
    if coords.atom_count() != original_atoms {
        return Err(HydrogenCoordinateError::CoordinateCount {
            expected: original_atoms,
            got: coords.atom_count(),
        });
    }

    let expanded = chematic_chem::add_hydrogens(mol);
    let view = chematic_smiles::rdkit_mol_view(&expanded)
        .map_err(|error| HydrogenCoordinateError::MoleculeView(error.to_string()))?;
    let mut out = Coords3D::new_zeroed(expanded.atom_count());
    out.points[..original_atoms].copy_from_slice(&coords.points);

    // `add_hydrogens()` appends H atoms in parent-atom order. RDKit places a
    // hydrogen immediately after appending it; filtering out later H atoms
    // reproduces the degree and neighbor set visible at that exact moment.
    for h in original_atoms..expanded.atom_count() {
        let active: Vec<usize> = view.neighbors(h).filter(|&n| n < h).collect();
        if active.len() != 1 {
            return Err(HydrogenCoordinateError::InvalidAddedHydrogen { atom: h });
        }
        let parent = active[0];
        let parent_pos = out.points[parent];
        let length = bond_length(view.atoms[parent].atomic_num);
        let Some(direction) = initial_hydrogen_direction(&view, &out, parent, h) else {
            out.points[h] = parent_pos;
            continue;
        };
        out.points[h] = parent_pos.add(&direction.scale(length));
    }
    Ok((expanded, out))
}

/// Add hydrogens to a heavy-atom conformer using RDKit-compatible initial
/// coordinates, then minimize the full molecule with MMFF94.
///
/// This is designed for conformer ranking: unlike a heavy-only MMFF pass, the
/// returned energy includes X-H bonded terms, hydrogen torsions, and all-atom
/// non-bonded interactions. The caller still receives coordinates aligned to
/// the original molecule atom order.
pub fn refine_mmff94_with_explicit_hydrogens(
    mol: &Molecule,
    coords: &Coords3D,
    max_steps: usize,
) -> Result<ExplicitHydrogenMmffResult, ExplicitHydrogenMmffError> {
    refine_mmff94_with_explicit_hydrogens_variant(mol, coords, max_steps, 0)
}

fn splitmix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9E37_79B9_7F4A_7C15);
    value = (value ^ (value >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    value ^ (value >> 31)
}

fn rotate_heteroatom_hydrogens(
    mol: &Molecule,
    coords: &mut Coords3D,
    first_hydrogen: usize,
    variant: u64,
) {
    if variant == 0 {
        return;
    }
    let golden = 0.618_033_988_749_894_9;
    for h in first_hydrogen..mol.atom_count() {
        let Some(parent) = mol
            .neighbors(AtomIdx(h as u32))
            .next()
            .map(|(atom, _)| atom.0 as usize)
        else {
            continue;
        };
        if !matches!(
            mol.atom(AtomIdx(parent as u32)).element.atomic_number(),
            7 | 8 | 15 | 16 | 34 | 52
        ) {
            continue;
        }
        let Some(anchor) = mol
            .neighbors(AtomIdx(parent as u32))
            .map(|(atom, _)| atom.0 as usize)
            .find(|&atom| atom != h && mol.atom(AtomIdx(atom as u32)).element.atomic_number() != 1)
        else {
            continue;
        };
        let axis = coords.points[anchor].sub(&coords.points[parent]);
        if norm_sq(axis) < ZERO_TOL_SQ {
            continue;
        }
        let phase = splitmix64(parent as u64) as f64 / u64::MAX as f64;
        let turns = (phase + variant as f64 * golden).fract();
        let relative = coords.points[h].sub(&coords.points[parent]);
        coords.points[h] =
            coords.points[parent].add(&rotate(relative, axis, turns * core::f64::consts::TAU));
    }
}

/// Variant of [`refine_mmff94_with_explicit_hydrogens`] that samples a
/// deterministic initial orientation for terminal heteroatom hydrogens.
/// Variant zero preserves RDKit's default placement; non-zero variants cover
/// otherwise hidden O-H/N-H/S-H rotamers without full explicit-H embedding.
pub fn refine_mmff94_with_explicit_hydrogens_variant(
    mol: &Molecule,
    coords: &Coords3D,
    max_steps: usize,
    hydrogen_variant: u64,
) -> Result<ExplicitHydrogenMmffResult, ExplicitHydrogenMmffError> {
    let heavy_atoms = mol.atom_count();
    let (expanded, mut expanded_coords) = add_hydrogens_with_coords(mol, coords)?;
    rotate_heteroatom_hydrogens(
        &expanded,
        &mut expanded_coords,
        heavy_atoms,
        hydrogen_variant,
    );
    let result = minimize_with_policy(
        &expanded,
        expanded_coords,
        ForceFieldPolicy::Mmff94BondAngleStrict,
        &MinimizeConfig {
            max_steps,
            ..MinimizeConfig::default()
        },
    )?;
    let energy = result.energy_after.total();
    let heavy_coords = Coords3D {
        points: result.coords.points[..heavy_atoms].to_vec(),
    };
    Ok(ExplicitHydrogenMmffResult {
        heavy_coords,
        all_atom_coords: result.coords,
        energy,
        converged: result.converged,
        iterations: result.iterations,
        max_residual_force: result.max_residual_force,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn typed_errors_keep_actionable_context() {
        let mol = chematic_smiles::parse("CC").unwrap();
        let error = match add_hydrogens_with_coords(&mol, &Coords3D::new_zeroed(1)) {
            Ok(_) => panic!("coordinate count mismatch must fail"),
            Err(error) => error,
        };
        assert_eq!(error.to_string(), "expected 2 atom coordinates, got 1");
        let wrapped: ExplicitHydrogenMmffError = error.into();
        assert!(wrapped.to_string().contains("expected 2 atom coordinates"));

        let view = HydrogenCoordinateError::MoleculeView("bad graph".into());
        assert_eq!(
            view.to_string(),
            "cannot build RDKit molecule view: bad graph"
        );
        let terminal = HydrogenCoordinateError::InvalidAddedHydrogen { atom: 7 };
        assert_eq!(terminal.to_string(), "added hydrogen 7 is not terminal");

        let force_field: ExplicitHydrogenMmffError =
            ForceFieldBridgeError::UnsupportedAtomType("test atom".into()).into();
        assert!(force_field.to_string().contains("test atom"));
    }

    #[test]
    fn vector_helpers_cover_degenerate_and_axis_aligned_inputs() {
        assert!(normalized(Point3::zero()).is_none());
        for vector in [
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 1.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::zero(),
        ] {
            let normal = perpendicular(vector);
            assert!((normal.norm() - 1.0).abs() < 1.0e-12);
            assert!(vector == Point3::zero() || normal.dot(&vector).abs() < 1.0e-12);
        }
    }

    #[test]
    fn placement_covers_sp_sp2_sp3_and_successive_hydrogens() {
        for smiles in ["C", "C=C", "C#C", "N", "O"] {
            let mol = chematic_smiles::parse(smiles).unwrap();
            let coords = crate::generate_coords(&mol);
            let (expanded, placed) = add_hydrogens_with_coords(&mol, &coords).unwrap();
            assert!(expanded.atom_count() >= mol.atom_count());
            assert!(placed.is_finite());
            for h in mol.atom_count()..expanded.atom_count() {
                let parent = expanded.neighbors(AtomIdx(h as u32)).next().unwrap().0;
                assert!(placed.points[h].distance(&placed.get(parent)) > 0.5);
            }
        }
    }

    #[test]
    fn coplanar_three_neighbor_direction_is_typed_as_degenerate() {
        let expanded = chematic_chem::add_hydrogens(&chematic_smiles::parse("C").unwrap());
        let view = chematic_smiles::rdkit_mol_view(&expanded).unwrap();
        let neighbors: Vec<_> = view.neighbors(0).take(3).collect();
        let mut coords = Coords3D::new_zeroed(expanded.atom_count());
        coords.points[neighbors[0]] = Point3::new(1.0, 0.0, 0.0);
        coords.points[neighbors[1]] = Point3::new(2.0, 0.0, 0.0);
        coords.points[neighbors[2]] = Point3::new(3.0, 0.0, 0.0);
        assert!(direction_from_three_neighbors(&view, &coords, 0, &neighbors).is_none());
    }

    #[test]
    fn unsupported_mmff_type_is_returned_as_force_field_error() {
        let mol = chematic_smiles::parse("C[U]").unwrap();
        let coords = crate::generate_coords(&mol);
        let error = refine_mmff94_with_explicit_hydrogens(&mol, &coords, 1).unwrap_err();
        assert!(
            matches!(
                &error,
                ExplicitHydrogenMmffError::ForceField(ForceFieldBridgeError::UnsupportedAtomType(
                    _
                ))
            ),
            "unexpected error: {error:?}"
        );
    }

    #[test]
    fn appends_finite_hydrogen_coordinates_without_moving_heavy_atoms() {
        let mol = chematic_smiles::parse("C[C@H](O)F").unwrap();
        let coords = Coords3D {
            points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.5, 0.0, 0.0),
                Point3::new(2.0, 1.3, 0.0),
                Point3::new(2.0, -0.6, 1.1),
            ],
        };
        let (expanded, with_h) = add_hydrogens_with_coords(&mol, &coords).unwrap();
        assert!(expanded.atom_count() > mol.atom_count());
        assert_eq!(&with_h.points[..mol.atom_count()], coords.points.as_slice());
        assert!(with_h.is_finite());
        let expected_h = [
            [-0.36666267303974204, 0.0, 1.0370913576919594],
            [-0.366668663480129, 0.8981450162261428, -0.5185456788459798],
            [-0.366668663480129, -0.8981450162261428, -0.5185456788459796],
            [1.8009665647921922, -0.5438354187481101, -0.907558353050786],
            [2.118461805906262, 1.608000695356281, 0.9333822219227633],
        ];
        assert_eq!(expanded.atom_count(), mol.atom_count() + expected_h.len());
        for (point, expected) in with_h.points[mol.atom_count()..].iter().zip(expected_h) {
            assert!((point.x - expected[0]).abs() < 1.0e-12);
            assert!((point.y - expected[1]).abs() < 1.0e-12);
            assert!((point.z - expected[2]).abs() < 1.0e-12);
        }
        for h in mol.atom_count()..expanded.atom_count() {
            let parent = expanded
                .neighbors(AtomIdx(h as u32))
                .next()
                .expect("terminal H")
                .0;
            assert!(with_h.points[h].distance(&with_h.get(parent)) > 0.9);
        }
    }

    #[test]
    fn full_atom_mmff_refinement_returns_original_atom_shape_and_finite_energy() {
        let mol = chematic_smiles::parse("CCO").unwrap();
        let coords = crate::generate_coords(&mol);
        let result = refine_mmff94_with_explicit_hydrogens(&mol, &coords, 200).unwrap();
        assert_eq!(result.heavy_coords.atom_count(), mol.atom_count());
        assert!(result.all_atom_coords.atom_count() > mol.atom_count());
        assert!(result.energy.is_finite());
        assert!(result.max_residual_force.is_finite());
    }

    #[test]
    fn hydrogen_variants_are_deterministic_and_preserve_heavy_atom_order() {
        let mol = chematic_smiles::parse("CC(=O)O").unwrap();
        let coords = crate::generate_coords(&mol);
        let a = refine_mmff94_with_explicit_hydrogens_variant(&mol, &coords, 0, 3).unwrap();
        let b = refine_mmff94_with_explicit_hydrogens_variant(&mol, &coords, 0, 3).unwrap();
        let base = refine_mmff94_with_explicit_hydrogens_variant(&mol, &coords, 0, 0).unwrap();
        assert_eq!(a.all_atom_coords, b.all_atom_coords);
        assert_eq!(a.heavy_coords.atom_count(), mol.atom_count());
        assert_ne!(a.all_atom_coords, base.all_atom_coords);
    }
}
