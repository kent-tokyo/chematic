//! User-supplied internal-coordinate constraints for distance geometry.
//!
//! Unlike [`crate::constraints`], which projects coordinates after they have
//! been generated, these constraints are converted to pair-distance bounds
//! before triangle smoothing.  This keeps the complete bounds matrix
//! self-consistent and makes the constraints part of the embedding itself.

use std::collections::HashMap;

use chematic_core::Molecule;

use super::dg::BoundsMatrix;

/// A constrained bond length, in angstrom.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BondLengthConstraint {
    pub atom1: usize,
    pub atom2: usize,
    pub lower: f64,
    pub upper: f64,
}

impl BondLengthConstraint {
    pub fn around(atom1: usize, atom2: usize, target: f64, tolerance: f64) -> Self {
        Self {
            atom1,
            atom2,
            lower: target - tolerance,
            upper: target + tolerance,
        }
    }
}

/// A constrained valence angle, in radians.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BondAngleConstraint {
    pub atom1: usize,
    pub center: usize,
    pub atom2: usize,
    pub lower: f64,
    pub upper: f64,
}

impl BondAngleConstraint {
    pub fn around_degrees(
        atom1: usize,
        center: usize,
        atom2: usize,
        target_degrees: f64,
        tolerance_degrees: f64,
    ) -> Self {
        Self {
            atom1,
            center,
            atom2,
            lower: (target_degrees - tolerance_degrees).to_radians(),
            upper: (target_degrees + tolerance_degrees).to_radians(),
        }
    }
}

/// A constrained dihedral angle, in radians.
///
/// Intervals may cross the conventional `-pi`/`pi` seam (for example
/// `180 +/- 5` degrees); they may span at most one full revolution.
/// A distance-bounds matrix encodes the cosine-equivalent 1-4 separation, so
/// `+phi` and `-phi` remain mirror-equivalent unless molecular stereo or the
/// later ETKDG torsion stage selects one sign.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TorsionAngleConstraint {
    pub atom1: usize,
    pub atom2: usize,
    pub atom3: usize,
    pub atom4: usize,
    pub lower: f64,
    pub upper: f64,
}

impl TorsionAngleConstraint {
    pub fn around_degrees(
        atom1: usize,
        atom2: usize,
        atom3: usize,
        atom4: usize,
        target_degrees: f64,
        tolerance_degrees: f64,
    ) -> Self {
        Self {
            atom1,
            atom2,
            atom3,
            atom4,
            lower: (target_degrees - tolerance_degrees).to_radians(),
            upper: (target_degrees + tolerance_degrees).to_radians(),
        }
    }
}

/// Internal coordinates to inject before distance-bound triangle smoothing.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct InternalCoordinateConstraints {
    pub bond_lengths: Vec<BondLengthConstraint>,
    pub bond_angles: Vec<BondAngleConstraint>,
    pub torsions: Vec<TorsionAngleConstraint>,
}

/// Invalid or mutually inconsistent user constraints.
#[derive(Clone, Debug, PartialEq)]
pub enum InternalCoordinateConstraintError {
    AtomOutOfRange {
        atom: usize,
        atom_count: usize,
    },
    RepeatedAtom {
        kind: &'static str,
    },
    MissingBond {
        atom1: usize,
        atom2: usize,
    },
    NonFinite {
        kind: &'static str,
    },
    InvalidInterval {
        kind: &'static str,
        lower: f64,
        upper: f64,
    },
    ConflictingPair {
        atom1: usize,
        atom2: usize,
    },
    CrossFragment {
        kind: &'static str,
    },
    TriangleSmoothingFailed,
}

impl core::fmt::Display for InternalCoordinateConstraintError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::AtomOutOfRange { atom, atom_count } => {
                write!(
                    f,
                    "atom index {atom} is outside molecule with {atom_count} atoms"
                )
            }
            Self::RepeatedAtom { kind } => write!(f, "{kind} constraint repeats an atom"),
            Self::MissingBond { atom1, atom2 } => {
                write!(f, "constraint requires missing bond {atom1}-{atom2}")
            }
            Self::NonFinite { kind } => write!(f, "{kind} constraint contains a non-finite value"),
            Self::InvalidInterval { kind, lower, upper } => write!(
                f,
                "{kind} constraint has invalid interval [{lower}, {upper}]"
            ),
            Self::ConflictingPair { atom1, atom2 } => {
                write!(
                    f,
                    "constraints for atom pair {atom1}-{atom2} do not overlap"
                )
            }
            Self::CrossFragment { kind } => {
                write!(f, "{kind} constraint spans disconnected fragments")
            }
            Self::TriangleSmoothingFailed => {
                write!(f, "constraints make the distance bounds inconsistent")
            }
        }
    }
}

impl std::error::Error for InternalCoordinateConstraintError {}

impl InternalCoordinateConstraints {
    pub(crate) fn validate(&self, mol: &Molecule) -> Result<(), InternalCoordinateConstraintError> {
        let n = mol.atom_count();
        for constraint in &self.bond_lengths {
            validate_atoms(&[constraint.atom1, constraint.atom2], n, "bond-length")?;
            validate_interval(
                constraint.lower,
                constraint.upper,
                "bond-length",
                0.0,
                f64::INFINITY,
            )?;
            require_bond(mol, constraint.atom1, constraint.atom2)?;
        }
        for constraint in &self.bond_angles {
            validate_atoms(
                &[constraint.atom1, constraint.center, constraint.atom2],
                n,
                "bond-angle",
            )?;
            validate_interval(
                constraint.lower,
                constraint.upper,
                "bond-angle",
                0.0,
                core::f64::consts::PI,
            )?;
            require_bond(mol, constraint.atom1, constraint.center)?;
            require_bond(mol, constraint.center, constraint.atom2)?;
        }
        for constraint in &self.torsions {
            validate_atoms(
                &[
                    constraint.atom1,
                    constraint.atom2,
                    constraint.atom3,
                    constraint.atom4,
                ],
                n,
                "torsion",
            )?;
            validate_torsion_interval(constraint.lower, constraint.upper)?;
            require_bond(mol, constraint.atom1, constraint.atom2)?;
            require_bond(mol, constraint.atom2, constraint.atom3)?;
            require_bond(mol, constraint.atom3, constraint.atom4)?;
        }
        Ok(())
    }

    pub(crate) fn remap_fragment(
        &self,
        source_atoms: &[usize],
    ) -> Result<Self, InternalCoordinateConstraintError> {
        let map: HashMap<usize, usize> = source_atoms
            .iter()
            .enumerate()
            .map(|(local, &source)| (source, local))
            .collect();
        let mut out = Self::default();
        for constraint in &self.bond_lengths {
            if let Some([atom1, atom2]) =
                remap_or_skip(&map, [constraint.atom1, constraint.atom2], "bond-length")?
            {
                out.bond_lengths.push(BondLengthConstraint {
                    atom1,
                    atom2,
                    ..*constraint
                });
            }
        }
        for constraint in &self.bond_angles {
            if let Some([atom1, center, atom2]) = remap_or_skip(
                &map,
                [constraint.atom1, constraint.center, constraint.atom2],
                "bond-angle",
            )? {
                out.bond_angles.push(BondAngleConstraint {
                    atom1,
                    center,
                    atom2,
                    ..*constraint
                });
            }
        }
        for constraint in &self.torsions {
            if let Some([atom1, atom2, atom3, atom4]) = remap_or_skip(
                &map,
                [
                    constraint.atom1,
                    constraint.atom2,
                    constraint.atom3,
                    constraint.atom4,
                ],
                "torsion",
            )? {
                out.torsions.push(TorsionAngleConstraint {
                    atom1,
                    atom2,
                    atom3,
                    atom4,
                    ..*constraint
                });
            }
        }
        Ok(out)
    }

    pub(crate) fn apply(
        &self,
        bounds: &mut BoundsMatrix,
    ) -> Result<(), InternalCoordinateConstraintError> {
        let mut user_pairs: HashMap<(usize, usize), (f64, f64)> = HashMap::new();
        for constraint in &self.bond_lengths {
            set_user_pair(
                bounds,
                &mut user_pairs,
                constraint.atom1,
                constraint.atom2,
                constraint.lower,
                constraint.upper,
            )?;
        }
        for constraint in &self.bond_angles {
            let r1 = midpoint(bounds, constraint.atom1, constraint.center);
            let r2 = midpoint(bounds, constraint.center, constraint.atom2);
            let lower = distance13(r1, r2, constraint.lower);
            let upper = distance13(r1, r2, constraint.upper);
            set_user_pair(
                bounds,
                &mut user_pairs,
                constraint.atom1,
                constraint.atom2,
                lower.min(upper),
                lower.max(upper),
            )?;
        }
        for constraint in &self.torsions {
            let r1 = midpoint(bounds, constraint.atom1, constraint.atom2);
            let r2 = midpoint(bounds, constraint.atom2, constraint.atom3);
            let r3 = midpoint(bounds, constraint.atom3, constraint.atom4);
            let angle1 =
                inferred_angle(bounds, constraint.atom1, constraint.atom2, constraint.atom3);
            let angle2 =
                inferred_angle(bounds, constraint.atom2, constraint.atom3, constraint.atom4);
            let candidates = torsion_candidates(constraint.lower, constraint.upper);
            let mut lower = f64::INFINITY;
            let mut upper: f64 = 0.0;
            for torsion in candidates {
                let distance = distance14(r1, r2, r3, angle1, angle2, torsion);
                lower = lower.min(distance);
                upper = upper.max(distance);
            }
            // A point torsion would otherwise create a zero-width interval,
            // which the distance-geometry engine intentionally rejects.
            let epsilon = 1.0e-6;
            set_user_pair(
                bounds,
                &mut user_pairs,
                constraint.atom1,
                constraint.atom4,
                (lower - epsilon).max(1.0e-6),
                upper + epsilon,
            )?;
        }
        Ok(())
    }
}

fn validate_atoms(
    atoms: &[usize],
    atom_count: usize,
    kind: &'static str,
) -> Result<(), InternalCoordinateConstraintError> {
    for &atom in atoms {
        if atom >= atom_count {
            return Err(InternalCoordinateConstraintError::AtomOutOfRange { atom, atom_count });
        }
    }
    for i in 0..atoms.len() {
        if atoms[i + 1..].contains(&atoms[i]) {
            return Err(InternalCoordinateConstraintError::RepeatedAtom { kind });
        }
    }
    Ok(())
}

fn validate_interval(
    lower: f64,
    upper: f64,
    kind: &'static str,
    minimum: f64,
    maximum: f64,
) -> Result<(), InternalCoordinateConstraintError> {
    if !lower.is_finite() || !upper.is_finite() {
        return Err(InternalCoordinateConstraintError::NonFinite { kind });
    }
    if lower < minimum || upper > maximum || lower > upper {
        return Err(InternalCoordinateConstraintError::InvalidInterval { kind, lower, upper });
    }
    Ok(())
}

fn validate_torsion_interval(
    lower: f64,
    upper: f64,
) -> Result<(), InternalCoordinateConstraintError> {
    if !lower.is_finite() || !upper.is_finite() {
        return Err(InternalCoordinateConstraintError::NonFinite { kind: "torsion" });
    }
    if lower > upper || upper - lower > 2.0 * core::f64::consts::PI {
        return Err(InternalCoordinateConstraintError::InvalidInterval {
            kind: "torsion",
            lower,
            upper,
        });
    }
    Ok(())
}

fn require_bond(
    mol: &Molecule,
    atom1: usize,
    atom2: usize,
) -> Result<(), InternalCoordinateConstraintError> {
    if mol
        .bond_between(
            chematic_core::AtomIdx(atom1 as u32),
            chematic_core::AtomIdx(atom2 as u32),
        )
        .is_none()
    {
        return Err(InternalCoordinateConstraintError::MissingBond { atom1, atom2 });
    }
    Ok(())
}

fn remap_or_skip<const N: usize>(
    map: &HashMap<usize, usize>,
    atoms: [usize; N],
    kind: &'static str,
) -> Result<Option<[usize; N]>, InternalCoordinateConstraintError> {
    let present = atoms.map(|atom| map.get(&atom).copied());
    if present.iter().all(Option::is_none) {
        return Ok(None);
    }
    if present.iter().any(Option::is_none) {
        return Err(InternalCoordinateConstraintError::CrossFragment { kind });
    }
    Ok(Some(present.map(Option::unwrap)))
}

fn set_user_pair(
    bounds: &mut BoundsMatrix,
    user_pairs: &mut HashMap<(usize, usize), (f64, f64)>,
    atom1: usize,
    atom2: usize,
    lower: f64,
    upper: f64,
) -> Result<(), InternalCoordinateConstraintError> {
    let pair = if atom1 < atom2 {
        (atom1, atom2)
    } else {
        (atom2, atom1)
    };
    let (lower, upper) = if let Some(&(old_lower, old_upper)) = user_pairs.get(&pair) {
        (old_lower.max(lower), old_upper.min(upper))
    } else {
        (lower, upper)
    };
    if !lower.is_finite() || !upper.is_finite() || lower <= 0.0 || lower > upper {
        return Err(InternalCoordinateConstraintError::ConflictingPair {
            atom1: pair.0,
            atom2: pair.1,
        });
    }
    user_pairs.insert(pair, (lower, upper));
    bounds.set_lower(pair.0, pair.1, lower);
    bounds.set_upper(pair.0, pair.1, upper.max(lower + 1.0e-12));
    Ok(())
}

fn midpoint(bounds: &BoundsMatrix, atom1: usize, atom2: usize) -> f64 {
    (bounds.lower(atom1, atom2) + bounds.upper(atom1, atom2)) * 0.5
}

fn distance13(r1: f64, r2: f64, angle: f64) -> f64 {
    (r1 * r1 + r2 * r2 - 2.0 * r1 * r2 * angle.cos())
        .max(0.0)
        .sqrt()
}

fn inferred_angle(bounds: &BoundsMatrix, atom1: usize, center: usize, atom2: usize) -> f64 {
    let r1 = midpoint(bounds, atom1, center);
    let r2 = midpoint(bounds, center, atom2);
    let opposite = midpoint(bounds, atom1, atom2);
    ((r1 * r1 + r2 * r2 - opposite * opposite) / (2.0 * r1 * r2))
        .clamp(-1.0, 1.0)
        .acos()
}

fn torsion_candidates(lower: f64, upper: f64) -> Vec<f64> {
    let mut values = vec![lower, upper, (lower + upper) * 0.5];
    let pi = core::f64::consts::PI;
    let first = (lower / pi).ceil() as i32;
    let last = (upper / pi).floor() as i32;
    for multiple in first..=last {
        let critical = f64::from(multiple) * pi;
        if lower <= critical && critical <= upper {
            values.push(critical);
        }
    }
    values
}

fn distance14(r1: f64, r2: f64, r3: f64, angle1: f64, angle2: f64, torsion: f64) -> f64 {
    let x1 = r1 * angle1.cos();
    let y1 = r1 * angle1.sin();
    let x4 = r2 - r3 * angle2.cos();
    let y4 = r3 * angle2.sin() * torsion.cos();
    let z4 = r3 * angle2.sin() * torsion.sin();
    ((x1 - x4).powi(2) + (y1 - y4).powi(2) + z4.powi(2)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rdkit_embed::{
        RdkitEmbedOptions, rdkit_bounds_matrix, rdkit_bounds_matrix_with_constraints,
        rdkit_embed_molecule_with_constraints,
    };

    fn molecule(smiles: &str) -> Molecule {
        chematic_smiles::parse(smiles).expect("valid test SMILES")
    }

    #[test]
    fn empty_constraints_preserve_the_rdkit_compatible_matrix() {
        let mol = molecule("CCCC");
        let ordinary = rdkit_bounds_matrix(&mol, true, true, false).unwrap();
        let constrained = rdkit_bounds_matrix_with_constraints(
            &mol,
            &InternalCoordinateConstraints::default(),
            true,
            true,
            false,
        )
        .unwrap();
        assert_eq!(ordinary, constrained);

        let options = RdkitEmbedOptions::default();
        let ordinary_coords = crate::rdkit_embed::rdkit_embed_molecule(&mol, &options).unwrap();
        let constrained_coords = rdkit_embed_molecule_with_constraints(
            &mol,
            &options,
            &InternalCoordinateConstraints::default(),
        )
        .unwrap();
        assert_eq!(ordinary_coords, constrained_coords);
    }

    #[test]
    fn custom_bond_angle_and_torsion_enter_the_smoothed_matrix() {
        let mol = molecule("CCCC");
        let constraints = InternalCoordinateConstraints {
            bond_lengths: vec![BondLengthConstraint::around(0, 1, 1.80, 0.01)],
            bond_angles: vec![BondAngleConstraint::around_degrees(0, 1, 2, 100.0, 1.0)],
            torsions: vec![TorsionAngleConstraint::around_degrees(
                0, 1, 2, 3, 60.0, 2.0,
            )],
        };
        let matrix =
            rdkit_bounds_matrix_with_constraints(&mol, &constraints, true, true, false).unwrap();
        assert!((1.79..=1.81).contains(&matrix[1][0]));
        assert!((1.79..=1.81).contains(&matrix[0][1]));
        assert!(matrix[2][0] <= matrix[0][2]);
        assert!(matrix[3][0] <= matrix[0][3]);

        let coordinates = rdkit_embed_molecule_with_constraints(
            &mol,
            &RdkitEmbedOptions::default(),
            &constraints,
        )
        .unwrap();
        let distance = coordinates[0]
            .iter()
            .zip(coordinates[1])
            .map(|(left, right)| (left - right).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!((1.75..=1.85).contains(&distance), "distance={distance}");
    }

    #[test]
    fn inconsistent_constraints_fail_closed_during_smoothing() {
        let mol = molecule("C1CC1");
        let constraints = InternalCoordinateConstraints {
            bond_lengths: vec![
                BondLengthConstraint::around(0, 1, 1.0, 0.001),
                BondLengthConstraint::around(1, 2, 1.0, 0.001),
                BondLengthConstraint::around(2, 0, 3.0, 0.001),
            ],
            ..Default::default()
        };
        let error = rdkit_bounds_matrix_with_constraints(&mol, &constraints, true, true, false)
            .unwrap_err();
        assert!(error.to_string().contains("bounds inconsistent"));
    }

    #[test]
    fn invalid_indices_and_missing_bonds_are_typed_errors() {
        let mol = molecule("CCC");
        let out_of_range = InternalCoordinateConstraints {
            bond_lengths: vec![BondLengthConstraint::around(0, 8, 1.5, 0.1)],
            ..Default::default()
        };
        assert_eq!(
            out_of_range.validate(&mol),
            Err(InternalCoordinateConstraintError::AtomOutOfRange {
                atom: 8,
                atom_count: 3,
            })
        );

        let missing_bond = InternalCoordinateConstraints {
            bond_lengths: vec![BondLengthConstraint::around(0, 2, 1.5, 0.1)],
            ..Default::default()
        };
        assert_eq!(
            missing_bond.validate(&mol),
            Err(InternalCoordinateConstraintError::MissingBond { atom1: 0, atom2: 2 })
        );
    }

    #[test]
    fn fragment_constraints_are_remapped_skipped_or_rejected_explicitly() {
        let constraints = InternalCoordinateConstraints {
            bond_lengths: vec![BondLengthConstraint::around(2, 4, 1.5, 0.1)],
            bond_angles: vec![BondAngleConstraint::around_degrees(2, 4, 6, 109.5, 1.0)],
            torsions: vec![TorsionAngleConstraint::around_degrees(
                2, 4, 6, 8, 180.0, 2.0,
            )],
        };
        let remapped = constraints.remap_fragment(&[8, 6, 4, 2]).unwrap();
        assert_eq!(remapped.bond_lengths[0].atom1, 3);
        assert_eq!(remapped.bond_lengths[0].atom2, 2);
        assert_eq!(remapped.bond_angles[0].atom1, 3);
        assert_eq!(remapped.bond_angles[0].center, 2);
        assert_eq!(remapped.bond_angles[0].atom2, 1);
        assert_eq!(remapped.torsions[0].atom1, 3);
        assert_eq!(remapped.torsions[0].atom2, 2);
        assert_eq!(remapped.torsions[0].atom3, 1);
        assert_eq!(remapped.torsions[0].atom4, 0);

        let skipped = constraints.remap_fragment(&[10, 11]).unwrap();
        assert_eq!(skipped, InternalCoordinateConstraints::default());

        for (fragment, kind) in [
            (vec![2], "bond-length"),
            (vec![2, 4], "bond-angle"),
            (vec![2, 4, 6], "torsion"),
        ] {
            let error = constraints.remap_fragment(&fragment).unwrap_err();
            assert_eq!(
                error,
                InternalCoordinateConstraintError::CrossFragment { kind }
            );
        }
    }

    #[test]
    fn validation_rejects_malformed_constraint_values() {
        let mol = molecule("CCCC");
        let cases = [
            InternalCoordinateConstraints {
                bond_lengths: vec![BondLengthConstraint::around(0, 0, 1.5, 0.1)],
                ..Default::default()
            },
            InternalCoordinateConstraints {
                bond_lengths: vec![BondLengthConstraint {
                    atom1: 0,
                    atom2: 1,
                    lower: f64::NAN,
                    upper: 1.6,
                }],
                ..Default::default()
            },
            InternalCoordinateConstraints {
                bond_angles: vec![BondAngleConstraint {
                    atom1: 0,
                    center: 1,
                    atom2: 2,
                    lower: 0.0,
                    upper: core::f64::consts::PI + 0.1,
                }],
                ..Default::default()
            },
            InternalCoordinateConstraints {
                torsions: vec![TorsionAngleConstraint {
                    atom1: 0,
                    atom2: 1,
                    atom3: 2,
                    atom4: 3,
                    lower: 0.0,
                    upper: 2.0 * core::f64::consts::PI + 0.1,
                }],
                ..Default::default()
            },
        ];
        assert!(matches!(
            cases[0].validate(&mol),
            Err(InternalCoordinateConstraintError::RepeatedAtom { .. })
        ));
        assert!(matches!(
            cases[1].validate(&mol),
            Err(InternalCoordinateConstraintError::NonFinite { .. })
        ));
        assert!(matches!(
            cases[2].validate(&mol),
            Err(InternalCoordinateConstraintError::InvalidInterval {
                kind: "bond-angle",
                ..
            })
        ));
        assert!(matches!(
            cases[3].validate(&mol),
            Err(InternalCoordinateConstraintError::InvalidInterval {
                kind: "torsion",
                ..
            })
        ));
    }

    #[test]
    fn conflicting_pair_constraints_and_error_messages_are_stable() {
        let constraints = InternalCoordinateConstraints {
            bond_lengths: vec![
                BondLengthConstraint::around(0, 1, 1.0, 0.01),
                BondLengthConstraint::around(1, 0, 2.0, 0.01),
            ],
            ..Default::default()
        };
        let mut bounds = BoundsMatrix::new(2);
        assert_eq!(
            constraints.apply(&mut bounds),
            Err(InternalCoordinateConstraintError::ConflictingPair { atom1: 0, atom2: 1 })
        );

        let errors = [
            InternalCoordinateConstraintError::AtomOutOfRange {
                atom: 4,
                atom_count: 2,
            },
            InternalCoordinateConstraintError::RepeatedAtom { kind: "torsion" },
            InternalCoordinateConstraintError::MissingBond { atom1: 0, atom2: 2 },
            InternalCoordinateConstraintError::NonFinite { kind: "bond-angle" },
            InternalCoordinateConstraintError::InvalidInterval {
                kind: "bond-length",
                lower: 2.0,
                upper: 1.0,
            },
            InternalCoordinateConstraintError::ConflictingPair { atom1: 0, atom2: 1 },
            InternalCoordinateConstraintError::CrossFragment { kind: "torsion" },
            InternalCoordinateConstraintError::TriangleSmoothingFailed,
        ];
        let messages: Vec<String> = errors.iter().map(ToString::to_string).collect();
        assert!(messages.iter().all(|message| !message.is_empty()));
        assert!(messages[0].contains("outside molecule"));
        assert!(messages[7].contains("distance bounds inconsistent"));
    }

    #[test]
    fn torsion_intervals_can_cross_the_pi_seam() {
        let mol = molecule("CCCC");
        let constraints = InternalCoordinateConstraints {
            torsions: vec![TorsionAngleConstraint::around_degrees(
                0, 1, 2, 3, 180.0, 5.0,
            )],
            ..Default::default()
        };
        let matrix =
            rdkit_bounds_matrix_with_constraints(&mol, &constraints, true, true, false).unwrap();
        assert!(matrix[3][0].is_finite());
        assert!(matrix[0][3].is_finite());
        assert!(matrix[3][0] <= matrix[0][3]);
    }
}
