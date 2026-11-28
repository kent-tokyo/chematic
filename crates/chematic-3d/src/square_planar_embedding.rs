//! Experimental square-planar angle constraints for a single ETKDG attempt.

use crate::rdkit_embed::{
    BondAngleConstraint, InternalCoordinateConstraints, RdkitEmbedError, RdkitEmbedOptions,
    rdkit_embed_molecule_with_constraints,
};
use crate::square_planar_geometry::{declared_neighbors, square_planar_angles};
use chematic_core::{Chirality, Coords3D, Molecule, Point3, STEREO_VACANT_SENTINEL};

/// Geometric policy for declared trans and cis ligand pairs, in degrees.
/// Defaults are an experimental screen, not chemically universal tolerances.
#[derive(Clone, Copy, Debug)]
pub struct SquarePlanarEmbedConfig {
    pub trans_min_degrees: f64,
    pub cis_tolerance_degrees: f64,
}

impl Default for SquarePlanarEmbedConfig {
    fn default() -> Self {
        Self {
            trans_min_degrees: 150.0,
            cis_tolerance_degrees: 30.0,
        }
    }
}

/// Failure to prepare, embed, or meet the chosen angle policy.
#[derive(Debug)]
pub enum SquarePlanarEmbedError {
    InvalidThreshold,
    NoDeclaration,
    UnsupportedCenter {
        center: u32,
        reason: &'static str,
    },
    ValidationUnavailable {
        center: u32,
        reason: &'static str,
    },
    Embed(RdkitEmbedError),
    InvalidCoordinates,
    AngleOutOfRange {
        center: u32,
        atoms: [u32; 2],
        degrees: f64,
        trans: bool,
    },
}

impl core::fmt::Display for SquarePlanarEmbedError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidThreshold => {
                write!(f, "invalid or overlapping square-planar angle thresholds")
            }
            Self::NoDeclaration => write!(f, "molecule has no square-planar declaration"),
            Self::UnsupportedCenter { center, reason } => {
                write!(f, "square-planar center {center}: {reason}")
            }
            Self::ValidationUnavailable { center, reason } => {
                write!(
                    f,
                    "square-planar center {center}: output validation unavailable: {reason}"
                )
            }
            Self::Embed(error) => write!(f, "square-planar embedding failed: {error}"),
            Self::InvalidCoordinates => {
                write!(f, "square-planar embedding returned invalid coordinates")
            }
            Self::AngleOutOfRange {
                center,
                atoms,
                degrees,
                trans,
            } => write!(
                f,
                "square-planar center {center}: {} pair {}-{} has angle {degrees:.3} degrees outside the selected limits",
                if *trans { "trans" } else { "cis" },
                atoms[0],
                atoms[1]
            ),
        }
    }
}
impl std::error::Error for SquarePlanarEmbedError {}

/// Generate once using explicit square-planar angle bounds, then remeasure
/// the returned coordinates. A failure returns no coordinates. Only real
/// ligand pairs are constrained. Other stereo and geometry are not validated.
pub fn rdkit_embed_square_planar(
    mol: &Molecule,
    embed_options: &RdkitEmbedOptions,
    config: SquarePlanarEmbedConfig,
) -> Result<Vec<[f64; 3]>, SquarePlanarEmbedError> {
    embed_with_bounds(mol, embed_options, config, config)
}

fn embed_with_bounds(
    mol: &Molecule,
    embed_options: &RdkitEmbedOptions,
    bounds: SquarePlanarEmbedConfig,
    acceptance: SquarePlanarEmbedConfig,
) -> Result<Vec<[f64; 3]>, SquarePlanarEmbedError> {
    let constraints = square_planar_constraints(mol, bounds)?;
    let coords = rdkit_embed_molecule_with_constraints(mol, embed_options, &constraints)
        .map_err(SquarePlanarEmbedError::Embed)?;
    validate_angles(mol, &coords, acceptance)?;
    Ok(coords)
}

/// A successful bounded retry. Coordinates meet the original angle policy.
#[derive(Debug)]
pub struct SquarePlanarRetryResult {
    pub coords: Vec<[f64; 3]>,
    pub attempts: u8,
    pub used_tightened_bounds: bool,
    /// Reason the first set of bounds was rejected, if a retry succeeded.
    pub initial_rejection: Option<String>,
}

/// The original failure, or both failures when a second attempt was made.
#[derive(Debug)]
pub enum SquarePlanarRetryError {
    Initial(SquarePlanarEmbedError),
    Both {
        initial: Box<SquarePlanarEmbedError>,
        retry: Box<SquarePlanarEmbedError>,
    },
}

impl core::fmt::Display for SquarePlanarRetryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Initial(error) => write!(f, "{error}"),
            Self::Both { initial, retry } => write!(
                f,
                "initial attempt: {initial}; tightened-bounds retry: {retry}"
            ),
        }
    }
}
impl std::error::Error for SquarePlanarRetryError {}

/// Run one attempt and, only after a final angle rejection, at most one more
/// with the same seed and tighter bounds. The original acceptance policy is
/// applied to both outputs. Setup and embedding errors are not retried.
pub fn rdkit_embed_square_planar_with_retry(
    mol: &Molecule,
    embed_options: &RdkitEmbedOptions,
    config: SquarePlanarEmbedConfig,
) -> Result<SquarePlanarRetryResult, SquarePlanarRetryError> {
    let initial = match rdkit_embed_square_planar(mol, embed_options, config) {
        Ok(coords) => {
            return Ok(SquarePlanarRetryResult {
                coords,
                attempts: 1,
                used_tightened_bounds: false,
                initial_rejection: None,
            });
        }
        Err(error @ SquarePlanarEmbedError::AngleOutOfRange { .. }) => error,
        Err(error) => return Err(SquarePlanarRetryError::Initial(error)),
    };
    let tighter = tightened_bounds(config);
    if tighter.trans_min_degrees == config.trans_min_degrees
        && tighter.cis_tolerance_degrees == config.cis_tolerance_degrees
    {
        return Err(SquarePlanarRetryError::Initial(initial));
    }
    match embed_with_bounds(mol, embed_options, tighter, config) {
        Ok(coords) => Ok(SquarePlanarRetryResult {
            coords,
            attempts: 2,
            used_tightened_bounds: true,
            initial_rejection: Some(initial.to_string()),
        }),
        Err(retry) => Err(SquarePlanarRetryError::Both {
            initial: Box::new(initial),
            retry: Box::new(retry),
        }),
    }
}

fn tightened_bounds(config: SquarePlanarEmbedConfig) -> SquarePlanarEmbedConfig {
    SquarePlanarEmbedConfig {
        trans_min_degrees: (config.trans_min_degrees + 20.0).min(180.0),
        cis_tolerance_degrees: (config.cis_tolerance_degrees - 10.0).max(0.0),
    }
}

fn square_planar_constraints(
    mol: &Molecule,
    config: SquarePlanarEmbedConfig,
) -> Result<InternalCoordinateConstraints, SquarePlanarEmbedError> {
    if !config.trans_min_degrees.is_finite()
        || !config.cis_tolerance_degrees.is_finite()
        || !(0.0..=180.0).contains(&config.trans_min_degrees)
        || !(0.0..90.0).contains(&config.cis_tolerance_degrees)
        || config.trans_min_degrees <= 90.0 + config.cis_tolerance_degrees
    {
        return Err(SquarePlanarEmbedError::InvalidThreshold);
    }
    let mut constraints = InternalCoordinateConstraints::default();
    let mut centers = 0;
    for (center, atom) in mol.atoms() {
        let permutation = match atom.chirality {
            Chirality::SquarePlanar(p) => p,
            Chirality::SquarePlanarUnnumbered => {
                return Err(SquarePlanarEmbedError::UnsupportedCenter {
                    center: center.0,
                    reason: "unnumbered_permutation",
                });
            }
            _ => continue,
        };
        centers += 1;
        let neighbors = declared_neighbors(mol, center).map_err(|reason| {
            SquarePlanarEmbedError::UnsupportedCenter {
                center: center.0,
                reason,
            }
        })?;
        let trans = permutation.trans_pairs();
        for i in 0..4 {
            for j in i + 1..4 {
                if neighbors[i] == STEREO_VACANT_SENTINEL || neighbors[j] == STEREO_VACANT_SENTINEL
                {
                    continue;
                }
                let (lower, upper) = if trans.contains(&(i as u8, j as u8)) {
                    (config.trans_min_degrees, 180.0)
                } else {
                    (
                        90.0 - config.cis_tolerance_degrees,
                        90.0 + config.cis_tolerance_degrees,
                    )
                };
                constraints.bond_angles.push(BondAngleConstraint {
                    atom1: neighbors[i] as usize,
                    center: center.0 as usize,
                    atom2: neighbors[j] as usize,
                    lower: lower.to_radians(),
                    upper: upper.to_radians(),
                });
            }
        }
    }
    if centers == 0 {
        return Err(SquarePlanarEmbedError::NoDeclaration);
    }
    Ok(constraints)
}

fn validate_angles(
    mol: &Molecule,
    coords: &[[f64; 3]],
    config: SquarePlanarEmbedConfig,
) -> Result<(), SquarePlanarEmbedError> {
    let geometry = Coords3D {
        points: coords
            .iter()
            .map(|&[x, y, z]| Point3::new(x, y, z))
            .collect(),
    };
    let reports =
        square_planar_angles(mol, &geometry).ok_or(SquarePlanarEmbedError::InvalidCoordinates)?;
    for report in reports {
        if let Some(reason) = report.unavailable {
            return Err(SquarePlanarEmbedError::ValidationUnavailable {
                center: report.center.0,
                reason,
            });
        }
        for (trans, angles) in [(true, report.trans), (false, report.cis)] {
            for angle in angles {
                let valid = if trans {
                    angle.degrees >= config.trans_min_degrees
                } else {
                    (angle.degrees - 90.0).abs() <= config.cis_tolerance_degrees
                };
                if !valid {
                    return Err(SquarePlanarEmbedError::AngleOutOfRange {
                        center: report.center.0,
                        atoms: angle.atoms,
                        degrees: angle.degrees,
                        trans,
                    });
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::square_planar_geometry::square_planar_angles;

    #[test]
    fn constructs_declared_pair_bounds_with_and_without_vacancy() {
        for (smiles, total) in [("[Pt@SP1](F)(Cl)(Br)I", 6), ("F[Pt@SP3](Cl)Br", 3)] {
            let mol = chematic_smiles::parse(smiles).unwrap();
            let constraints =
                square_planar_constraints(&mol, SquarePlanarEmbedConfig::default()).unwrap();
            assert_eq!(constraints.bond_angles.len(), total);
            assert_eq!(
                constraints
                    .bond_angles
                    .iter()
                    .filter(|a| (a.upper.to_degrees() - 180.0).abs() < 1e-10)
                    .count(),
                if total == 6 { 2 } else { 1 }
            );
        }
    }

    #[test]
    fn final_angle_check_does_not_accept_a_distorted_result() {
        let mol = chematic_smiles::parse("F[Pt@SP3](Cl)Br").unwrap();
        let ideal = [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
        ];
        let config = SquarePlanarEmbedConfig::default();
        assert!(validate_angles(&mol, &ideal, config).is_ok());
        let distorted = [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        assert!(matches!(
            validate_angles(&mol, &distorted, config),
            Err(SquarePlanarEmbedError::AngleOutOfRange { trans: true, .. })
        ));
        let bad_cis = [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
        ];
        assert!(matches!(
            validate_angles(&mol, &bad_cis, config),
            Err(SquarePlanarEmbedError::AngleOutOfRange { trans: false, .. })
        ));
        let zero = [[0.0; 3]; 4];
        assert!(matches!(
            validate_angles(&mol, &zero, config),
            Err(SquarePlanarEmbedError::ValidationUnavailable {
                reason: "degenerate_ligand_vector",
                ..
            })
        ));
        assert!(
            square_planar_angles(
                &mol,
                &Coords3D {
                    points: ideal
                        .iter()
                        .map(|p| Point3::new(p[0], p[1], p[2]))
                        .collect()
                }
            )
            .is_some()
        );
    }

    #[test]
    fn rejects_missing_unsupported_and_ambiguous_input() {
        let simple = chematic_smiles::parse("CCO").unwrap();
        assert!(matches!(
            square_planar_constraints(&simple, SquarePlanarEmbedConfig::default()),
            Err(SquarePlanarEmbedError::NoDeclaration)
        ));
        let partial = chematic_smiles::parse("F[Pt@SP3]Cl").unwrap();
        assert!(matches!(
            square_planar_constraints(&partial, SquarePlanarEmbedConfig::default()),
            Err(SquarePlanarEmbedError::UnsupportedCenter { .. })
        ));
        let unnumbered = chematic_smiles::parse("[Pt@SP](F)(Cl)(Br)I").unwrap();
        assert!(matches!(
            square_planar_constraints(&unnumbered, SquarePlanarEmbedConfig::default()),
            Err(SquarePlanarEmbedError::UnsupportedCenter { .. })
        ));
        let mol = chematic_smiles::parse("F[Pt@SP3](Cl)Br").unwrap();
        for config in [
            SquarePlanarEmbedConfig {
                trans_min_degrees: 120.0,
                cis_tolerance_degrees: 30.0,
            },
            SquarePlanarEmbedConfig {
                trans_min_degrees: f64::NAN,
                cis_tolerance_degrees: 30.0,
            },
            SquarePlanarEmbedConfig {
                trans_min_degrees: 150.0,
                cis_tolerance_degrees: -1.0,
            },
        ] {
            assert!(matches!(
                square_planar_constraints(&mol, config),
                Err(SquarePlanarEmbedError::InvalidThreshold)
            ));
        }
    }

    #[test]
    fn retry_tightening_is_bounded_and_initial_errors_are_preserved() {
        let tighter = tightened_bounds(SquarePlanarEmbedConfig::default());
        assert_eq!(
            (tighter.trans_min_degrees, tighter.cis_tolerance_degrees),
            (170.0, 20.0)
        );
        let clamped = tightened_bounds(SquarePlanarEmbedConfig {
            trans_min_degrees: 179.0,
            cis_tolerance_degrees: 4.0,
        });
        assert_eq!(
            (clamped.trans_min_degrees, clamped.cis_tolerance_degrees),
            (180.0, 0.0)
        );
        let mol = chematic_smiles::parse("CCO").unwrap();
        let err = rdkit_embed_square_planar_with_retry(
            &mol,
            &RdkitEmbedOptions::default(),
            SquarePlanarEmbedConfig::default(),
        )
        .unwrap_err();
        assert!(matches!(
            err,
            SquarePlanarRetryError::Initial(SquarePlanarEmbedError::NoDeclaration)
        ));
    }

    #[test]
    fn multiple_fragments_keep_each_centers_constraints_local() {
        use chematic_core::AtomIdx;
        let mol = chematic_smiles::parse("[Pt@SP1](F)(Cl)(Br)I.F[Pt@SP3](Cl)Br").unwrap();
        for config in [
            SquarePlanarEmbedConfig::default(),
            tightened_bounds(SquarePlanarEmbedConfig::default()),
        ] {
            let constraints = square_planar_constraints(&mol, config).unwrap();
            assert_eq!(constraints.bond_angles.len(), 9);
            for angle in constraints.bond_angles {
                let center = AtomIdx(angle.center as u32);
                assert!(
                    mol.bond_between(center, AtomIdx(angle.atom1 as u32))
                        .is_some()
                );
                assert!(
                    mol.bond_between(center, AtomIdx(angle.atom2 as u32))
                        .is_some()
                );
            }
        }
    }
}
