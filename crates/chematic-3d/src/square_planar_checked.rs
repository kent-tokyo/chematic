//! Optional square-planar shape and other-declared-stereo screening.

use crate::rdkit_embed::RdkitEmbedOptions;
use crate::square_planar_embedding::{
    SquarePlanarEmbedConfig, SquarePlanarRetryError, SquarePlanarRetryResult,
    rdkit_embed_square_planar_with_retry,
};
use crate::square_planar_geometry::square_planar_angles;
use crate::stereo_constraints::{
    StereoElement, StereoRejectionReason, StereoStatus, verify_stereo,
};
use chematic_core::{Chirality, Coords3D, Molecule, Point3};

/// Additional, explicitly selected planarity threshold. This dimensionless
/// limit is not a universal physical definition of square-planar geometry.
#[derive(Clone, Copy, Debug)]
pub struct SquarePlanarCheckedConfig {
    pub angles: SquarePlanarEmbedConfig,
    pub max_unit_triple_product: f64,
}

#[derive(Debug)]
pub enum SquarePlanarCheckedError {
    InvalidPlanarityLimit,
    Embedding(SquarePlanarRetryError),
    InvalidCoordinates,
    PlanarityUnavailable {
        center: u32,
    },
    PlanarityExceeded {
        center: u32,
        measured: f64,
        limit: f64,
    },
    OtherStereo {
        violated: Vec<StereoElement>,
        unevaluable: Vec<(StereoElement, StereoRejectionReason)>,
    },
}

impl core::fmt::Display for SquarePlanarCheckedError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidPlanarityLimit => write!(
                f,
                "max_unit_triple_product must be finite and between 0 and 1"
            ),
            Self::Embedding(error) => write!(f, "{error}"),
            Self::InvalidCoordinates => write!(f, "invalid generated coordinates"),
            Self::PlanarityUnavailable { center } => {
                write!(f, "square-planar center {center}: planarity unavailable")
            }
            Self::PlanarityExceeded {
                center,
                measured,
                limit,
            } => write!(
                f,
                "square-planar center {center}: unit triple product {measured:.6} exceeds {limit:.6}"
            ),
            Self::OtherStereo {
                violated,
                unevaluable,
            } => write!(
                f,
                "other declared stereo rejected: violated={violated:?}; unevaluable={unevaluable:?}"
            ),
        }
    }
}
impl std::error::Error for SquarePlanarCheckedError {}

#[derive(Debug)]
pub struct SquarePlanarCheckedResult {
    pub embedding: SquarePlanarRetryResult,
    pub max_unit_triple_product: f64,
    pub other_stereo_checked: usize,
}

/// Generate with the bounded square-planar retry, then apply the explicit
/// planarity limit and the library's supported non-square-planar stereo checks.
/// Unknown stereo is rejected. Bond lengths, energy and unimplemented stereo
/// classes outside the verifier's stated coverage are not certified.
pub fn rdkit_embed_square_planar_checked(
    mol: &Molecule,
    options: &RdkitEmbedOptions,
    config: SquarePlanarCheckedConfig,
) -> Result<SquarePlanarCheckedResult, SquarePlanarCheckedError> {
    if !config.max_unit_triple_product.is_finite()
        || !(0.0..=1.0).contains(&config.max_unit_triple_product)
    {
        return Err(SquarePlanarCheckedError::InvalidPlanarityLimit);
    }
    let embedding = rdkit_embed_square_planar_with_retry(mol, options, config.angles)
        .map_err(SquarePlanarCheckedError::Embedding)?;
    let coords = Coords3D {
        points: embedding
            .coords
            .iter()
            .map(|&[x, y, z]| Point3::new(x, y, z))
            .collect(),
    };
    let max = check_planarity(mol, &coords, config.max_unit_triple_product)?;
    let other_stereo_checked = check_other_stereo(mol, &coords)?;
    Ok(SquarePlanarCheckedResult {
        embedding,
        max_unit_triple_product: max,
        other_stereo_checked,
    })
}

fn check_planarity(
    mol: &Molecule,
    coords: &Coords3D,
    limit: f64,
) -> Result<f64, SquarePlanarCheckedError> {
    let reports =
        square_planar_angles(mol, coords).ok_or(SquarePlanarCheckedError::InvalidCoordinates)?;
    let mut maximum = 0.0_f64;
    for report in reports {
        let value = report.max_unit_triple_product.ok_or(
            SquarePlanarCheckedError::PlanarityUnavailable {
                center: report.center.0,
            },
        )?;
        if value > limit {
            return Err(SquarePlanarCheckedError::PlanarityExceeded {
                center: report.center.0,
                measured: value,
                limit,
            });
        }
        maximum = maximum.max(value);
    }
    Ok(maximum)
}

fn check_other_stereo(
    mol: &Molecule,
    coords: &Coords3D,
) -> Result<usize, SquarePlanarCheckedError> {
    let report = verify_stereo(mol, coords);
    let mut checked = 0;
    let mut violated = Vec::new();
    let mut unevaluable = Vec::new();
    for item in report.tetrahedral {
        if matches!(mol.atom(item.atom).chirality, Chirality::SquarePlanar(_))
            && matches!(
                item.status,
                StereoStatus::Unevaluable(StereoRejectionReason::UnsupportedCoordination)
            )
        {
            continue;
        }
        checked += 1;
        match item.status {
            StereoStatus::Satisfied => (),
            StereoStatus::Violated => violated.push(StereoElement::Tetrahedral(item.atom)),
            StereoStatus::Unevaluable(reason) => {
                unevaluable.push((StereoElement::Tetrahedral(item.atom), reason))
            }
        }
    }
    for item in report.double_bond {
        checked += 1;
        match item.status {
            StereoStatus::Satisfied => (),
            StereoStatus::Violated => violated.push(StereoElement::DoubleBond(item.bond)),
            StereoStatus::Unevaluable(reason) => {
                unevaluable.push((StereoElement::DoubleBond(item.bond), reason))
            }
        }
    }
    if !violated.is_empty() || !unevaluable.is_empty() {
        return Err(SquarePlanarCheckedError::OtherStereo {
            violated,
            unevaluable,
        });
    }
    Ok(checked)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rdkit_embed::{RdkitEmbedOptions, etkdg_embed_multiple_confs_corrected};

    fn geometry(values: &[[f64; 3]]) -> Coords3D {
        Coords3D {
            points: values
                .iter()
                .map(|&[x, y, z]| Point3::new(x, y, z))
                .collect(),
        }
    }

    #[test]
    fn planarity_rejects_nonplanar_geometry_that_passes_angle_screen() {
        let mol = chematic_smiles::parse("[Pt@SP1](F)(Cl)(Br)I").unwrap();
        let theta = 20.0_f64.to_radians();
        let coords = geometry(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-theta.cos(), 0.0, theta.sin()],
            [0.0, -theta.cos(), theta.sin()],
        ]);
        let report = square_planar_angles(&mol, &coords).unwrap().remove(0);
        assert!(report.trans.iter().all(|a| a.degrees >= 150.0));
        assert!(
            report
                .cis
                .iter()
                .all(|a| (60.0..=120.0).contains(&a.degrees))
        );
        assert!(matches!(
            check_planarity(&mol, &coords, 0.2),
            Err(SquarePlanarCheckedError::PlanarityExceeded { center: 0, .. })
        ));
        assert!((check_planarity(&mol, &coords, 0.4).unwrap() - theta.sin()).abs() < 1e-12);
    }

    #[test]
    fn other_tetrahedral_stereo_rejects_mirror_and_unevaluable_coordinates() {
        let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse("C[C@H](O)F").unwrap());
        let coords = etkdg_embed_multiple_confs_corrected(&mol, 1, &RdkitEmbedOptions::default())
            .unwrap()
            .remove(0);
        assert_eq!(check_other_stereo(&mol, &geometry(&coords)).unwrap(), 1);
        let mirror: Vec<_> = coords.iter().map(|&[x, y, z]| [-x, y, z]).collect();
        assert!(matches!(check_other_stereo(&mol, &geometry(&mirror)),
            Err(SquarePlanarCheckedError::OtherStereo { violated, .. }) if violated.len() == 1));
        let zero = Coords3D::new_zeroed(mol.atom_count());
        assert!(matches!(check_other_stereo(&mol, &zero),
            Err(SquarePlanarCheckedError::OtherStereo { unevaluable, .. }) if unevaluable.len() == 1));
    }

    #[test]
    fn other_double_bond_stereo_is_checked() {
        let mol = chematic_smiles::parse("F/C=C/F").unwrap();
        let trans = geometry(&[[0., 1., 0.], [0., 0., 0.], [1., 0., 0.], [1., -1., 0.]]);
        let cis = geometry(&[[0., 1., 0.], [0., 0., 0.], [1., 0., 0.], [1., 1., 0.]]);
        assert_eq!(check_other_stereo(&mol, &trans).unwrap(), 1);
        assert!(matches!(check_other_stereo(&mol, &cis),
            Err(SquarePlanarCheckedError::OtherStereo { violated, .. }) if violated.len() == 1));
    }

    #[test]
    fn only_numbered_square_planar_centers_are_exempted_from_other_stereo() {
        let numbered = chematic_smiles::parse("[Pt@SP1](F)(Cl)(Br)I").unwrap();
        let coords = geometry(&[
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [-1., 0., 0.],
            [0., -1., 0.],
        ]);
        assert_eq!(check_other_stereo(&numbered, &coords).unwrap(), 0);
        let unnumbered = chematic_smiles::parse("[Pt@SP](F)(Cl)(Br)I").unwrap();
        assert!(matches!(check_other_stereo(&unnumbered, &coords),
            Err(SquarePlanarCheckedError::OtherStereo { unevaluable, .. }) if unevaluable.len() == 1));
        let bad_limit = SquarePlanarCheckedConfig {
            angles: SquarePlanarEmbedConfig::default(),
            max_unit_triple_product: f64::NAN,
        };
        assert!(matches!(
            rdkit_embed_square_planar_checked(&numbered, &RdkitEmbedOptions::default(), bad_limit),
            Err(SquarePlanarCheckedError::InvalidPlanarityLimit)
        ));
    }

    #[test]
    fn unsupported_other_selenium_stereo_is_not_silently_accepted() {
        use chematic_core::AtomIdx;
        let smiles = "Brc1ccc2c3c(cccc13)[Se@@](c1ccccc1)=[Se@SP3]2Br";
        let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse(smiles).unwrap());
        let coords = Coords3D::new_zeroed(mol.atom_count());
        let error = check_other_stereo(&mol, &coords).unwrap_err();
        let SquarePlanarCheckedError::OtherStereo { unevaluable, .. } = error else {
            panic!("expected unevaluable selenium stereo")
        };
        assert!(unevaluable.contains(&(
            StereoElement::Tetrahedral(AtomIdx(11)),
            StereoRejectionReason::UnsupportedCoordination,
        )));
        assert!(
            !unevaluable
                .iter()
                .any(|(element, _)| *element == StereoElement::Tetrahedral(AtomIdx(18)))
        );
    }
}
