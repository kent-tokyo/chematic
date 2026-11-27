//! Measured angles at explicitly declared square-planar centers.
//!
//! This is a diagnostic, not a stereo acceptance test. It neither changes
//! embedding nor replaces the existing conservative stereo verifier.

use chematic_core::{
    AtomIdx, Chirality, Coords3D, Molecule, STEREO_VACANT_SENTINEL, SquarePlanarPermutation,
};

/// One ligand pair and its angle at the center, in degrees.
#[derive(Clone, Debug, PartialEq)]
pub struct LigandAngle {
    pub atoms: [u32; 2],
    pub degrees: f64,
}

/// Declared trans and cis pair angles, in parser-recorded neighbor order.
#[derive(Clone, Debug, PartialEq)]
pub struct SquarePlanarAngles {
    pub center: AtomIdx,
    pub permutation: Option<SquarePlanarPermutation>,
    /// Four declared positions; `None` is a vacant coordination site.
    pub neighbors: Option<[Option<u32>; 4]>,
    pub trans: Vec<LigandAngle>,
    pub cis: Vec<LigandAngle>,
    /// `None` means angles were measured, not that the geometry passed a check.
    pub unavailable: Option<&'static str>,
}

/// Measure ligand-center-ligand angles for numbered SP declarations with four
/// explicit ligands, or three explicit ligands and one recorded vacant site.
///
/// `None` indicates globally invalid dimensions or nonfinite coordinates.
/// Missing/invalid neighbor order, unnumbered SP and coincident or numerically
/// unrepresentable vectors are reported per center without inferring geometry.
pub fn square_planar_angles(mol: &Molecule, coords: &Coords3D) -> Option<Vec<SquarePlanarAngles>> {
    if coords.atom_count() != mol.atom_count() || !coords.is_finite() {
        return None;
    }
    let mut reports = Vec::new();
    for (center, atom) in mol.atoms() {
        let permutation = match atom.chirality {
            Chirality::SquarePlanar(p) => Some(p),
            Chirality::SquarePlanarUnnumbered => None,
            _ => continue,
        };
        reports.push(measure_center(mol, coords, center, permutation));
    }
    Some(reports)
}

fn measure_center(
    mol: &Molecule,
    coords: &Coords3D,
    center: AtomIdx,
    permutation: Option<SquarePlanarPermutation>,
) -> SquarePlanarAngles {
    let mut result = SquarePlanarAngles {
        center,
        permutation,
        neighbors: None,
        trans: Vec::new(),
        cis: Vec::new(),
        unavailable: None,
    };
    let Some(permutation) = permutation else {
        result.unavailable = Some("unnumbered_permutation");
        return result;
    };
    let Some(order) = mol.stereo_neighbor_order(center) else {
        result.unavailable = Some("missing_neighbor_order");
        return result;
    };
    let Ok(neighbors): Result<[u32; 4], _> = order.try_into() else {
        result.unavailable = Some("unsupported_neighbor_order");
        return result;
    };
    let vacancies = neighbors
        .iter()
        .filter(|&&nb| nb == STEREO_VACANT_SENTINEL)
        .count();
    if vacancies > 1
        || mol.degree(center) + vacancies != 4
        || neighbors.iter().enumerate().any(|(i, &nb)| {
            nb != STEREO_VACANT_SENTINEL
                && (nb as usize >= mol.atom_count()
                    || nb == center.0
                    || neighbors[..i].contains(&nb)
                    || mol.bond_between(center, AtomIdx(nb)).is_none())
        })
    {
        result.unavailable = Some("unsupported_neighbor_order");
        return result;
    }
    result.neighbors = Some(neighbors.map(|nb| (nb != STEREO_VACANT_SENTINEL).then_some(nb)));
    let origin = coords.get(center);
    let mut unit = [[0.0; 3]; 4];
    for (i, &nb) in neighbors.iter().enumerate() {
        if nb == STEREO_VACANT_SENTINEL {
            continue;
        }
        let p = coords.get(AtomIdx(nb));
        let v = [p.x - origin.x, p.y - origin.y, p.z - origin.z];
        let norm = v[0].hypot(v[1]).hypot(v[2]);
        if !norm.is_finite() || norm <= 1e-12 {
            result.unavailable = Some("degenerate_ligand_vector");
            return result;
        }
        unit[i] = v.map(|x| x / norm);
    }
    let trans = permutation.trans_pairs();
    for i in 0..4 {
        for j in i + 1..4 {
            if neighbors[i] == STEREO_VACANT_SENTINEL || neighbors[j] == STEREO_VACANT_SENTINEL {
                continue;
            }
            let cosine: f64 = (0..3).map(|k| unit[i][k] * unit[j][k]).sum();
            let angle = LigandAngle {
                atoms: [neighbors[i], neighbors[j]],
                degrees: cosine.clamp(-1.0, 1.0).acos().to_degrees(),
            };
            if trans.contains(&(i as u8, j as u8)) {
                result.trans.push(angle);
            } else {
                result.cis.push(angle);
            }
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::Point3;

    #[test]
    fn all_permutations_measure_declared_trans_and_cis_pairs() {
        for tag in 1..=3 {
            let mol = chematic_smiles::parse(&format!("[Pt@SP{tag}](F)(Cl)(Br)I")).unwrap();
            let center = AtomIdx(0);
            let mut coords = Coords3D::new_zeroed(mol.atom_count());
            let Chirality::SquarePlanar(permutation) = mol.atom(center).chirality else {
                panic!()
            };
            let neighbors = mol.stereo_neighbor_order(center).unwrap();
            let pairs = permutation.trans_pairs();
            for (position, vector) in [
                (pairs[0].0, [1.0, 0.0, 0.0]),
                (pairs[0].1, [-1.0, 0.0, 0.0]),
                (pairs[1].0, [0.0, 2.0, 0.0]),
                (pairs[1].1, [0.0, -2.0, 0.0]),
            ] {
                coords.set(
                    AtomIdx(neighbors[position as usize]),
                    Point3::new(vector[0], vector[1], vector[2]),
                );
            }
            let report = square_planar_angles(&mol, &coords).unwrap().remove(0);
            assert_eq!(report.permutation, Some(permutation));
            assert_eq!(
                report.neighbors.unwrap(),
                [
                    Some(neighbors[0]),
                    Some(neighbors[1]),
                    Some(neighbors[2]),
                    Some(neighbors[3])
                ]
            );
            assert_eq!(report.unavailable, None);
            assert_eq!((report.trans.len(), report.cis.len()), (2, 4));
            assert!(
                report
                    .trans
                    .iter()
                    .all(|a| (a.degrees - 180.0).abs() < 1e-10)
            );
            assert!(report.cis.iter().all(|a| (a.degrees - 90.0).abs() < 1e-10));
            for point in &mut coords.points {
                *point = Point3::new(-point.y + 4.0, point.x - 9.0, point.z + 2.0);
            }
            let rotated = square_planar_angles(&mol, &coords).unwrap().remove(0);
            assert_eq!(report, rotated);
        }
    }

    #[test]
    fn invalid_inputs_and_unavailable_centers_are_explicit() {
        let mol = chematic_smiles::parse("[Pt@SP1](F)(Cl)(Br)I").unwrap();
        assert!(square_planar_angles(&mol, &Coords3D::new_zeroed(1)).is_none());
        let mut coords = Coords3D::new_zeroed(mol.atom_count());
        let report = square_planar_angles(&mol, &coords).unwrap().remove(0);
        assert_eq!(report.unavailable, Some("degenerate_ligand_vector"));
        assert!(report.trans.is_empty());
        coords.points[0].x = -1e308;
        coords.points[1].x = 1e308;
        let report = square_planar_angles(&mol, &coords).unwrap().remove(0);
        assert_eq!(report.unavailable, Some("degenerate_ligand_vector"));
        assert!(report.trans.is_empty());
        coords.points[0].x = f64::NAN;
        assert!(square_planar_angles(&mol, &coords).is_none());
        let unnumbered = chematic_smiles::parse("[Pt@SP](F)(Cl)(Br)I").unwrap();
        let report =
            square_planar_angles(&unnumbered, &Coords3D::new_zeroed(unnumbered.atom_count()))
                .unwrap();
        assert_eq!(report[0].unavailable, Some("unnumbered_permutation"));
        let two_ligands = chematic_smiles::parse("F[Pt@SP3]Cl").unwrap();
        let report = square_planar_angles(
            &two_ligands,
            &Coords3D::new_zeroed(two_ligands.atom_count()),
        )
        .unwrap();
        assert_eq!(report[0].unavailable, Some("unsupported_neighbor_order"));
        let ordinary = chematic_smiles::parse("CCO").unwrap();
        assert!(
            square_planar_angles(&ordinary, &Coords3D::new_zeroed(ordinary.atom_count()))
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn vacant_site_retains_only_measurable_pairs() {
        let mol = chematic_smiles::parse("F[Pt@SP3](Cl)Br").unwrap();
        let mut coords = Coords3D::new_zeroed(mol.atom_count());
        coords.set(AtomIdx(0), Point3::new(1.0, 0.0, 0.0));
        coords.set(AtomIdx(2), Point3::new(0.0, 1.0, 0.0));
        coords.set(AtomIdx(3), Point3::new(-1.0, 0.0, 0.0));
        let report = square_planar_angles(&mol, &coords).unwrap().remove(0);
        assert_eq!(report.neighbors, Some([Some(0), None, Some(2), Some(3)]));
        assert_eq!(report.unavailable, None);
        assert_eq!(report.trans.len(), 1);
        assert_eq!(report.cis.len(), 2);
        assert_eq!(report.trans[0].atoms, [0, 3]);
        assert!((report.trans[0].degrees - 180.0).abs() < 1e-10);
    }
}
