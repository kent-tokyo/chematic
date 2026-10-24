//! Bounded torsion moves around the best converged MMFF94 conformer.

use chematic_core::{AtomIdx, Coords3D, Molecule, Point3};

use crate::mol_transforms::{get_dihedral, try_set_dihedral};
use crate::stereo_constraints::verify_stereo;

use super::{RdkitMmffConformer, optimize_mmff_conformer, select_quality_candidates};

pub(super) fn refine_best_with_torsion_moves(
    mol: &Molecule,
    candidates: Vec<RdkitMmffConformer>,
    max_iters: u32,
) -> Vec<RdkitMmffConformer> {
    if max_iters == 0 {
        return candidates;
    }
    let Some(source) = candidates
        .iter()
        .filter(|candidate| candidate.status == 0 && candidate.energy.is_finite())
        .min_by(|left, right| left.energy.total_cmp(&right.energy))
    else {
        return candidates;
    };
    let source = source.clone();
    let Ok(field) = chematic_ff::rdkit_mmff::RdkitMmffField::new(mol, &source.coords, 100.0, true)
    else {
        return candidates;
    };
    let source_coords = Coords3D {
        points: source
            .coords
            .iter()
            .map(|&[x, y, z]| Point3::new(x, y, z))
            .collect(),
    };
    let mut best = source;
    for (b, c) in chematic_chem::rotatable_bond_atom_pairs(mol)
        .into_iter()
        .take(12)
    {
        let Some(a) = heavy_neighbor(mol, b, c) else {
            continue;
        };
        let Some(d) = heavy_neighbor(mol, c, b) else {
            continue;
        };
        let Some(angle) = get_dihedral(&source_coords, a, b, c, d) else {
            continue;
        };
        for delta in [-std::f64::consts::TAU / 3.0, std::f64::consts::TAU / 3.0] {
            let Ok(rotated) = try_set_dihedral(&source_coords, mol, a, b, c, d, angle + delta)
            else {
                continue;
            };
            let coords: Vec<[f64; 3]> = rotated.points.iter().map(|p| [p.x, p.y, p.z]).collect();
            let trial = optimize_mmff_conformer(&field, coords, max_iters);
            if trial.status != 0
                || !trial.energy.is_finite()
                || trial.energy >= best.energy - 1e-6
                || !trial.coords.iter().flatten().all(|value| value.is_finite())
            {
                continue;
            }
            let checked = Coords3D {
                points: trial
                    .coords
                    .iter()
                    .map(|&[x, y, z]| Point3::new(x, y, z))
                    .collect(),
            };
            let stereo = verify_stereo(mol, &checked);
            if stereo.n_violations() == 0 && stereo.n_unevaluable() == 0 {
                best = trial;
            }
        }
    }
    if best.energy
        >= candidates
            .iter()
            .filter(|candidate| candidate.status == 0 && candidate.energy.is_finite())
            .map(|candidate| candidate.energy)
            .min_by(f64::total_cmp)
            .unwrap_or(f64::INFINITY)
            - 1e-6
    {
        return candidates;
    }
    let limit = candidates.len();
    select_quality_candidates(candidates, vec![best], limit)
}

fn heavy_neighbor(mol: &Molecule, atom: AtomIdx, excluded: AtomIdx) -> Option<AtomIdx> {
    mol.neighbors(atom)
        .map(|(neighbor, _)| neighbor)
        .filter(|&neighbor| neighbor != excluded && mol.atom(neighbor).element.atomic_number() > 1)
        .min_by_key(|neighbor| neighbor.0)
}
