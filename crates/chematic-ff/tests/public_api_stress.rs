//! Broad force-field public-API coverage over the checked-in stress corpus.

use chematic_ff::rdkit_mmff::{MmffVariant, RdkitMmffField};
use chematic_ff::rdkit_uff::RdkitUffField;
use chematic_ff::*;

const STRESS: &str = include_str!("../../../validation/cosmolkit_comparison/corpora/stress-v1.smi");

fn linear_coords(atom_count: usize) -> Vec<[f64; 3]> {
    (0..atom_count)
        .map(|i| {
            let i = i as f64;
            [1.35 * i, 0.2 * (i % 3.0), 0.15 * (i % 5.0)]
        })
        .collect()
}

#[test]
fn force_field_surface_handles_special_chemistry() {
    let mut exercised = 0usize;
    let mut energetic = 0usize;
    for line in STRESS.lines() {
        let smiles = line.split_whitespace().next().unwrap_or("");
        let Ok(mol) = chematic_smiles::parse(smiles) else {
            continue;
        };
        if mol.atom_count() == 0 || mol.atom_count() > 120 {
            continue;
        }
        exercised += 1;

        let dreiding = assign_dreiding_types(&mol);
        for &ty in &dreiding {
            let _ = dreiding_vdw(ty);
            let _ = dreiding_angle(ty);
            let _ = dreiding_torsion_barrier(ty, ty);
        }
        for bond in mol.bonds().map(|(_, bond)| bond) {
            let first = dreiding[bond.atom1.0 as usize];
            let second = dreiding[bond.atom2.0 as usize];
            let _ = dreiding_bond_len(first, second, bond.order);
        }

        let numeric = assign_mmff94_numeric_types(&mol);
        let rings = chematic_perception::find_symmetrized_sssr(&mol);
        let _ = compute_mmff94_aromatic_view(&mol, rings.rings());
        let _ = assign_mmff94_numeric_types_with_view(&mol);
        let _ = mmff94_charges_numeric(&mol);
        let _ = mmff94_charges_bci(&mol);
        if let Ok(types) = assign_mmff94_types(&mol) {
            let matrix = ElectrostaticMatrix::new(&mol, &types);
            let _ = matrix.get_scaling(0, mol.atom_count().saturating_sub(1));
            let _ = matrix.get_scaling(usize::MAX, 0);
            let _ = matrix.is_1_4(0, mol.atom_count().saturating_sub(1));
            let _ = matrix.is_1_4(usize::MAX, 0);
            for &ty in &types {
                let _ = mmff94_formal_charge(ty);
                let _ = mmff94_charge_params(ty);
                let _ = mmff94_vdw_params(ty);
                let _ = bci(ty, ty);
            }
            for bond in mol.bonds().map(|(_, bond)| bond) {
                let first = types[bond.atom1.0 as usize];
                let second = types[bond.atom2.0 as usize];
                let _ = mmff94_bond_params(first, second, bond.order);
                let _ = mmff94_bond_dipole(first, second, 1.5);
                let _ = mmff94_electrostatic_scaling_1_4(first, second);
            }
        }

        if let Ok(types) = &numeric {
            for &ty in types {
                let _ = mmff94_numeric_type_info(ty);
                let _ = pbci_for(ty);
                let _ = mmff94_vdw_energy(ty);
                let _ = mmff94_vdw_combined(ty, ty);
                let _ = mmff94_bond_energy(0, ty, ty);
            }
        }

        let uff = assign_uff_types(&mol);
        let labels = chematic_ff::rdkit_uff::rdkit_uff_atom_labels(&mol);
        let _ = chematic_ff::rdkit_uff::rdkit_uff_has_all_params(&mol);
        for pair in labels.windows(2) {
            let _ = chematic_ff::rdkit_uff::rdkit_uff_bond_rest_length(1.0, &pair[0], &pair[1]);
        }

        if energetic < 256 && mol.atom_count() <= 64 {
            energetic += 1;
            let coords = linear_coords(mol.atom_count());
            let _ = uff_total_energy(&mol, &uff, &coords);
            let _ = minimize_uff(&mol, &uff, coords.clone(), 1);
            let _ = minimize_uff_with_constraint(&mol, &uff, coords.clone(), 1, |_| true);

            if let Ok(model) = Mmff94EnergyModel::new(&mol) {
                let _ = model.energy(&coords);
                let _ = model.energy_breakdown(&coords);
                let _ = model.nonbonded_gradient(&coords);
                let _ = model.cutoff_nonbonded_energy(&coords);
                let _ = model.cutoff_nonbonded_gradient(&coords);
                let _ = model.bond_angle_gradient(&coords);
                let _ = model.stretch_bend_gradient(&coords);
                let _ = model.torsion_gradient(&coords);
                let _ = model.oop_gradient(&coords);
                let _ = model.bounded_analytic_gradient(&coords);
                let _ = model.bounded_analytic_gradient_cutoff(&coords);
            }

            let flat: Vec<f64> = coords.iter().flatten().copied().collect();
            if let Ok(field) = RdkitMmffField::with_terms(
                &mol,
                &coords,
                MmffVariant::Mmff94s,
                100.0,
                true,
                [true; 7],
            ) {
                let _ = field.num_atoms();
                let _ = field.energy(&flat);
                let _ = field.energy_with_distances(&flat, &flat);
                let mut gradient = vec![0.0; flat.len()];
                field.gradient(&flat, &mut gradient);
            }
            let field = RdkitUffField::new(&mol, &coords, 100.0, true);
            let _ = field.energy(&flat);
            let _ = field.energy_with_distances(&flat, &flat);
            let mut gradient = vec![0.0; flat.len()];
            field.gradient(&flat, &mut gradient);
        }
    }

    assert!(
        exercised >= 4_000,
        "stress corpus unexpectedly shrank: {exercised}"
    );
    assert_eq!(energetic, 256);
}

#[test]
fn batch_and_registry_boundaries_are_total() {
    let mut batch = MMFF94BatchProperties::new();
    assert!(batch.is_empty());
    assert!(batch.average_charges().is_empty());
    assert!(batch.charge_variance().is_empty());

    let mol = chematic_smiles::parse("CCO").unwrap();
    let types = assign_mmff94_types(&mol).unwrap();
    let matrix = ElectrostaticMatrix::new(&mol, &types);
    batch.push(vec![0.1, -0.2, 0.1], vec![(0, 1, 0.2)], matrix.clone());
    batch.push(vec![0.2, -0.4, 0.2], vec![(1, 2, 0.3)], matrix);
    assert_eq!(batch.len(), 2);
    assert_eq!(batch.average_charges().len(), 3);
    assert_eq!(batch.charge_variance().len(), 3);

    for id in 0..=u8::MAX {
        let _ = mmff94_numeric_type_info(id);
        let _ = mmff94_vdw_energy(id);
        let _ = mmff94_bond_energy(0, id, id);
        let _ = mmff94_angle_energy(0, id, id, id);
        let _ = mmff94_torsion_energy(0, id, id, id, id);
        let _ = mmff94_oop(id, id, id, id);
        let _ = mmff94_stbn_type_only(0, id, id, id);
    }
}
