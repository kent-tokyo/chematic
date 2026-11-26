use super::*;

fn graph(atoms: &[(u32, i32)], bonds: &[(usize, usize, BondType)]) -> Mol {
    let mut mol = Mol::default();
    for &(anum, charge) in atoms {
        let mut atom = Atom::new(anum);
        atom.charge = charge;
        atom.no_implicit = true;
        mol.add_atom(atom);
    }
    for &(a, b, bt) in bonds {
        mol.add_bond(Bond::new(a, b, bt));
    }
    mol
}

#[test]
fn protonated_oxygen_transfer_keeps_charge_and_hydrogen() {
    // The corresponding protonated nitromethane InChI is read by RDKit
    // 2025.09.3 as C[N+](=O)O, with one retained proton and charge +1.
    let mut molecule = graph(
        &[(7, 0), (8, 0), (8, 1), (6, 0)],
        &[(0, 1, D), (0, 2, D), (0, 3, S)],
    );
    molecule.atoms[2].num_explicit_hs = 1;
    assert!(v5n_5(&mut molecule, 0, 8));
    assert_eq!(
        (
            molecule.atoms[0].charge,
            molecule.atoms[1].charge,
            molecule.atoms[2].charge
        ),
        (1, 0, 0)
    );
    assert_eq!(
        (
            molecule.atoms[1].num_explicit_hs,
            molecule.atoms[2].num_explicit_hs
        ),
        (1, 0)
    );
    assert_eq!((molecule.bonds[0].bt, molecule.bonds[1].bt), (S, D));
}
