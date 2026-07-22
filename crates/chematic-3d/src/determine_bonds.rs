//! Bond order inference from 3D coordinates.
//!
//! Provides [`determine_bonds`], a safe alternative to RDKit's `DetermineBonds`
//! (issue #7922, which can freeze the CPU on molecules with 50+ atoms via an
//! exhaustive backtracking search).
//!
//! # Algorithm
//! Three bounded phases — no recursion, no backtracking, O(n²) worst case:
//!
//! 1. **Connectivity**: single bond if `dist(a, b) ≤ cov_radius[a] + cov_radius[b] + tol`.
//! 2. **Bond order**: greedily upgrade Single → Double → Triple to satisfy atom valences.
//! 3. **Aromaticity**: delegate to [`chematic_perception::apply_aromaticity`] (SSSR + Hückel).
//!
//! # Precision
//! Bond order assignment is accurate only when **explicit hydrogen atoms** are included
//! in `atoms`. Without H, the algorithm cannot distinguish C=O from C-O (both look like
//! a single-bond O with remaining valence 1). With H, oxygen in ethanol has degree 2
//! (C-O + O-H) → remaining 0 → correctly stays Single.

use crate::coords::Point3;
use chematic_core::{Atom, AtomIdx, BondIdx, BondOrder, Element, Molecule, MoleculeBuilder};
use chematic_perception::apply_aromaticity;

/// Maximum number of atoms accepted by [`determine_bonds`].
///
/// Above this threshold the function returns [`DetermineError::TooManyAtoms`]
/// rather than running. This protects against the O(n²) connectivity phase
/// becoming slow (300 atoms → 44,850 pair checks, well within 1 ms).
pub const MAX_ATOMS: usize = 300;

/// Errors returned by [`determine_bonds`].
#[derive(Debug, Clone, PartialEq)]
pub enum DetermineError {
    /// Input contained more than [`MAX_ATOMS`] atoms.
    TooManyAtoms(usize),
    /// Input contained zero atoms.
    EmptyMolecule,
}

impl std::fmt::Display for DetermineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DetermineError::TooManyAtoms(n) => {
                write!(f, "molecule has {} atoms; maximum is {}", n, MAX_ATOMS)
            }
            DetermineError::EmptyMolecule => write!(f, "empty molecule"),
        }
    }
}

/// Infer bond connectivity and bond orders from 3D coordinates.
///
/// # Arguments
/// - `atoms`: slice of (element, position) pairs, **including explicit hydrogen**.
/// - `tolerance_angstrom`: extra tolerance added to covalent radius sum (default 0.40 Å).
///
/// # Returns
/// A kekulized `Molecule` with Single, Double, Triple, and Aromatic bonds.
/// Includes all input atoms (including H). Call [`chematic_chem::remove_hydrogens`]
/// before generating SMILES if implicit H is desired.
///
/// # Safety
/// Never freezes. All loops are O(n²). No recursion. No backtracking.
pub fn determine_bonds(
    atoms: &[(Element, Point3)],
    tolerance_angstrom: f64,
) -> Result<Molecule, DetermineError> {
    let n = atoms.len();
    if n == 0 {
        return Err(DetermineError::EmptyMolecule);
    }
    if n > MAX_ATOMS {
        return Err(DetermineError::TooManyAtoms(n));
    }

    // ── Phase 1: Build connectivity (all bonds initially Single) ──────────────
    let mut builder = MoleculeBuilder::new();
    for (element, _) in atoms {
        builder.add_atom(Atom::new(*element));
    }
    for i in 0..n {
        let (elem_i, pos_i) = &atoms[i];
        let r_i = elem_i.covalent_radius() as f64;
        for j in (i + 1)..n {
            let (elem_j, pos_j) = &atoms[j];
            let r_j = elem_j.covalent_radius() as f64;
            let threshold = r_i + r_j + tolerance_angstrom;
            if pos_i.distance(pos_j) <= threshold {
                let _ = builder.add_bond(AtomIdx(i as u32), AtomIdx(j as u32), BondOrder::Single);
            }
        }
    }
    let mol = builder.build();

    // ── Phase 2: Upgrade bond orders to satisfy atom valences ─────────────────
    //
    // remaining_valence[i] = lowest_normal_valence_≥_current_degree − current_degree
    // H always has remaining=0 (degree=1, valence=1); it never participates in upgrades.
    // Transition metals return empty normal_valences() → remaining=0 → skipped.
    let bond_count = mol.bond_count();

    let mut remaining: Vec<i32> = (0..n)
        .map(|i| {
            let idx = AtomIdx(i as u32);
            let atom = mol.atom(idx);
            let degree = mol.neighbors(idx).count() as i32;
            let valences = atom.element.normal_valences();
            if valences.is_empty() {
                return 0;
            }
            let best = valences
                .iter()
                .find(|&&v| v as i32 >= degree)
                .copied()
                .unwrap_or(*valences.last().unwrap());
            (best as i32) - degree
        })
        .collect();

    // Track bond orders in a mutable Vec; rebuild molecule once at the end.
    let mut bond_orders: Vec<BondOrder> = (0..bond_count)
        .map(|i| mol.bond(BondIdx(i as u32)).order)
        .collect();

    // Precompute bond endpoints: they are immutable throughout Phase 2,
    // so fetching once avoids repeated mol.bond() calls in the hot loop.
    let endpoints: Vec<(AtomIdx, AtomIdx)> = (0..bond_count)
        .map(|i| {
            let b = mol.bond(BondIdx(i as u32));
            (b.atom1, b.atom2)
        })
        .collect();

    let mut changed = true;
    while changed {
        changed = false;
        for (bi, &(a_idx, b_idx)) in endpoints.iter().enumerate() {
            let a = a_idx.0 as usize;
            let b = b_idx.0 as usize;
            match bond_orders[bi] {
                BondOrder::Single if remaining[a] > 0 && remaining[b] > 0 => {
                    bond_orders[bi] = BondOrder::Double;
                    remaining[a] -= 1;
                    remaining[b] -= 1;
                    changed = true;
                }
                BondOrder::Double if remaining[a] > 0 && remaining[b] > 0 => {
                    bond_orders[bi] = BondOrder::Triple;
                    remaining[a] -= 1;
                    remaining[b] -= 1;
                    changed = true;
                }
                _ => {}
            }
        }
    }

    // Rebuild with upgraded bond orders (endpoints reused — no extra mol.bond() calls)
    let mut final_builder = MoleculeBuilder::new();
    for (_, atom) in mol.atoms() {
        final_builder.add_atom(atom.clone());
    }
    for (bi, &(a_idx, b_idx)) in endpoints.iter().enumerate() {
        let _ = final_builder.add_bond(a_idx, b_idx, bond_orders[bi]);
    }
    let mol_kekulized = final_builder.build();

    // ── Phase 3: Aromaticity perception ───────────────────────────────────────
    // Delegates to chematic_perception::apply_aromaticity (SSSR + Hückel 4n+2).
    // No-op for non-aromatic molecules; handles benzene, pyridine, pyrrole, furan, etc.
    let mol_final = apply_aromaticity(&mol_kekulized);

    Ok(mol_final)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::BondOrder;

    #[test]
    fn test_methane_all_single() {
        // Tetrahedral CH4: C at origin + 4 H at ~1.09 Å
        let atoms = vec![
            (Element::C, Point3::new(0.000, 0.000, 0.000)),
            (Element::H, Point3::new(0.629, 0.629, 0.629)),
            (Element::H, Point3::new(-0.629, -0.629, 0.629)),
            (Element::H, Point3::new(-0.629, 0.629, -0.629)),
            (Element::H, Point3::new(0.629, -0.629, -0.629)),
        ];
        let mol = determine_bonds(&atoms, 0.40).unwrap();
        assert_eq!(mol.atom_count(), 5);
        assert_eq!(mol.bond_count(), 4, "CH4 has 4 bonds");
        for (_, bond) in mol.bonds() {
            assert_eq!(bond.order, BondOrder::Single, "all CH4 bonds are single");
        }
    }

    #[test]
    fn test_ethanol_all_single() {
        // CCO with explicit H: no spurious double bonds from C or O.
        // Approximate geometry: C-C=1.54Å, C-O=1.43Å, O-H=0.96Å
        let atoms = vec![
            (Element::C, Point3::new(0.000, 0.000, 0.000)), // CH3
            (Element::C, Point3::new(1.540, 0.000, 0.000)), // CH2
            (Element::O, Point3::new(2.060, 1.200, 0.000)), // OH
            (Element::H, Point3::new(-0.380, 1.030, 0.000)),
            (Element::H, Point3::new(-0.380, -0.515, 0.890)),
            (Element::H, Point3::new(-0.380, -0.515, -0.890)),
            (Element::H, Point3::new(1.920, -0.515, 0.890)),
            (Element::H, Point3::new(1.920, -0.515, -0.890)),
            (Element::H, Point3::new(2.940, 1.170, 0.000)), // H on O
        ];
        let mol = determine_bonds(&atoms, 0.40).unwrap();
        assert_eq!(mol.atom_count(), 9);
        for (_, bond) in mol.bonds() {
            assert_eq!(
                bond.order,
                BondOrder::Single,
                "ethanol has only single bonds"
            );
        }
    }

    #[test]
    fn test_benzene_aromatic_bonds() {
        // Regular hexagon: C-C ≈ 1.40 Å, C-H ≈ 1.09 Å
        use std::f64::consts::PI;
        let mut atoms: Vec<(Element, Point3)> = Vec::new();
        for i in 0..6 {
            let angle = i as f64 * PI / 3.0;
            atoms.push((
                Element::C,
                Point3::new(1.40 * angle.cos(), 1.40 * angle.sin(), 0.0),
            ));
        }
        for i in 0..6 {
            let angle = i as f64 * PI / 3.0;
            atoms.push((
                Element::H,
                Point3::new(2.49 * angle.cos(), 2.49 * angle.sin(), 0.0),
            ));
        }
        let mol = determine_bonds(&atoms, 0.40).unwrap();
        assert_eq!(mol.atom_count(), 12);
        assert_eq!(mol.bond_count(), 12, "6 C-C + 6 C-H");
        let aromatic_count = mol
            .bonds()
            .filter(|(_, b)| b.order == BondOrder::Aromatic)
            .count();
        assert_eq!(aromatic_count, 6, "benzene: 6 aromatic C-C bonds");
        let single_count = mol
            .bonds()
            .filter(|(_, b)| b.order == BondOrder::Single)
            .count();
        assert_eq!(single_count, 6, "benzene: 6 C-H single bonds");
    }

    #[test]
    fn test_glycine_carbonyl_double_bond() {
        // Glycine NCC(=O)O — carboxyl group with explicit H on OH only.
        // The =O (no H) should be assigned Double; the O-H gets Single.
        // N-C-C(=O)-O with approximate geometry:
        let atoms = vec![
            (Element::N, Point3::new(-1.600, 0.000, 0.000)), // NH2
            (Element::C, Point3::new(-0.500, 0.600, 0.000)), // alpha-C
            (Element::C, Point3::new(0.800, 0.000, 0.000)),  // carboxyl C
            (Element::O, Point3::new(0.900, -1.200, 0.000)), // =O (no H)
            (Element::O, Point3::new(1.800, 0.800, 0.000)),  // -OH
            (Element::H, Point3::new(2.700, 0.400, 0.000)),  // H on -OH
            (Element::H, Point3::new(-2.000, 0.600, 0.000)), // H on N
            (Element::H, Point3::new(-2.000, -0.600, 0.000)), // H on N
            (Element::H, Point3::new(-0.600, 1.200, 0.890)), // H on alpha-C
            (Element::H, Point3::new(-0.600, 1.200, -0.890)), // H on alpha-C
        ];
        let mol = determine_bonds(&atoms, 0.40).unwrap();

        // Find carboxyl C (index 2) and its bonds to O
        let carboxyl_c = AtomIdx(2);
        let o_bond_orders: Vec<BondOrder> = mol
            .neighbors(carboxyl_c)
            .filter(|(nb, _)| mol.atom(*nb).element == Element::O)
            .map(|(_, bond_idx)| mol.bond(bond_idx).order)
            .collect();

        assert!(
            o_bond_orders.contains(&BondOrder::Double),
            "carboxyl C=O must be Double; got {:?}",
            o_bond_orders
        );
        assert!(
            o_bond_orders.contains(&BondOrder::Single),
            "carboxyl C-OH must be Single; got {:?}",
            o_bond_orders
        );
    }

    #[test]
    fn test_too_many_atoms_error() {
        let atoms: Vec<(Element, Point3)> = (0..(MAX_ATOMS + 1))
            .map(|i| (Element::C, Point3::new(i as f64 * 2.0, 0.0, 0.0)))
            .collect();
        match determine_bonds(&atoms, 0.40) {
            Err(DetermineError::TooManyAtoms(n)) => assert_eq!(n, MAX_ATOMS + 1),
            other => panic!("expected TooManyAtoms, got {:?}", other.err()),
        }
    }

    #[test]
    fn test_empty_molecule_error() {
        match determine_bonds(&[], 0.40) {
            Err(DetermineError::EmptyMolecule) => {}
            other => panic!("expected EmptyMolecule, got {:?}", other.err()),
        }
    }

    #[test]
    fn test_exact_max_atoms_ok() {
        // MAX_ATOMS atoms placed far apart (no bonds) should succeed.
        let atoms: Vec<(Element, Point3)> = (0..MAX_ATOMS)
            .map(|i| (Element::C, Point3::new(i as f64 * 10.0, 0.0, 0.0)))
            .collect();
        let mol = determine_bonds(&atoms, 0.40).unwrap();
        assert_eq!(mol.atom_count(), MAX_ATOMS);
        assert_eq!(mol.bond_count(), 0, "far-apart atoms have no bonds");
    }

    #[test]
    fn test_formaldehyde_double_bond() {
        // H2C=O: C should have Double bond to O.
        // Without H: C degree=1 (only O), remaining=3; O degree=1, remaining=1 → upgrades to Double
        //   (O's remaining hits 0, blocking further upgrade — accidentally correct for C=O but
        //    misleading for C≡N where N remaining=2 would allow Triple).
        // With H: C degree=3 (2×H + 1×O), remaining=1; O degree=1, remaining=1 → upgrades to Double (correct)
        let atoms = vec![
            (Element::C, Point3::new(0.000, 0.000, 0.000)),
            (Element::O, Point3::new(0.000, 1.200, 0.000)),
            (Element::H, Point3::new(0.940, -0.540, 0.000)),
            (Element::H, Point3::new(-0.940, -0.540, 0.000)),
        ];
        let mol = determine_bonds(&atoms, 0.40).unwrap();
        let c_idx = AtomIdx(0);
        let o_bond_order = mol
            .neighbors(c_idx)
            .find(|(nb, _)| mol.atom(*nb).element == Element::O)
            .map(|(_, bi)| mol.bond(bi).order)
            .unwrap();
        assert_eq!(o_bond_order, BondOrder::Double, "H2C=O: C-O must be Double");
    }

    #[test]
    fn test_te_c_pair_now_upgrades_to_double_dependent_path_change() {
        // Dependent-path finding for the Te normal_valences() fix (chematic-core
        // element.rs, atomic number 52 [2, 4, 6]): this crate's Phase 2 upgrade
        // loop (lines ~102-118 above) treats an element with an *empty*
        // normal_valences() as always having remaining_valence=0, so it never
        // participates in a Single->Double upgrade regardless of geometry (the
        // algorithm is valence-driven, not distance-driven, for order upgrades —
        // see test_formaldehyde_double_bond's comment on the same behavior).
        //
        // BEFORE this fix (empirically confirmed via `git stash` on element.rs
        // alone, then re-running this exact fixture): a bare Te-C pair connects
        // (Phase 1, distance within covalent-radius-sum+tolerance) but the bond
        // stays Single forever, because remaining[Te] was hard-coded to 0.
        //
        // AFTER this fix: Te gets valences=[2,4,6], degree=1, remaining=2-1=1>0,
        // matching C's remaining=4-1=3>0, so Phase 2's greedy loop upgrades the
        // bond to Double (then stops, since remaining[Te] hits 0).
        //
        // This is a real, reportable behavior change for 3D Te fixtures with no
        // other bonds to Te — flagged here rather than silently accepted, per
        // the task's explicit instruction. It only fires for isolated/low-degree
        // Te atoms; Te atoms with degree >= 2 that already satisfy valence 2
        // (e.g. a bent Te with two single bonds, the common case) are unaffected.
        let atoms = vec![
            (Element::TE, Point3::new(0.000, 0.000, 0.000)),
            (Element::C, Point3::new(2.000, 0.000, 0.000)),
        ];
        let mol = determine_bonds(&atoms, 0.40).unwrap();
        assert_eq!(mol.bond_count(), 1, "Te-C pair must connect");
        let order = mol.bond(BondIdx(0)).order;
        assert_eq!(
            order,
            BondOrder::Double,
            "AFTER the Te valence-table fix, a bare degree-1 Te-C pair upgrades \
             to Double (BEFORE the fix it stayed Single — Te's remaining valence \
             was always 0)."
        );
    }
}
