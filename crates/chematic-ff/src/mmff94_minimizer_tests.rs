use super::*;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mmff94_energy::mmff94_stbn_type_only;
    use crate::mmff94_energy::mmff94_torsion_energy;
    use crate::mmff94_energy::{mmff94_angle_energy, mmff94_bond_energy};
    use crate::mmff94_numeric::{
        assign_mmff94_numeric_types, assign_mmff94_numeric_types_with_view,
    };
    use chematic_core::molecule::MoleculeBuilder;
    use chematic_core::{Atom, BondOrder, Element};

    fn methane_mol() -> (Molecule, Vec<[f64; 3]>) {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let h1 = b.add_atom(Atom::new(Element::H));
        let h2 = b.add_atom(Atom::new(Element::H));
        let h3 = b.add_atom(Atom::new(Element::H));
        let h4 = b.add_atom(Atom::new(Element::H));
        b.add_bond(c, h1, BondOrder::Single).unwrap();
        b.add_bond(c, h2, BondOrder::Single).unwrap();
        b.add_bond(c, h3, BondOrder::Single).unwrap();
        b.add_bond(c, h4, BondOrder::Single).unwrap();
        let mol = b.build();
        let coords = vec![
            [0.0, 0.0, 0.0],
            [0.630, 0.630, 0.630],
            [-0.630, -0.630, 0.630],
            [-0.630, 0.630, -0.630],
            [0.630, -0.630, -0.630],
        ];
        (mol, coords)
    }

    #[test]
    fn prepared_energy_model_matches_one_shot_api() {
        let (mol, coords) = methane_mol();
        let model = Mmff94EnergyModel::new(&mol).expect("prepare MMFF94 model");
        let one_shot = mmff94_energy_breakdown(&mol, &coords).expect("one-shot energy");
        let prepared = model.energy_breakdown(&coords);
        assert!((prepared.total - one_shot.total).abs() < 1e-12);
        assert!((model.energy(&coords) - one_shot.total).abs() < 1e-12);
    }

    fn butane_backbone() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let c0 = b.add_atom(Atom::new(Element::C));
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c0, c1, BondOrder::Single).unwrap();
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.add_bond(c2, c3, BondOrder::Single).unwrap();
        b.build()
    }

    #[test]
    fn energy_is_finite_for_methane() {
        let (mol, coords) = methane_mol();
        let e = mmff94_total_energy(&mol, &coords).expect("energy");
        assert!(e.is_finite(), "energy={}", e);
    }

    #[test]
    fn torsion_differs_by_conformation() {
        let mol = butane_backbone();
        let types = assign_mmff94_numeric_types(&mol).expect("types");
        // Gauche: ~60° central dihedral
        let coords_gauche = vec![
            [0.0, 0.0, 0.0_f64],
            [1.508, 0.0, 0.0],
            [2.016, 1.192, 0.688],
            [3.524, 1.192, 0.688],
        ];
        // Anti: ~180° central dihedral
        let coords_anti = vec![
            [0.0, 0.0, 0.0_f64],
            [1.508, 0.0, 0.0],
            [3.016, 0.0, 0.0],
            [4.524, 0.0, 0.0],
        ];
        let e_gauche = torsion_energy(&mol, &coords_gauche, &types);
        let e_anti = torsion_energy(&mol, &coords_anti, &types);
        assert!(e_gauche.is_finite());
        assert!(e_anti.is_finite());
        assert!(
            (e_gauche - e_anti).abs() > 1e-6,
            "torsion must differ: gauche={}, anti={}",
            e_gauche,
            e_anti
        );
    }

    #[test]
    fn vdw_more_repulsive_at_short_range() {
        let mol = butane_backbone();
        let types = assign_mmff94_numeric_types(&mol).expect("types");
        // Atoms 0 and 3 are 1-4 (not excluded from vdW)
        let coords_close = vec![
            [0.0, 0.0, 0.0_f64],
            [1.5, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.5, 0.0, 0.0], // atom 3 very close to atom 0
        ];
        let coords_far = vec![
            [0.0, 0.0, 0.0_f64],
            [1.5, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [8.0, 0.0, 0.0],
        ];
        let e_close = vdw_energy(&mol, &coords_close, &types);
        let e_far = vdw_energy(&mol, &coords_far, &types);
        assert!(e_close.is_finite());
        assert!(e_far.is_finite());
        assert!(e_close > e_far, "close={} should > far={}", e_close, e_far);
    }

    #[test]
    fn buffered_14_7_matches_rdkit_reference_pair() {
        // RDKit MMFF94 type 37/37 at this shared-coordinate contact reports
        // 4.2464853537 kcal/mol. The former t^7 * (t^7 - 2) shortcut returned
        // 11.9247 and inflated aromatic non-bonded energies.
        let energy =
            mmff94_vdw_energy_value(2.7745809246687023, 4.193078986609192, 0.06779699304291385);
        assert!((energy - 4.246485353735497).abs() < 1e-10);
    }

    /// First SSSR ring of exactly `size` atoms, in ring (bonded-consecutive)
    /// order, or `None` if the molecule has no ring of that size.
    fn first_ring_of_size(mol: &Molecule, size: usize) -> Option<Vec<usize>> {
        let ring_set = find_sssr(mol);
        ring_set
            .rings()
            .iter()
            .find(|r| r.len() == size)
            .map(|r| r.iter().map(|a| a.0 as usize).collect())
    }

    // ── FF-1: bond_type_for (#173) ─────────────────────────────────────────

    #[test]
    fn bond_type_for_ethene_double_bond_is_type_0_not_1() {
        // Issue #173's primary repro: a real C=C double bond must never be
        // routed to the sbmb bond_type=1 row (r0=1.430, wrong) instead of the
        // correct bond_type=0 row (r0=1.333).
        let mol = chematic_smiles::parse("C=C").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let (_, bond) = mol.bonds().next().unwrap();
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        assert_eq!(types[i], 2, "vinylic carbon should be type 2");
        assert_eq!(types[j], 2);
        let bt = bond_type_for(types[i], types[j], bond.order);
        assert_eq!(bt, 0, "C=C double bond must classify as bond_type 0");
        let p = mmff94_bond_energy(bt, types[i], types[j]).expect("(0,2,2) row must exist");
        assert!((p.r0 - 1.333).abs() < 1e-6, "r0={} should be 1.333", p.r0);
    }

    #[test]
    fn bond_type_for_acetaldehyde_sp3_carbonyl_single_bond_is_type_0() {
        // Issue #173's secondary repro: CC=O's Cα(sp3)-C(carbonyl) single
        // bond must resolve to the real (0,1,3) row, not silently miss via a
        // wrongly-routed bond_type=1 lookup (which has no (1,1,3) row at all).
        let mol = chematic_smiles::parse("CC=O").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        // atom 0 = CH3 (type 1), atom 1 = C=O carbon (type 3)
        assert_eq!(types[0], 1);
        assert_eq!(types[1], 3);
        let order = bond_order_between(&mol, 0, 1);
        let bt = bond_type_for(types[0], types[1], order);
        assert_eq!(
            bt, 0,
            "sp3-carbonyl single bond must classify as bond_type 0"
        );
        let p = mmff94_bond_energy(bt, types[0], types[1]).expect("(0,1,3) row must exist");
        assert!((p.r0 - 1.492).abs() < 1e-6, "r0={} should be 1.492", p.r0);
    }

    #[test]
    fn ethene_energy_minimum_is_near_1_333_angstrom() {
        // Energy-perturbation reproduction of issue #173's own measurement:
        // mmff94_total_energy("C=C")'s minimum used to sit at r≈1.430 Å; it
        // must now sit near the chemically correct 1.333 Å.
        let mol = chematic_smiles::parse("C=C").unwrap();
        let scan = |r: f64| {
            let coords = vec![[0.0, 0.0, 0.0], [r, 0.0, 0.0]];
            mmff94_total_energy(&mol, &coords).unwrap()
        };
        let e_correct = scan(1.333);
        let e_old_wrong_minimum = scan(1.430);
        let e_far = scan(1.6);
        assert!(
            e_correct < e_old_wrong_minimum,
            "energy at true r0=1.333 ({e_correct}) should be lower than at the old wrong \
             minimum 1.430 ({e_old_wrong_minimum})"
        );
        assert!(
            e_correct < e_far,
            "energy at 1.333 ({e_correct}) should be lower than further out at 1.6 ({e_far})"
        );
    }

    #[test]
    fn acetaldehyde_bond_stretch_now_changes_energy() {
        // Issue #173's exact repro: stretching the Cα-C(carbonyl) bond by a
        // full 1 Å used to produce dE = 0.0 (silent missing parameter).
        let mol = chematic_smiles::parse("CC=O").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        // Build coords: C0 at origin, C1 (carbonyl C) along +x, O at a fixed
        // offset from C1 so the C=O bond itself is untouched by the scan.
        let coords_at = |r: f64| vec![[0.0, 0.0, 0.0], [r, 0.0, 0.0], [r + 1.2, 0.9, 0.0]];
        let e_short = bond_energy(&mol, &coords_at(1.5), &types);
        let e_long = bond_energy(&mol, &coords_at(2.5), &types);
        assert!(
            (e_long - e_short).abs() > 1e-6,
            "bond stretch must now change bond energy: short={e_short} long={e_long} (was \
             identical pre-fix per issue #173)"
        );
    }

    #[test]
    fn bond_gradient_restores_ethene_toward_1_333() {
        // Gradient test using this file's own existing finite-difference
        // pattern (`compute_gradient`, already used by both minimizers).
        let mol = chematic_smiles::parse("C=C").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let charges = mmff94_charges_numeric(&mol).unwrap_or_else(|_| vec![0.0; 2]);
        let rings = find_sssr(&mol);
        // Stretched well past the true r0=1.333 minimum.
        let coords = vec![[0.0, 0.0, 0.0_f64], [1.6, 0.0, 0.0]];
        let grad = compute_gradient(&mol, &coords, &types, &charges, rings.rings(), 1e-4);
        // Restoring force on atom 1 should point back toward atom 0 (-x),
        // i.e. dE/dx > 0 at this atom (energy increases as x increases further).
        assert!(
            grad[1][0] > 0.0,
            "gradient on stretched atom should point back toward equilibrium: grad_x={}",
            grad[1][0]
        );
    }

    #[test]
    fn benzene_ring_bonds_type_37_resolve_to_the_rdkit_verified_row() {
        // Issue #227 Phase 1B-0: atom typing fixed to the real Halgren/RDKit
        // numbering (type 37 = CB, not the old wrong 63 = C5A), AND
        // `bond_type_for` fixed to treat `BondOrder::Aromatic` like
        // double/triple (BT=0 unconditionally), not like `Single` run
        // through the MLTB_TYPES-AND check. Both fixes were required
        // together -- fixing only the atom typer would have made every
        // benzene ring bond resolve to the WRONG (1, 37, 37, kb=5.178,
        // r0=1.436) row (biphenyl's *non-aromatic* inter-ring bond value)
        // instead of the correct (0, 37, 37, kb=5.573, r0=1.374) row --
        // still "coverage success," still silently wrong, the same failure
        // shape as the `furan` finding this whole PR exists to close.
        // Cross-checked directly against a live RDKit oracle
        // (`props.GetMMFFBondStretchParams`), not assumed:
        // `AllChem.MMFFGetMoleculeProperties` on benzene returns exactly
        // `(bondType=0, kb=5.573, r0=1.374)` for every ring bond.
        let mol = chematic_smiles::parse("c1ccccc1").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        assert_eq!(mol.bonds().count(), 6);
        for (_, bond) in mol.bonds() {
            let ti = types[bond.atom1.0 as usize];
            let tj = types[bond.atom2.0 as usize];
            assert_eq!(
                ti, 37,
                "benzene carbon must type as 37 (CB), matching RDKit"
            );
            assert_eq!(tj, 37);
            assert_eq!(bond.order, BondOrder::Aromatic);
            let bt = bond_type_for(ti, tj, bond.order);
            assert_eq!(
                bt, 0,
                "aromatic ring bond between two type-37 atoms must resolve to bond_type=0, \
                 matching RDKit's GetMMFFBondStretchParams for benzene"
            );
            let params = mmff94_bond_energy(bt, ti, tj)
                .expect("every benzene ring bond must resolve to a real row, not silently miss");
            assert!(
                (params.r0 - 1.374).abs() < 1e-6,
                "benzene bond r0 must match RDKit's 1.374 A, got {}",
                params.r0
            );
            assert!(
                (params.kb - 5.573).abs() < 1e-6,
                "benzene bond kb must match RDKit's 5.573, got {}",
                params.kb
            );
        }
    }

    #[test]
    fn furan_c_c_bond_no_longer_collides_with_a_nitrogen_row() {
        // The exact regression this whole Phase 1B-0 PR exists to close
        // (issue #227's audit PR #235): before the fix, furan's ring
        // carbons were mistyped as 38/37 (real 38 = NPYD, pyridine-type
        // NITROGEN) and its C-C bond silently resolved to
        // `Some(BondEnergyParams { kb: 5.002, r0: 1.246 })` -- a real table
        // row, but one that belongs to a nitrogen-involved bond, not a
        // furan C-C bond. Pinning the old wrong value here so a future
        // reader can see exactly what this test prevents. After the fix,
        // furan's ring carbons must type as 63/64 (C5A/C5B, both carbon),
        // and their mutual bond must resolve to a genuinely carbon-carbon
        // row (or, honestly, `None` -- anything but a repeat of the old
        // nitrogen-row collision).
        let mol = chematic_smiles::parse("c1ccoc1").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        for (_, bond) in mol.bonds() {
            let a1 = mol.atom(bond.atom1);
            let a2 = mol.atom(bond.atom2);
            if a1.element != chematic_core::Element::C || a2.element != chematic_core::Element::C {
                continue;
            }
            let ti = types[bond.atom1.0 as usize];
            let tj = types[bond.atom2.0 as usize];
            assert_ne!(
                (ti, tj),
                (38, 37),
                "must not regress to the old wrong furan C-C typing"
            );
            if let Some(params) = mmff94_bond_energy(bond_type_for(ti, tj, bond.order), ti, tj) {
                assert!(
                    (params.r0 - 1.246).abs() > 1e-6,
                    "furan C-C bond must not resolve to the old nitrogen-row r0=1.246, \
                     types=({ti},{tj})"
                );
            }
        }
    }

    // ── FF-1: torsion_type_for ring types 4/5 (#175) ───────────────────────

    #[test]
    fn torsion_type_for_cyclopentane_ring_is_type_5() {
        let mol = chematic_smiles::parse("C1CCCC1").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let ring = first_ring_of_size(&mol, 5).expect("cyclopentane must have a 5-ring");
        let (i, j, k, l) = (ring[0], ring[1], ring[2], ring[3]);
        let tt = torsion_type_for(&mol, i, j, k, l, types[i], types[j], types[k], types[l]);
        assert_eq!(
            tt, 5,
            "in-ring torsion of a 5-membered ring must classify as type 5"
        );
    }

    #[test]
    fn torsion_type_for_cyclobutane_ring_is_type_4() {
        let mol = chematic_smiles::parse("C1CCC1").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let ring = first_ring_of_size(&mol, 4).expect("cyclobutane must have a 4-ring");
        let (i, j, k, l) = (ring[0], ring[1], ring[2], ring[3]);
        let tt = torsion_type_for(&mol, i, j, k, l, types[i], types[j], types[k], types[l]);
        assert_eq!(
            tt, 4,
            "in-ring torsion of a 4-membered ring must classify as type 4"
        );
    }

    #[test]
    fn torsion_type_for_ring5_override_requires_a_type_1_atom() {
        // RDKit's real 5-ring override additionally requires at least one of
        // the four torsion atoms to be MMFF numeric type 1 (`ti==1 || tj==1
        // || tk==1 || tl==1`) -- a condition the prior SSSR-based ring check
        // (which this fix replaced entirely) had no equivalent of at all.
        // Same 5-ring topology as `torsion_type_for_cyclopentane_ring_is_type_5`
        // above (where all four real chematic-assigned types happen to be 1,
        // so the condition is trivially satisfied there) but with the four
        // types forced away from 1 here -- purely local ring-adjacency
        // geometry (`ring_size_4_or_5` returning 5) must NOT be sufficient on
        // its own.
        let mol = chematic_smiles::parse("C1CCCC1").unwrap();
        let ring = first_ring_of_size(&mol, 5).expect("cyclopentane must have a 5-ring");
        let (i, j, k, l) = (ring[0], ring[1], ring[2], ring[3]);
        let tt = torsion_type_for(&mol, i, j, k, l, 20, 20, 20, 20);
        assert_ne!(
            tt, 5,
            "5-ring override must not fire when none of the four atoms is MMFF type 1"
        );
    }

    // ── issue #227: torsion_type_for's corrected bond-type-based base
    // classification and empirically-required type-2 override ─────────────

    #[test]
    fn torsion_type_for_aromatic_jk_bond_is_type_0_despite_mltb_membership() {
        // The core of issue #227's torsion bug: a prior version of this
        // function classified the non-ring base case purely from atom-type
        // MLTB_TYPES membership, `(MLTB(tj), MLTB(tk)) -> 0/1/2`, completely
        // ignoring the j-k bond's own real MMFF bond order/type. Benzene's
        // ring carbons are all type 37 (in MLTB_TYPES), so the prior formula
        // would have classified any benzene torsion as type 2 (both
        // "sp2-sp2") -- RDKit's real `getMMFFTorsionType` instead uses
        // `bond_type_for`'s own result as the base code, and an aromatic
        // bond always gets bond_type 0 (same as a real double/triple bond),
        // regardless of the endpoints' own MLTB membership.
        let mol = chematic_smiles::parse("c1ccccc1").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        assert_eq!(types[0], 37, "benzene carbon should be type 37");
        assert!(
            MLTB_TYPES.binary_search(&37).is_ok(),
            "type 37 is in MLTB_TYPES -- the prior atom-type-membership formula would wrongly see this as conjugated"
        );
        let tt = torsion_type_for(&mol, 0, 1, 2, 3, types[0], types[1], types[2], types[3]);
        assert_eq!(
            tt, 0,
            "an aromatic j-k bond must classify as torsion type 0 via bond_type_for, not type 2 via MLTB membership"
        );
    }

    #[test]
    fn torsion_type_for_type_2_override_fires_when_flanking_bond_is_sbmb() {
        // Synthetic 4-atom chain i-j-k-l, all single bonds (topology only
        // matters for bond order / ring adjacency here -- the four MMFF
        // types are passed explicitly to isolate the classification
        // formula). ti=tj=2 (both in MLTB_TYPES) makes the i-j bond an sbmb
        // bond (`bond_type_ij` = 1, a single bond directly linking two
        // independently conjugation-capable atoms). tk=tl=1 (not in
        // MLTB_TYPES) makes the j-k bond's own base code 0 (neither a
        // multiple bond nor doubly-conjugated). RDKit's empirically-required
        // override (`bond_type_jk==0 && order_jk==Single && (bond_type_ij==1
        // || bond_type_kl==1)`, needed to pass RDKit's own CYGUAN01
        // regression test, not derivable from Halgren's MMFF.IV page 609
        // formula alone) must bump this to type 2.
        let mut b = MoleculeBuilder::new();
        let i = b.add_atom(Atom::new(Element::C));
        let j = b.add_atom(Atom::new(Element::C));
        let k = b.add_atom(Atom::new(Element::C));
        let l = b.add_atom(Atom::new(Element::C));
        b.add_bond(i, j, BondOrder::Single).unwrap();
        b.add_bond(j, k, BondOrder::Single).unwrap();
        b.add_bond(k, l, BondOrder::Single).unwrap();
        let mol = b.build();
        let tt = torsion_type_for(
            &mol,
            i.0 as usize,
            j.0 as usize,
            k.0 as usize,
            l.0 as usize,
            2,
            2,
            1,
            1,
        );
        assert_eq!(
            tt, 2,
            "base type 0 at j-k with a flanking sbmb (bond_type=1) i-j bond must override to type 2"
        );
    }

    #[test]
    fn torsion_type_for_no_override_when_neither_flank_is_sbmb() {
        // Same topology as the override test above, but with NEITHER flank
        // an sbmb bond (ti=tl=1, not in MLTB_TYPES) -- the override's
        // `bond_type_ij==1 || bond_type_kl==1` condition must not be
        // satisfied trivially just because bond_type_jk==0.
        let mut b = MoleculeBuilder::new();
        let i = b.add_atom(Atom::new(Element::C));
        let j = b.add_atom(Atom::new(Element::C));
        let k = b.add_atom(Atom::new(Element::C));
        let l = b.add_atom(Atom::new(Element::C));
        b.add_bond(i, j, BondOrder::Single).unwrap();
        b.add_bond(j, k, BondOrder::Single).unwrap();
        b.add_bond(k, l, BondOrder::Single).unwrap();
        let mol = b.build();
        let tt = torsion_type_for(
            &mol,
            i.0 as usize,
            j.0 as usize,
            k.0 as usize,
            l.0 as usize,
            1,
            2,
            1,
            1,
        );
        assert_eq!(
            tt, 0,
            "base type 0 must stand when neither flanking bond is an sbmb (bond_type=1) bond"
        );
    }

    #[test]
    fn torsion_type_for_exocyclic_substituent_is_not_ring_typed() {
        // Regression guard for the advisor-flagged failure mode: a torsion
        // whose central j-k bond is a ring bond but whose *terminal* atom is
        // an exocyclic substituent (methylcyclopentane's exocyclic methyl)
        // must NOT be misclassified as a ring torsion.
        let mol = chematic_smiles::parse("CC1CCCC1").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        // atom 0 = exocyclic methyl C, atom 1 = ring C bonded to atom 0.
        let ring = first_ring_of_size(&mol, 5).expect("methylcyclopentane must have a 5-ring");
        assert!(
            ring.contains(&1),
            "atom 1 should be the substituted ring atom"
        );
        let j = 1usize;
        let k = *ring
            .iter()
            .find(|&&a| a != 1 && mol_bonded(&mol, j, a))
            .unwrap();
        let i = 0usize; // exocyclic methyl carbon, NOT a ring member
        let l = *ring
            .iter()
            .find(|&&a| a != j && a != k && mol_bonded(&mol, k, a))
            .unwrap();
        let tt = torsion_type_for(&mol, i, j, k, l, types[i], types[j], types[k], types[l]);
        assert_ne!(
            tt, 5,
            "torsion with an exocyclic terminal atom must not be classified as a ring torsion"
        );
    }

    fn mol_bonded(mol: &Molecule, a: usize, b: usize) -> bool {
        mol.bond_between(AtomIdx(a as u32), AtomIdx(b as u32))
            .is_some()
    }

    // ── torsion_no_term_by_design / NoTermByDesign (issue #227 Phase 1) ─────

    #[test]
    fn torsion_no_term_by_design_true_for_linear_central_atom_either_side() {
        // Type 4 (CSP, acetylenic/nitrile C) and 53 (=N=, cumulated
        // azide/diazo N) are the two `linear: true` types this corpus
        // exercises; type 61 is also linear per the registry but unreachable
        // in this corpus.
        assert!(torsion_no_term_by_design(4, 37), "linear on the j side");
        assert!(torsion_no_term_by_design(37, 53), "linear on the k side");
        assert!(torsion_no_term_by_design(4, 53), "linear on both sides");
    }

    #[test]
    fn torsion_no_term_by_design_false_when_neither_central_atom_is_linear() {
        assert!(!torsion_no_term_by_design(1, 1)); // sp3 C - sp3 C
        assert!(!torsion_no_term_by_design(37, 3)); // aromatic C - C=O (caffeine's own case)
    }

    #[test]
    fn torsion_no_term_by_design_unknown_type_fails_closed_to_false() {
        // A type absent from the registry must never be treated as linear
        // (would incorrectly suppress a real coverage-gap report).
        assert!(!torsion_no_term_by_design(250, 1));
    }

    #[test]
    fn nitrile_torsion_has_no_table_row_and_is_flagged_no_term_by_design() {
        // The exact chembl_tier_b_0001 shape (issue #227's 2 real
        // table_gap instances): Ar-Ar-C#N, central k = type 4 (CSP). Oracle-
        // confirmed RDKit's own GetMMFFTorsionParams also returns None here
        // (rdkit==2026.03.4) -- this is genuinely absent from the table on
        // both sides, not a routing bug.
        assert!(mmff94_torsion_energy(1, 37, 37, 4, 42).is_none());
        assert!(torsion_no_term_by_design(37, 4));
    }

    #[test]
    fn cumulated_azide_torsion_has_no_table_row_and_is_flagged_no_term_by_design() {
        // chembl_tier_b_0080's real shape: Ar-N=[N+]=[N-], central k = type
        // 53 (=N=). Oracle-confirmed None as above.
        assert!(mmff94_torsion_energy(0, 37, 9, 53, 47).is_none());
        assert!(torsion_no_term_by_design(9, 53));
    }

    #[test]
    fn caffeine_dione_ring_torsion_resolves_end_to_end_via_the_reperceived_view() {
        // Issue #227 Phase 1's actual root-cause fix, exercised through the
        // real production call shape (assign_mmff94_numeric_types_with_view
        // + torsion_type_for + mmff94_torsion_energy on the SAME reperceived
        // molecule), not just the diagnostic tool. Expected value is the
        // live RDKit oracle's own GetMMFFTorsionParams(mol,4,5,6,7) result
        // for this exact molecule/atom-quadruple (rdkit==2026.03.4):
        // (torsionType=1, 0.0, 2.5, 0.0).
        let m = chematic_smiles::parse("Cn1cnc2c1c(=O)n(C)c(=O)n2C").unwrap();
        let (types, view) = assign_mmff94_numeric_types_with_view(&m).unwrap();
        let (i, j, k, l) = (4, 5, 6, 7);
        let tt = torsion_type_for(&view, i, j, k, l, types[i], types[j], types[k], types[l]);
        assert_eq!(tt, 1, "must classify as type 1, not type 0");
        let params = mmff94_torsion_energy(tt, types[i], types[j], types[k], types[l])
            .expect("a real table row must resolve once fed the correct bond order");
        assert!((params.v1 - 0.0).abs() < 1e-9);
        assert!((params.v2 - 2.5).abs() < 1e-9, "v2={}", params.v2);
        assert!((params.v3 - 0.0).abs() < 1e-9);
    }

    #[test]
    fn torsion_lookup_is_symmetric_under_reversal_for_the_caffeine_case() {
        let m = chematic_smiles::parse("Cn1cnc2c1c(=O)n(C)c(=O)n2C").unwrap();
        let (types, view) = assign_mmff94_numeric_types_with_view(&m).unwrap();
        let (i, j, k, l) = (4usize, 5usize, 6usize, 7usize);
        let fwd = torsion_type_for(&view, i, j, k, l, types[i], types[j], types[k], types[l]);
        let rev = torsion_type_for(&view, l, k, j, i, types[l], types[k], types[j], types[i]);
        let p_fwd = mmff94_torsion_energy(fwd, types[i], types[j], types[k], types[l]).unwrap();
        let p_rev = mmff94_torsion_energy(rev, types[l], types[k], types[j], types[i]).unwrap();
        assert!((p_fwd.v1 - p_rev.v1).abs() < 1e-12);
        assert!((p_fwd.v2 - p_rev.v2).abs() < 1e-12);
        assert!((p_fwd.v3 - p_rev.v3).abs() < 1e-12);
    }

    #[test]
    fn cyclopentane_ring_torsion_uses_distinct_ring_specific_parameters() {
        // Confirms the classification fix is not inert: the real MMFF94
        // table has a dedicated (5, 1, 1, 1, 1) row distinct from the
        // generic (0, 0, 1, 1, 0) wildcard row old code would have used.
        let generic = mmff94_torsion_energy(0, 1, 1, 1, 1).expect("generic type-0 row");
        let ring5 = mmff94_torsion_energy(5, 1, 1, 1, 1).expect("ring type-5 row");
        assert!(
            (generic.v1, generic.v2, generic.v3) != (ring5.v1, ring5.v2, ring5.v3),
            "ring-specific torsion params must differ from the generic fallback"
        );
    }

    // ── FF-2: angle_type_for (#174) ─────────────────────────────────────────

    #[test]
    fn angle_type_for_cyclopropane_ring_angle_is_type_3() {
        let mol = chematic_smiles::parse("C1CC1").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let ring = first_ring_of_size(&mol, 3).expect("cyclopropane must have a 3-ring");
        let (i, j, k) = (ring[0], ring[1], ring[2]);
        let rings = find_sssr(&mol);
        let at = angle_type_for(&mol, rings.rings(), i, j, k, &types);
        assert_eq!(
            at, 3,
            "3-ring angle with two sp3 (bond_type 0) flanking bonds must be type 3"
        );
        // The specific (3,1,1,1) row doesn't exist for this crate's current
        // atom typer (which doesn't yet distinguish ring-size-specific sp3
        // carbon types 20/22 from plain type 1) — the type-0 fallback added
        // to `mmff94_angle_energy` must still return a usable value rather
        // than silently dropping the angle term.
        assert!(
            mmff94_angle_energy(at, types[i], types[j], types[k]).is_some(),
            "angle energy lookup must fall back to type 0, not silently miss"
        );
    }

    #[test]
    fn angle_type_for_cyclobutane_ring_angle_is_type_4() {
        let mol = chematic_smiles::parse("C1CCC1").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let ring = first_ring_of_size(&mol, 4).expect("cyclobutane must have a 4-ring");
        let (i, j, k) = (ring[0], ring[1], ring[2]);
        let rings = find_sssr(&mol);
        let at = angle_type_for(&mol, rings.rings(), i, j, k, &types);
        assert_eq!(
            at, 4,
            "4-ring angle with two sp3 (bond_type 0) flanking bonds must be type 4"
        );
    }

    // ── FF-2b: angle_type_for's corrected bt_sum=1/2 ring-offset formula
    // (issue #227 Priority 2C) ─────────────────────────────────────────────
    // RDKit's real `getMMFFAngleType`: 3-ring bt_sum=2 -> 6 (was wrongly 8);
    // 4-ring bt_sum=1 -> 7 (was wrongly 6), bt_sum=2 -> 8 (was wrongly 7).
    // Measured LATENT on the 265-molecule corpus (0/113 reachable), so these
    // fixtures are constructed molecules, not corpus-derived.

    #[test]
    fn angle_type_for_radialene_cyclopropane_bt_sum_2_ring_angle_is_type_6_not_8() {
        // [3]radialene (trimethylenecyclopropane): all three cyclopropane
        // ring carbons carry an exocyclic =CH2, so every ring bond is a
        // formal single bond between two MLTB (sp2, type 2) atoms ->
        // bond_type 1 on both flanking bonds of any ring angle -> bt_sum=2.
        let mol = chematic_smiles::parse("C1(=C)C(=C)C1=C").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let ring = first_ring_of_size(&mol, 3).expect("radialene core must have a 3-ring");
        let (i, j, k) = (ring[0], ring[1], ring[2]);
        let rings = find_sssr(&mol);
        let bt_ij = bond_type_for(types[i], types[j], bond_order_between(&mol, i, j));
        let bt_jk = bond_type_for(types[j], types[k], bond_order_between(&mol, j, k));
        assert_eq!(
            bt_ij + bt_jk,
            2,
            "both ring-flanking bonds must be bond_type 1"
        );
        let at = angle_type_for(&mol, rings.rings(), i, j, k, &types);
        assert_eq!(
            at, 6,
            "3-ring angle with bt_sum=2 must be type 6 (RDKit formula), not the old wrong 8"
        );
    }

    /// Real (molecule-adjacent, not just ring-list-adjacent) neighbors of
    /// `j` that also belong to `ring` -- avoids treating a ring's diagonal
    /// (non-bonded) atom pairs as flanking bonds.
    fn ring_bonded_neighbors(mol: &Molecule, ring: &[usize], j: usize) -> Vec<usize> {
        mol.neighbors(AtomIdx(j as u32))
            .map(|(n, _)| n.0 as usize)
            .filter(|n| ring.contains(n))
            .collect()
    }

    #[test]
    fn angle_type_for_4_ring_bt_sum_1_ring_angle_is_type_7_not_6() {
        // Two adjacent ring carbons carry exocyclic =CH2 (MLTB), the other
        // two are plain sp3 CH2 -- the ring angle centered on one
        // MLTB-substituted carbon has exactly one MLTB-MLTB flanking bond
        // (bond_type 1) and one MLTB-sp3 flanking bond (bond_type 0).
        let mol = chematic_smiles::parse("C1(=C)C(=C)CC1").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let ring = first_ring_of_size(&mol, 4).expect("must have a 4-ring");
        let rings = find_sssr(&mol);
        let mut found = None;
        for &j in &ring {
            let ring_neighbors = ring_bonded_neighbors(&mol, &ring, j);
            assert_eq!(
                ring_neighbors.len(),
                2,
                "each 4-ring atom has exactly 2 ring neighbors"
            );
            let (i, k) = (ring_neighbors[0], ring_neighbors[1]);
            if !atoms_share_ring_of_size(rings.rings(), &[i, j, k], 4) {
                continue;
            }
            let bt_ij = bond_type_for(types[i], types[j], bond_order_between(&mol, i, j));
            let bt_jk = bond_type_for(types[j], types[k], bond_order_between(&mol, j, k));
            if bt_ij + bt_jk == 1 {
                found = Some((i, j, k));
                break;
            }
        }
        let (i, j, k) = found.expect("must find a bt_sum=1 ring angle in this fixture");
        let at = angle_type_for(&mol, rings.rings(), i, j, k, &types);
        assert_eq!(
            at, 7,
            "4-ring angle with bt_sum=1 must be type 7 (RDKit formula), not the old wrong 6"
        );
    }

    #[test]
    fn angle_type_for_4_ring_bt_sum_2_ring_angle_is_type_8_not_7() {
        // Three consecutive ring carbons all carry exocyclic =CH2 (MLTB);
        // the angle centered on the middle one has both flanking ring bonds
        // MLTB-MLTB (bond_type 1 each) -> bt_sum=2.
        let mol = chematic_smiles::parse("C1(=C)C(=C)C(=C)C1").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let ring = first_ring_of_size(&mol, 4).expect("must have a 4-ring");
        let rings = find_sssr(&mol);
        let mut found = None;
        for &j in &ring {
            let ring_neighbors = ring_bonded_neighbors(&mol, &ring, j);
            assert_eq!(
                ring_neighbors.len(),
                2,
                "each 4-ring atom has exactly 2 ring neighbors"
            );
            let (i, k) = (ring_neighbors[0], ring_neighbors[1]);
            if !atoms_share_ring_of_size(rings.rings(), &[i, j, k], 4) {
                continue;
            }
            let bt_ij = bond_type_for(types[i], types[j], bond_order_between(&mol, i, j));
            let bt_jk = bond_type_for(types[j], types[k], bond_order_between(&mol, j, k));
            if bt_ij + bt_jk == 2 {
                found = Some((i, j, k));
                break;
            }
        }
        let (i, j, k) = found.expect("must find a bt_sum=2 ring angle in this fixture");
        let at = angle_type_for(&mol, rings.rings(), i, j, k, &types);
        assert_eq!(
            at, 8,
            "4-ring angle with bt_sum=2 must be type 8 (RDKit formula), not the old wrong 7"
        );
    }

    // ── FF-3: stretch_bend_type_for (issue #227 Priority 2C) ────────────────
    // Reproduces the diagnostic's own self-consistency proof
    // (`mmff94_stbn_equivalence_diagnostic_227.rs`'s top doc comment,
    // ~lines 38-56): `MMFF94_STBN`'s frozen data has exactly one row keyed
    // 5 -- an all-CR3R (type 22, cyclopropane ring carbon) triple -- and its
    // 11 key-4 rows are all CR4R (type 20) triples. A 3-ring/bt_sum=0 angle
    // is angle_type 3, never 5 under `angle_type_for`'s own table, so if the
    // STBN key column really were `angle_type`, the key-5 row would be
    // unreachable garbage. `getMMFFStretchBendType` resolves this: angle
    // type 3 -> stretch-bend type 5 (not 3), and angle type 4 -> stretch-
    // bend type 4 (the one case where the numeral is unchanged).

    #[test]
    fn stretch_bend_type_for_cyclopropane_all_cr3r_is_type_5_not_3() {
        assert_eq!(
            stretch_bend_type_for(3, 22, 22, 0, 0),
            5,
            "angle_type 3 (3-ring, bt_sum=0) must map to stretch-bend type 5"
        );
        // End-to-end: the resulting key 5 must actually resolve the real
        // MMFF94_STBN row for the CR3R triple -- confirming the column
        // really is stretch_bend_type, not angle_type.
        let sbt = stretch_bend_type_for(3, 22, 22, 0, 0);
        assert_eq!(
            mmff94_stbn_type_only(sbt, 22, 22, 22),
            Some((0.0, 0.0)),
            "stretch-bend type 5 must resolve the (5,22,22,22,0.0,0.0) MMFF94_STBN row"
        );
    }

    #[test]
    fn stretch_bend_type_for_cr4r_4_ring_is_type_4() {
        assert_eq!(
            stretch_bend_type_for(4, 20, 20, 0, 0),
            4,
            "angle_type 4 (4-ring, bt_sum=0) must map to stretch-bend type 4 (unchanged numeral)"
        );
    }

    #[test]
    fn angle_type_for_butadiene_sp2_single_bond_is_type_1_and_finds_dedicated_row() {
        // Concrete non-ring demonstration that the fix is not inert: the
        // flanking C1-C2 single bond between two vinylic (type 2) carbons is
        // a real MMFF94 bond_type=1 (sbmb) bond, giving angle_type=1, which
        // has a dedicated (1,2,2,2) row (theta0=121.55) distinct from the
        // old hardcoded angle_type=0's row for the same atom-type triple
        // (theta0=118.043, reachable at angle_type=0 via issue #227 Stage
        // B's eqLevel equivalence ladder -- type 2's own eqLevel table
        // substitutes it to type 1 at Level 4, and a real (0,1,2,1) row
        // exists). Getting `angle_type` wrong here is a silently-DIFFERENT
        // wrong value, not a clean miss -- the same "silent wrong parameter"
        // failure class as the #236 furan collision, which is exactly why
        // `angle_type_for`'s own correctness (asserted below) matters.
        let mol = chematic_smiles::parse("C=CC=C").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        assert_eq!(types[0], 2);
        assert_eq!(types[1], 2);
        assert_eq!(types[2], 2);
        let rings = find_sssr(&mol);
        let at = angle_type_for(&mol, rings.rings(), 0, 1, 2, &types);
        assert_eq!(at, 1, "C0=C1-C2 angle must classify as angle_type 1");
        let wrong = mmff94_angle_energy(0, 2, 2, 2).expect("angle_type=0 now resolves too");
        let p = mmff94_angle_energy(at, 2, 2, 2).expect("(1,2,2,2) row must exist");
        assert!(
            (p.theta0 - 121.55).abs() < 1e-6,
            "theta0={} should be 121.55",
            p.theta0
        );
        assert!(
            (p.theta0 - wrong.theta0).abs() > 1.0,
            "angle_type=0's wrong theta0={} should differ meaningfully from the correct {}",
            wrong.theta0,
            p.theta0
        );
    }

    #[test]
    fn butadiene_angle_energy_minimum_is_near_correct_theta0() {
        // Energy-perturbation test: with the fix, this angle now has a real
        // restoring force toward its correct 121.55° equilibrium (pre-fix it
        // silently contributed zero energy everywhere — no minimum at all).
        let mol = chematic_smiles::parse("C=CC=C").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let rings = find_sssr(&mol);
        let scan_energy = |theta_deg: f64| {
            let r = 1.45_f64;
            let half = theta_deg.to_radians() / 2.0;
            // j (atom1) at origin; i (atom0) and k (atom2) placed to form
            // the target angle at j.
            let coords = vec![
                [r * half.dcos(), r * half.dsin(), 0.0],
                [0.0, 0.0, 0.0],
                [r * half.dcos(), -r * half.dsin(), 0.0],
                [r * half.dcos() + 1.3, -r * half.dsin() - 0.9, 0.0],
            ];
            angle_energy(&mol, &coords, &types, rings.rings())
        };
        let e_correct = scan_energy(121.55);
        let e_distorted = scan_energy(100.0);
        assert!(
            e_distorted > e_correct,
            "energy away from the true 121.55° minimum ({e_distorted}) should exceed energy \
             at the minimum ({e_correct})"
        );
    }

    #[test]
    fn angle_gradient_restores_butadiene_toward_theta0() {
        let mol = chematic_smiles::parse("C=CC=C").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let charges = mmff94_charges_numeric(&mol).unwrap_or_else(|_| vec![0.0; mol.atom_count()]);
        let rings = find_sssr(&mol);
        let r = 1.45_f64;
        let half = 100.0_f64.to_radians() / 2.0; // distorted away from 121.55°
        let coords = vec![
            [r * half.dcos(), r * half.dsin(), 0.0],
            [0.0, 0.0, 0.0],
            [r * half.dcos(), -r * half.dsin(), 0.0],
            [r * half.dcos() + 1.3, -r * half.dsin() - 0.9, 0.0],
        ];
        let grad = compute_gradient(&mol, &coords, &types, &charges, rings.rings(), 1e-4);
        let grad_norm: f64 = grad
            .iter()
            .flat_map(|g| g.iter())
            .map(|x| x * x)
            .sum::<f64>()
            .sqrt();
        assert!(
            grad_norm > 1e-6,
            "gradient should be nonzero away from equilibrium angle"
        );
    }

    #[test]
    fn dihedral_anti_is_pi() {
        let i = [0.0_f64, 0.0, 0.0];
        let j = [1.0, 0.0, 0.0];
        let k = [2.0, 0.0, 1.0];
        let l = [3.0, 0.0, 0.0];
        let phi = dihedral(i, j, k, l);
        assert!(phi.abs() > 2.5, "anti dihedral ≈ π: {}", phi);
    }

    #[test]
    fn dihedral_syn_is_zero() {
        let i = [0.0_f64, 1.0, 0.0];
        let j = [0.0, 0.0, 0.0];
        let k = [1.0, 0.0, 0.0];
        let l = [1.0, 1.0, 0.0];
        let phi = dihedral(i, j, k, l);
        assert!(phi.abs() < 0.1, "syn dihedral ≈ 0: {}", phi);
    }

    #[test]
    fn minimize_reduces_energy_for_methane() {
        let (mol, _) = methane_mol();
        let mut coords = vec![
            [0.0, 0.0, 0.0_f64],
            [1.5, 1.5, 1.5],
            [-1.5, -1.5, 1.5],
            [-1.5, 1.5, -1.5],
            [1.5, -1.5, -1.5],
        ];
        let e_before = mmff94_total_energy(&mol, &coords).expect("energy before");
        let result = minimize_mmff94_full(&mol, &mut coords, 300).expect("minimize");
        assert!(
            result.energy <= e_before,
            "minimize should reduce energy: {} → {}",
            e_before,
            result.energy
        );
        assert!(result.energy.is_finite());
        assert!(result.iterations > 0);
    }

    #[test]
    fn lbfgs_reduces_energy_for_methane() {
        let (mol, _) = methane_mol();
        let mut coords = vec![
            [0.0, 0.0, 0.0_f64],
            [1.5, 1.5, 1.5],
            [-1.5, -1.5, 1.5],
            [-1.5, 1.5, -1.5],
            [1.5, -1.5, -1.5],
        ];
        let e_before = mmff94_total_energy(&mol, &coords).expect("energy before");
        let result = minimize_mmff94_lbfgs(&mol, &mut coords, 300).expect("lbfgs");
        assert!(
            result.energy <= e_before,
            "L-BFGS should reduce energy: {} → {}",
            e_before,
            result.energy
        );
        assert!(result.energy.is_finite());
    }

    #[test]
    fn bounded_analytic_lbfgs_reduces_energy_for_butane() {
        let mol = chematic_smiles::parse("CCCC").expect("butane");
        let model = Mmff94EnergyModel::new(&mol).expect("MMFF94 preparation");
        let mut coords = vec![
            [0.0, 0.0, 0.0],
            [1.54, 0.1, 0.0],
            [2.95, 0.85, 0.2],
            [4.25, 0.15, 1.05],
            [0.0, 1.05, 0.8],
            [1.55, -0.9, -0.75],
            [1.55, 0.25, 1.05],
            [2.95, 1.9, -0.55],
            [2.95, 1.05, 1.25],
            [4.2, -0.8, 1.45],
            [4.2, 0.0, 2.05],
            [4.9, 0.35, 0.35],
            [4.0, 0.95, 0.05],
            [4.6, 0.35, 1.55],
        ];
        let before = model.energy(&coords);
        let result = model
            .minimize_lbfgs_bounded_analytic(&mut coords, 100)
            .expect("bounded analytic L-BFGS");
        assert!(result.energy.is_finite());
        assert!(
            result.energy <= before,
            "analytic L-BFGS should reduce energy: {before} -> {}",
            result.energy
        );
    }

    #[test]
    fn bounded_analytic_constraint_never_commits_a_rejected_step() {
        let mol = chematic_smiles::parse("CCCC").expect("butane");
        let model = Mmff94EnergyModel::new(&mol).expect("MMFF94 preparation");
        let mut coords = vec![
            [0.0, 0.0, 0.0],
            [1.54, 0.1, 0.0],
            [2.95, 0.85, 0.2],
            [4.25, 0.15, 1.05],
            [0.0, 1.05, 0.8],
            [1.55, -0.9, -0.75],
            [1.55, 0.25, 1.05],
            [2.95, 1.9, -0.55],
            [2.95, 1.05, 1.25],
            [4.2, -0.8, 1.45],
            [4.2, 0.0, 2.05],
            [4.9, 0.35, 0.35],
            [4.0, 0.95, 0.05],
            [4.6, 0.35, 1.55],
        ];
        let initial = coords.clone();
        let result = model
            .minimize_lbfgs_bounded_analytic_with_constraint(&mut coords, 100, |_| false)
            .expect("constraint rejection is ordinary non-convergence");
        assert_eq!(coords, initial);
        assert!(!result.converged);
        assert_eq!(result.iterations, 1);
        assert_eq!(
            result.termination,
            Mmff94TerminationReason::ConstraintRejectedFallback
        );
    }

    #[test]
    fn bounded_analytic_reports_iteration_limit_separately() {
        let mol = chematic_smiles::parse("CCCC").expect("butane");
        let model = Mmff94EnergyModel::new(&mol).expect("MMFF94 preparation");
        let mut coords = vec![
            [0.0, 0.0, 0.0],
            [1.54, 0.1, 0.0],
            [2.95, 0.85, 0.2],
            [4.25, 0.15, 1.05],
        ];
        let before = coords.clone();
        let result = model
            .minimize_lbfgs_bounded_analytic(&mut coords, 0)
            .expect("zero-budget minimization");
        assert_eq!(coords, before);
        assert!(!result.converged);
        assert_eq!(result.iterations, 0);
        assert_eq!(result.termination, Mmff94TerminationReason::IterationLimit);
    }

    #[test]
    fn lbfgs_converges_in_fewer_iters_than_sd() {
        let (mol, _) = methane_mol();
        // Moderately distorted — both should converge but L-BFGS faster
        let base_coords = vec![
            [0.0, 0.0, 0.0_f64],
            [1.2, 1.2, 1.2],
            [-1.2, -1.2, 1.2],
            [-1.2, 1.2, -1.2],
            [1.2, -1.2, -1.2],
        ];
        let mut coords_sd = base_coords.clone();
        let mut coords_lbfgs = base_coords;
        let sd = minimize_mmff94_full(&mol, &mut coords_sd, 500).expect("sd");
        let lb = minimize_mmff94_lbfgs(&mol, &mut coords_lbfgs, 500).expect("lbfgs");
        // Both should converge; L-BFGS should need ≤ SD iterations
        assert!(
            lb.iterations <= sd.iterations || lb.converged,
            "L-BFGS iters={} SD iters={}",
            lb.iterations,
            sd.iterations
        );
        assert!(lb.energy.is_finite());
    }

    #[test]
    fn energy_breakdown_sums_to_total() {
        let (mol, coords) = methane_mol();
        let bd = mmff94_energy_breakdown(&mol, &coords).expect("breakdown");
        let sum =
            bd.bond + bd.angle + bd.stretch_bend + bd.torsion + bd.oop + bd.vdw + bd.electrostatic;
        assert!(
            (sum - bd.total).abs() < 1e-10,
            "sum={} total={}",
            sum,
            bd.total
        );
        assert!(bd.total.is_finite());
    }

    /// Frozen 58-molecule corpus, copied from `scripts/etkdg_vs_rdkit_gap.py::CORPUS`
    /// (the same corpus PR #169's `mmff94_bridge_coverage_report` example measured
    /// coverage against). Used here to count MMFF94 parameter-lookup *misses*
    /// per term class (bond/angle/torsion) before vs. after the classification
    /// fixes in this file — the metric that actually matters, since a
    /// "corrected" classification that routes to a still-missing row is a
    /// silent zero-energy regression, not an improvement.
    const CORPUS_58: &[(&str, &str)] = &[
        ("benzene", "c1ccccc1"),
        ("naphthalene", "c1ccc2ccccc2c1"),
        ("pyridine", "c1ccncc1"),
        ("furan", "c1ccoc1"),
        ("thiophene", "c1ccsc1"),
        ("adamantane", "C1CC2CC3CC1CC(C2)C3"),
        ("cubane", "C1C2C3C1C4C2C3C4"),
        ("cyclohexane", "C1CCCCC1"),
        ("cyclopentane", "C1CCCC1"),
        ("indole", "c1ccc2[nH]ccc2c1"),
        ("purine", "c1ncc2[nH]cnc2n1"),
        ("quinoline", "c1ccc2ncccc2c1"),
        ("anthracene", "c1ccc2cc3ccccc3cc2c1"),
        ("pyrene", "c1cc2ccc3cccc4ccc(c1)c2c34"),
        ("biphenyl", "c1ccc(-c2ccccc2)cc1"),
        ("butane", "CCCC"),
        ("hexane", "CCCCCC"),
        ("decane", "CCCCCCCCCC"),
        ("triethylene_glycol", "OCCOCCOCCO"),
        ("hexanediol", "OCCCCCCO"),
        ("hexadecane", "CCCCCCCCCCCCCCCC"),
        ("cyclododecane", "C1CCCCCCCCCCC1"),
        ("crown_12_4", "O1CCOCCOCCOCC1"),
        ("cyclooctadecane", "C1CCCCCCCCCCCCCCCCC1"),
        ("l_alanine", "N[C@@H](C)C(=O)O"),
        ("d_alanine", "N[C@H](C)C(=O)O"),
        ("l_serine", "N[C@@H](CO)C(=O)O"),
        ("l_threonine", "C[C@H](O)[C@@H](N)C(=O)O"),
        ("2_butanol_R", "C[C@H](O)CC"),
        ("2_butanol_S", "C[C@@H](O)CC"),
        ("2_chlorobutane_R", "C[C@H](Cl)CC"),
        ("ibuprofen_S", "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O"),
        ("naproxen_S", "COc1ccc2cc([C@H](C)C(=O)O)ccc2c1"),
        ("menthol", "C[C@@H]1CC[C@@H](C(C)C)C[C@H]1O"),
        ("chfclbr_R", "[C@H](F)(Cl)Br"),
        ("chfclbr_S", "[C@@H](F)(Cl)Br"),
        ("quaternary_1_R", "[C@](F)(Cl)(Br)I"),
        ("quaternary_1_S", "[C@@](F)(Cl)(Br)I"),
        ("quaternary_2_R", "[C@](C)(N)(O)F"),
        ("quaternary_2_S", "[C@@](C)(N)(O)F"),
        ("but2ene_E", "C/C=C/C"),
        ("but2ene_Z", r"C/C=C\C"),
        ("chloropropene_E", "C(/C=C/C)Cl"),
        ("chloropropene_Z", r"C(/C=C\C)Cl"),
        ("cinnamic_acid_E", "OC(=O)/C=C/c1ccccc1"),
        ("cinnamic_acid_Z", r"OC(=O)/C=C\c1ccccc1"),
        ("pent2ene_E", "CC/C=C/C"),
        ("pent2ene_Z", r"CC/C=C\C"),
        ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
        ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
        ("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
        ("diphenhydramine", "CN(C)CCOC(c1ccccc1)c1ccccc1"),
        (
            "penicillin_core",
            "CC1(C)S[C@@H]2[C@H](NC(=O)C)C(=O)N2[C@H]1C(=O)O",
        ),
        (
            "testosterone",
            "C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O",
        ),
        (
            "cholesterol",
            "C[C@H](CCCC(C)C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
        ),
        (
            "atorvastatin_fragment",
            "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",
        ),
        ("gly_ala_gly", "NCC(=O)N[C@@H](C)C(=O)NCC(=O)O"),
    ];

    /// Per-term-class MMFF94 parameter-lookup miss counts across [`CORPUS_58`].
    struct CoverageCounts {
        bond_total: usize,
        bond_miss: usize,
        angle_total: usize,
        angle_miss: usize,
        torsion_total: usize,
        torsion_miss: usize,
    }

    fn corpus_coverage() -> CoverageCounts {
        let mut c = CoverageCounts {
            bond_total: 0,
            bond_miss: 0,
            angle_total: 0,
            angle_miss: 0,
            torsion_total: 0,
            torsion_miss: 0,
        };
        for (name, smiles) in CORPUS_58 {
            let mol = chematic_smiles::parse(smiles).unwrap_or_else(|e| {
                panic!("corpus molecule {name} ({smiles}) failed to parse: {e}")
            });
            let types = match assign_mmff94_numeric_types(&mol) {
                Ok(t) => t,
                Err(_) => continue, // unsupported element/type — not this file's concern
            };
            let rings = find_sssr(&mol);

            for (_, bond) in mol.bonds() {
                let i = bond.atom1.0 as usize;
                let j = bond.atom2.0 as usize;
                c.bond_total += 1;
                let bt = bond_type_for(types[i], types[j], bond.order);
                if mmff94_bond_energy(bt, types[i], types[j]).is_none() {
                    c.bond_miss += 1;
                }
            }

            for j_idx in 0..mol.atom_count() {
                let j = AtomIdx(j_idx as u32);
                let neighbors: Vec<usize> = mol.neighbors(j).map(|(nb, _)| nb.0 as usize).collect();
                if neighbors.len() < 2 {
                    continue;
                }
                for (ii, &i) in neighbors.iter().enumerate() {
                    for &k in &neighbors[ii + 1..] {
                        c.angle_total += 1;
                        let at = angle_type_for(&mol, rings.rings(), i, j_idx, k, &types);
                        if mmff94_angle_energy(at, types[i], types[j_idx], types[k]).is_none() {
                            c.angle_miss += 1;
                        }
                    }
                }
            }

            for (_, bond) in mol.bonds() {
                let j = bond.atom1.0 as usize;
                let k = bond.atom2.0 as usize;
                let nbrs_j: Vec<usize> = mol
                    .neighbors(bond.atom1)
                    .map(|(nb, _)| nb.0 as usize)
                    .collect();
                let nbrs_k: Vec<usize> = mol
                    .neighbors(bond.atom2)
                    .map(|(nb, _)| nb.0 as usize)
                    .collect();
                for &i in &nbrs_j {
                    if i == k {
                        continue;
                    }
                    for &l in &nbrs_k {
                        if l == j {
                            continue;
                        }
                        c.torsion_total += 1;
                        let tt = torsion_type_for(
                            &mol, i, j, k, l, types[i], types[j], types[k], types[l],
                        );
                        if mmff94_torsion_energy(tt, types[i], types[j], types[k], types[l])
                            .is_none()
                        {
                            c.torsion_miss += 1;
                        }
                    }
                }
            }
        }
        c
    }

    #[test]
    fn corpus_58_parameter_coverage_does_not_regress() {
        // Numbers on THIS 58-molecule corpus (CORPUS_58, distinct from the
        // 265-molecule Wave 1 corpus used elsewhere in issue #227) have
        // drifted across several since-merged fixes; re-measured fresh
        // rather than trusted from stale history:
        //   original bond_type_for/angle_type_for/torsion_type_for:
        //     bond 129/585 miss, angle 292/737 miss, torsion 305/841 miss
        //   immediately pre-this-fix (post atom-type-parity + stretch-bend
        //   fixes, i.e. this branch's parent commit):
        //     bond   1/585 miss, angle  24/737 miss, torsion  22/841 miss
        //   post this fix (issue #227 torsion classification: bond-type-jk-
        //   based base case + local ring-4/5 override, replacing the old
        //   atom-type-membership + SSSR-based formula):
        //     bond   1/585 miss, angle  24/737 miss, torsion   4/841 miss
        // Torsion coverage DOES improve now (unlike the original bond/angle-
        // only fix this comment used to describe): the corrected
        // classification routes several of this corpus's torsions to a
        // table row chematic's own existing fallback chain can reach that
        // the old, wrong classification code could not.
        let c = corpus_coverage();
        eprintln!(
            "bond: {}/{} miss, angle: {}/{} miss, torsion: {}/{} miss",
            c.bond_miss, c.bond_total, c.angle_miss, c.angle_total, c.torsion_miss, c.torsion_total
        );
        // Loose tripwires (see MEMORY.md's note on avoiding over-tuned
        // pseudo-gates), not tight equality: a future change may legitimately
        // shift these by a few rows (e.g. an atom-typer improvement), but
        // must never regress back toward the pre-fix counts above.
        assert!(
            c.bond_miss <= 45,
            "bond misses regressed toward pre-fix (129): {}",
            c.bond_miss
        );
        assert!(
            c.angle_miss <= 290,
            "angle misses regressed toward pre-fix (292): {}",
            c.angle_miss
        );
        assert!(
            c.torsion_miss <= 20,
            "torsion misses regressed toward pre-this-fix (22): {}",
            c.torsion_miss
        );
    }

    #[test]
    fn energy_breakdown_bond_term_positive_for_distorted() {
        let (mol, _) = methane_mol();
        // Very distorted C-H bonds → high bond energy
        let stretched = vec![
            [0.0, 0.0, 0.0_f64],
            [2.0, 2.0, 2.0],
            [-2.0, -2.0, 2.0],
            [-2.0, 2.0, -2.0],
            [2.0, -2.0, -2.0],
        ];
        let bd = mmff94_energy_breakdown(&mol, &stretched).expect("breakdown");
        assert!(
            bd.bond > 0.0,
            "stretched bond energy should be positive: {}",
            bd.bond
        );
    }

    #[test]
    fn prepared_nonbonded_gradient_matches_finite_difference() {
        let mol = chematic_smiles::parse("CCC").expect("propane");
        let model = Mmff94EnergyModel::new(&mol).expect("MMFF94 preparation");
        let coords = vec![
            [0.0, 0.0, 0.0],
            [1.54, 0.1, 0.0],
            [3.08, 0.0, 0.2],
            [-0.4, 0.9, 0.8],
            [-0.4, -0.9, -0.8],
            [1.54, 1.0, 0.9],
            [1.54, -1.0, -0.9],
            [3.48, 0.8, 0.8],
            [3.48, -0.8, -0.8],
            [3.08, 0.0, 1.2],
            [3.08, 0.0, -0.8],
        ];
        let actual = model.nonbonded_gradient(&coords);
        let vdw_only = prepared_nonbonded_gradient(&coords, &model.vdw_pairs, &[]);
        let electrostatic_only =
            prepared_nonbonded_gradient(&coords, &[], &model.electrostatic_pairs);
        let delta = 1e-5;
        let mut work = coords.clone();
        for (i, row) in actual.iter().enumerate() {
            for (axis, &value) in row.iter().enumerate() {
                work[i][axis] += delta;
                let ep = model.vdw_energy(&work) + model.electrostatic_energy(&work);
                work[i][axis] -= 2.0 * delta;
                let em = model.vdw_energy(&work) + model.electrostatic_energy(&work);
                work[i][axis] += delta;
                let expected = (ep - em) / (2.0 * delta);
                assert!((value - expected).abs() < 1e-6 * (1.0 + expected.abs()));

                let mut plus = coords.clone();
                let mut minus = coords.clone();
                plus[i][axis] += delta;
                minus[i][axis] -= delta;
                let expected_vdw =
                    (model.vdw_energy(&plus) - model.vdw_energy(&minus)) / (2.0 * delta);
                let expected_electrostatic = (model.electrostatic_energy(&plus)
                    - model.electrostatic_energy(&minus))
                    / (2.0 * delta);
                assert!(
                    (vdw_only[i][axis] - expected_vdw).abs() < 1e-6 * (1.0 + expected_vdw.abs()),
                    "vdW gradient mismatch at atom {i}, axis {axis}: actual={}, expected={expected_vdw}",
                    vdw_only[i][axis]
                );
                assert!(
                    (electrostatic_only[i][axis] - expected_electrostatic).abs()
                        < 1e-6 * (1.0 + expected_electrostatic.abs()),
                    "electrostatic gradient mismatch at atom {i}, axis {axis}: actual={}, expected={expected_electrostatic}",
                    electrostatic_only[i][axis]
                );
            }
        }
    }

    #[test]
    fn coordinate_dependent_vdw_neighbor_list_matches_cutoff_reference() {
        let mol = chematic_smiles::parse("CCCC").expect("butane");
        let model = Mmff94EnergyModel::new(&mol).expect("MMFF94 preparation");
        let coords = vec![
            [0.0, 0.0, 0.0],
            [1.54, 0.1, 0.0],
            [2.95, 0.85, 0.2],
            [4.25, 0.15, 1.05],
            [0.0, 1.05, 0.8],
            [1.55, -0.9, -0.75],
            [1.55, 0.25, 1.05],
            [2.95, 1.9, -0.55],
            [2.95, 1.05, 1.25],
            [4.2, -0.8, 1.45],
            [4.2, 0.0, 2.05],
            [4.9, 0.35, 0.35],
            [4.0, 0.95, 0.05],
            [4.6, 0.35, 1.55],
        ];
        let active = prepared_vdw_neighbor_list(&coords, &model.vdw_pairs);
        let actual: std::collections::HashSet<_> =
            active.iter().map(|pair| (pair.i, pair.j)).collect();
        let expected: std::collections::HashSet<_> = model
            .vdw_pairs
            .iter()
            .filter(|pair| dist(coords[pair.i], coords[pair.j]) <= 10.0)
            .map(|pair| (pair.i, pair.j))
            .collect();
        assert_eq!(actual, expected);
        assert_eq!(
            vdw_energy_pairs(&coords, &model.vdw_pairs),
            model.vdw_pairs.iter().fold(0.0, |total, pair| {
                let r = dist(coords[pair.i], coords[pair.j]);
                if r > 10.0 || r <= 0.01 {
                    total
                } else {
                    total + mmff94_vdw_energy_value(r, pair.r_star, pair.epsilon)
                }
            })
        );
    }

    #[test]
    fn coordinate_dependent_vdw_neighbor_list_handles_boundaries_and_negative_coordinates() {
        let coords = vec![
            [-10.0, -0.0, 0.0],
            [-0.00000000001, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [9.99999999999, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [10.00000000001, 0.0, 0.0],
            [19.99999999999, 0.0, 0.0],
            [20.0, 0.0, 0.0],
            [-5.0, 8.0, -6.0],
        ];
        let pairs: Vec<_> = (0..coords.len())
            .flat_map(|i| {
                (i + 1..coords.len()).map(move |j| PreparedVdwPair {
                    i,
                    j,
                    r_star: 1.0,
                    epsilon: 1.0,
                })
            })
            .collect();
        let active: std::collections::HashSet<_> = prepared_vdw_neighbor_list(&coords, &pairs)
            .iter()
            .map(|pair| (pair.i, pair.j))
            .collect();
        let expected: std::collections::HashSet<_> = pairs
            .iter()
            .filter(|pair| dist(coords[pair.i], coords[pair.j]) <= 10.0)
            .map(|pair| (pair.i, pair.j))
            .collect();
        assert_eq!(active, expected);
    }

    #[test]
    fn coordinate_dependent_vdw_neighbor_list_falls_back_for_non_finite_coordinates() {
        let coords = vec![[0.0, 0.0, 0.0], [f64::NAN, 0.0, 0.0]];
        let pairs = vec![PreparedVdwPair {
            i: 0,
            j: 1,
            r_star: 1.0,
            epsilon: 1.0,
        }];
        let active = prepared_vdw_neighbor_list(&coords, &pairs);
        assert_eq!(active.len(), pairs.len());
        assert_eq!((active[0].i, active[0].j), (0, 1));
    }

    #[test]
    fn coordinate_dependent_electrostatic_neighbor_list_matches_cutoff_reference() {
        let coords = vec![
            [-10.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [10.0 + 1e-10, 0.0, 0.0],
            [-5.0, 8.0, -6.0],
        ];
        let pairs: Vec<_> = (0..coords.len())
            .flat_map(|i| {
                (i + 1..coords.len()).map(move |j| PreparedElectrostaticPair {
                    i,
                    j,
                    scale_charge_product: 1.0,
                })
            })
            .collect();
        let active: std::collections::HashSet<_> =
            prepared_electrostatic_neighbor_list(&coords, &pairs)
                .iter()
                .map(|pair| (pair.i, pair.j))
                .collect();
        let expected: std::collections::HashSet<_> = pairs
            .iter()
            .filter(|pair| dist(coords[pair.i], coords[pair.j]) <= 10.0)
            .map(|pair| (pair.i, pair.j))
            .collect();
        assert_eq!(active, expected);
        let reference: f64 = pairs
            .iter()
            .filter(|pair| dist(coords[pair.i], coords[pair.j]) <= 10.0)
            .map(|pair| 332.0716 / (dist(coords[pair.i], coords[pair.j]) + 0.05))
            .sum();
        assert_eq!(
            electrostatic_energy_cutoff_pairs(&coords, &pairs),
            reference
        );
    }

    #[test]
    fn cutoff_nonbonded_gradient_matches_cutoff_energy_finite_difference() {
        let mol = chematic_smiles::parse("CCO").expect("ethanol");
        let model = Mmff94EnergyModel::new(&mol).expect("MMFF94 preparation");
        let coords = vec![
            [0.0, 0.0, 0.0],
            [1.52, 0.1, 0.0],
            [2.9, 0.0, 0.2],
            [-0.4, 0.9, 0.8],
            [-0.4, -0.9, -0.8],
            [1.52, 1.0, 0.9],
            [1.52, -1.0, -0.9],
            [2.9, 0.8, 0.9],
            [2.9, -0.8, -0.7],
        ];
        let actual = model.cutoff_nonbonded_gradient(&coords);
        let delta = 1e-5;
        for (i, row) in actual.iter().enumerate() {
            for (axis, &value) in row.iter().enumerate() {
                let mut plus = coords.clone();
                let mut minus = coords.clone();
                plus[i][axis] += delta;
                minus[i][axis] -= delta;
                let ep = model.cutoff_nonbonded_energy(&plus);
                let em = model.cutoff_nonbonded_energy(&minus);
                let expected = ((ep.0 + ep.1) - (em.0 + em.1)) / (2.0 * delta);
                assert!(
                    (value - expected).abs() < 1e-5 * (1.0 + expected.abs()),
                    "cutoff nonbonded gradient mismatch at atom {i}, axis {axis}: actual={value}, expected={expected}"
                );
            }
        }
    }

    #[test]
    fn prepared_bond_angle_gradient_matches_finite_difference() {
        let (mol, coords) = methane_mol();
        let model = Mmff94EnergyModel::new(&mol).expect("MMFF94 preparation");
        let actual = model.bond_angle_gradient(&coords);
        let delta = 1e-5;
        let mut work = coords.clone();
        for (i, row) in actual.iter().enumerate() {
            for (axis, &value) in row.iter().enumerate() {
                work[i][axis] += delta;
                let ep = model.bond_energy(&work) + model.angle_energy(&work);
                work[i][axis] -= 2.0 * delta;
                let em = model.bond_energy(&work) + model.angle_energy(&work);
                work[i][axis] += delta;
                let expected = (ep - em) / (2.0 * delta);
                assert!(
                    (value - expected).abs() < 2e-5 * (1.0 + expected.abs()),
                    "bond/angle gradient mismatch at atom {i}, axis {axis}: actual={value}, expected={expected}"
                );
            }
        }
    }

    #[test]
    fn prepared_stretch_bend_gradient_matches_finite_difference() {
        let mol = chematic_smiles::parse("CCC").expect("propane");
        let model = Mmff94EnergyModel::new(&mol).expect("MMFF94 preparation");
        let coords = vec![
            [0.0, 0.0, 0.0],
            [1.54, 0.1, 0.0],
            [3.08, 0.0, 0.2],
            [-0.4, 0.9, 0.8],
            [-0.4, -0.9, -0.8],
            [1.54, 1.0, 0.9],
            [1.54, -1.0, -0.9],
            [3.48, 0.8, 0.8],
            [3.48, -0.8, -0.8],
            [3.08, 0.0, 1.2],
            [3.08, 0.0, -0.8],
        ];
        let actual = model.stretch_bend_gradient(&coords);
        let delta = 1e-5;
        let mut nonzero = false;
        for (i, row) in actual.iter().enumerate() {
            for (axis, &value) in row.iter().enumerate() {
                let mut plus = coords.clone();
                let mut minus = coords.clone();
                plus[i][axis] += delta;
                minus[i][axis] -= delta;
                let expected = (model.stretch_bend_energy(&plus)
                    - model.stretch_bend_energy(&minus))
                    / (2.0 * delta);
                assert!(
                    (value - expected).abs() < 2e-5 * (1.0 + expected.abs()),
                    "stretch-bend gradient mismatch at atom {i}, axis {axis}: actual={value}, expected={expected}"
                );
                nonzero |= value.abs() > 1e-8;
            }
        }

        assert!(nonzero, "propane stretch-bend probe must exercise a term");
    }

    #[test]
    fn prepared_torsion_gradient_matches_finite_difference() {
        let mol = chematic_smiles::parse("CCCC").expect("butane");
        let model = Mmff94EnergyModel::new(&mol).expect("MMFF94 preparation");
        assert!(!model.torsions.is_empty(), "butane should prepare torsions");
        let coords = vec![
            [0.0, 0.0, 0.0],
            [1.54, 0.1, 0.0],
            [2.95, 0.85, 0.2],
            [4.25, 0.15, 1.05],
            [0.0, 1.05, 0.8],
            [1.55, -0.9, -0.75],
            [1.55, 0.25, 1.05],
            [2.95, 1.9, -0.55],
            [2.95, 1.05, 1.25],
            [4.2, -0.8, 1.45],
            [4.2, 0.0, 2.05],
            [4.9, 0.35, 0.35],
            [4.0, 0.95, 0.05],
            [4.6, 0.35, 1.55],
        ];
        let actual = model.torsion_gradient(&coords);
        let delta = 1e-5;
        let mut work = coords.clone();
        for (i, row) in actual.iter().enumerate() {
            for (axis, &value) in row.iter().enumerate() {
                work[i][axis] += delta;
                let ep = model.torsion_energy(&work);
                work[i][axis] -= 2.0 * delta;
                let em = model.torsion_energy(&work);
                work[i][axis] += delta;
                let expected = (ep - em) / (2.0 * delta);
                assert!(
                    (value - expected).abs() < 2e-5 * (1.0 + expected.abs()),
                    "torsion gradient mismatch at atom {i}, axis {axis}: actual={value}, expected={expected}"
                );
            }
        }
    }

    #[test]
    fn prepared_oop_gradient_matches_finite_difference() {
        let mol = chematic_smiles::parse("C=C(C)C").expect("isobutene");
        let model = Mmff94EnergyModel::new(&mol).expect("MMFF94 preparation");
        assert!(!model.oops.is_empty(), "trigonal center should prepare OOP");
        let coords = vec![
            [0.0, 0.0, 0.0],
            [1.34, 0.1, 0.0],
            [1.85, 1.35, 0.45],
            [1.95, -1.15, 0.25],
        ];
        let actual = model.oop_gradient(&coords);
        let delta = 1e-5;
        let mut work = coords.clone();
        for (i, row) in actual.iter().enumerate() {
            for (axis, &value) in row.iter().enumerate() {
                work[i][axis] += delta;
                let ep = model.oop_energy(&work);
                work[i][axis] -= 2.0 * delta;
                let em = model.oop_energy(&work);
                work[i][axis] += delta;
                let expected = (ep - em) / (2.0 * delta);
                assert!(
                    (value - expected).abs() < 3e-5 * (1.0 + expected.abs()),
                    "OOP gradient mismatch at atom {i}, axis {axis}: actual={value}, expected={expected}"
                );
            }
        }
        let term = model.oops[0];
        let mut singular = vec![[0.0; 3]; mol.atom_count()];
        singular[term.j] = [0.0, 0.0, 0.0];
        singular[term.i] = [1.0, 0.0, 0.0];
        singular[term.k] = [0.0, 1.0, 0.0];
        singular[term.l] = [0.0, 0.0, 1.0];
        let guarded = model.oop_gradient(&singular);
        assert!(
            guarded.iter().flatten().all(|value| value.abs() < 1e-12),
            "OOP gradient must fail closed at the asin branch point"
        );
    }

    #[test]
    fn bounded_analytic_gradient_matches_prepared_total_energy() {
        let mol = chematic_smiles::parse("CCCC").expect("butane");
        let model = Mmff94EnergyModel::new(&mol).expect("MMFF94 preparation");
        let coords = vec![
            [0.0, 0.0, 0.0],
            [1.54, 0.1, 0.0],
            [2.95, 0.85, 0.2],
            [4.25, 0.15, 1.05],
            [0.0, 1.05, 0.8],
            [1.55, -0.9, -0.75],
            [1.55, 0.25, 1.05],
            [2.95, 1.9, -0.55],
            [2.95, 1.05, 1.25],
            [4.2, -0.8, 1.45],
            [4.2, 0.0, 2.05],
            [4.9, 0.35, 0.35],
            [4.0, 0.95, 0.05],
            [4.6, 0.35, 1.55],
        ];
        let actual = model.bounded_analytic_gradient(&coords);
        let delta = 1e-5;
        let mut work = coords.clone();
        for (i, row) in actual.iter().enumerate() {
            for (axis, &value) in row.iter().enumerate() {
                work[i][axis] += delta;
                let ep = model.energy(&work);
                work[i][axis] -= 2.0 * delta;
                let em = model.energy(&work);
                work[i][axis] += delta;
                let expected = (ep - em) / (2.0 * delta);
                assert!(
                    (value - expected).abs() < 4e-5 * (1.0 + expected.abs()),
                    "combined prepared gradient mismatch at atom {i}, axis {axis}: actual={value}, expected={expected}"
                );
            }
        }
    }

    #[test]
    fn cutoff_bounded_analytic_gradient_matches_cutoff_total_energy() {
        let mol = chematic_smiles::parse("CCCC").expect("butane");
        let model = Mmff94EnergyModel::new(&mol).expect("MMFF94 preparation");
        let coords = vec![
            [0.0, 0.0, 0.0],
            [1.54, 0.1, 0.0],
            [2.95, 0.85, 0.2],
            [4.25, 0.15, 1.05],
            [0.0, 1.05, 0.8],
            [1.55, -0.9, -0.75],
            [1.55, 0.25, 1.05],
            [2.95, 1.9, -0.55],
            [2.95, 1.05, 1.25],
            [4.2, -0.8, 1.45],
            [4.2, 0.0, 2.05],
            [4.9, 0.35, 0.35],
            [4.0, 0.95, 0.05],
            [4.6, 0.35, 1.55],
        ];
        let actual = model.bounded_analytic_gradient_cutoff(&coords);
        let delta = 1e-5;
        let mut work = coords.clone();
        for (i, row) in actual.iter().enumerate() {
            for (axis, &value) in row.iter().enumerate() {
                work[i][axis] += delta;
                let ep = model.energy_cutoff_nonbonded(&work);
                work[i][axis] -= 2.0 * delta;
                let em = model.energy_cutoff_nonbonded(&work);
                work[i][axis] += delta;
                let expected = (ep - em) / (2.0 * delta);
                assert!(
                    (value - expected).abs() < 4e-5 * (1.0 + expected.abs()),
                    "cutoff combined gradient mismatch at atom {i}, axis {axis}: actual={value}, expected={expected}"
                );
            }
        }
    }
}

#[cfg(test)]
mod oop_boundary_contract_tests {
    use super::*;

    #[test]
    fn direct_and_prepared_out_of_plane_terms_agree_on_distorted_geometries() {
        for text in [
            "C=C(C)C",
            "CC(=O)N",
            "C[N+](=O)[O-]",
            "NC(=[NH2+])N",
            "c1ccccc1",
        ] {
            let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse(text).unwrap());
            let (types, view) = assign_mmff94_numeric_types_with_view(&mol).unwrap();
            let model = Mmff94EnergyModel::new(&mol).unwrap();
            assert!(!model.oops.is_empty(), "{text}");
            // A deterministic nonsymmetric geometry avoids accidentally planar centres.
            let coords: Vec<_> = (0..mol.atom_count())
                .map(|i| {
                    let t = i as f64 + 0.25;
                    [t.dcos(), t.dsin(), (t * 0.7).dsin()]
                })
                .collect();
            let direct = oop_energy(&view, &coords, &types);
            let prepared = model.oop_energy(&coords);
            // Fitted MMFF OOP constants can be negative (e.g. amide N).
            assert!(direct.is_finite() && direct.abs() > 0.0, "{text}: {direct}");
            assert!(
                (direct - prepared).abs() < 1e-10 * (1.0 + prepared.abs()),
                "{text}: {direct} != {prepared}"
            );
            let moved: Vec<_> = coords
                .iter()
                .map(|p| [-p[1] + 3.0, p[0] - 2.0, p[2] + 1.0])
                .collect();
            let transformed = oop_energy(&view, &moved, &types);
            assert!(
                (direct - transformed).abs() < 1e-10 * (1.0 + direct.abs()),
                "{text}"
            );
            assert_eq!(
                oop_energy(&view, &vec![[0.0; 3]; mol.atom_count()], &types),
                0.0
            );
        }
    }
}

#[cfg(test)]
mod optimizer_recovery_contract_tests {
    use super::*;

    #[test]
    fn recovery_step_backtracks_instead_of_overshooting_a_nearby_minimum() {
        let coords = [[1e-5, 0.0, 0.0]];
        let direction = [[-0.001, 0.0, 0.0]];
        let energy = |c: &[[f64; 3]]| 50.0 * c[0][0] * c[0][0];
        let before = energy(&coords);
        assert!(energy(&[[-0.00999, 0.0, 0.0]]) > before);
        let (after, value) =
            armijo_backtracking(&coords, &direction, (before, -1e-6), 10.0, &energy, &|_| {
                true
            })
            .expect("a smaller descending recovery step exists");
        assert!(value < before);
        assert!(after[0][0].abs() < coords[0][0].abs());
    }

    #[test]
    fn line_search_rejects_nonfinite_energy_and_respects_geometry_boundaries() {
        let coords = [[1.0, 0.0, 0.0]];
        let direction = [[-1.0, 0.0, 0.0]];
        assert!(
            armijo_backtracking(
                &coords,
                &direction,
                (1.0, -1.0),
                1.0,
                &|_| f64::NEG_INFINITY,
                &|_| true
            )
            .is_none()
        );
        assert!(
            armijo_backtracking(
                &coords,
                &direction,
                (1.0, -1.0),
                1.0,
                &|c| c[0][0] * c[0][0],
                &|_| false
            )
            .is_none()
        );
        assert_eq!(finite_max_gradient(&[[0.0, f64::NAN, 0.0]]), f64::INFINITY);
    }

    #[test]
    fn last_allowed_step_reports_the_final_gradient_convergence() {
        let mol = chematic_smiles::parse("CC").unwrap();
        let model = Mmff94EnergyModel::new(&mol).unwrap();
        let bond = model.bonds[0].params;
        let offset = 1.1e-4 / (MDYNE_A_TO_KCAL_MOL * bond.kb);
        let mut coords = vec![[0.0, 0.0, 0.0], [bond.r0 + offset, 0.0, 0.0]];
        assert!(finite_max_gradient(&model.bounded_analytic_gradient(&coords)) > 1e-4);
        let result = model
            .minimize_lbfgs_bounded_analytic(&mut coords, 1)
            .unwrap();
        assert!(finite_max_gradient(&model.bounded_analytic_gradient(&coords)) < 1e-4);
        assert!(result.converged);
        assert_eq!(result.iterations, 1);
        assert_eq!(
            result.termination,
            Mmff94TerminationReason::GradientConverged
        );
    }

    #[test]
    fn zero_iteration_budget_reports_stationarity_of_the_returned_geometry() {
        let mol = chematic_smiles::parse("CC").unwrap();
        let model = Mmff94EnergyModel::new(&mol).unwrap();
        let r0 = model.bonds[0].params.r0;
        let mut coords = vec![[0.0, 0.0, 0.0], [r0, 0.0, 0.0]];
        let stationary = model
            .minimize_lbfgs_bounded_analytic(&mut coords, 0)
            .unwrap();
        assert!(stationary.converged);
        assert_eq!(stationary.iterations, 0);
        assert_eq!(
            stationary.termination,
            Mmff94TerminationReason::GradientConverged
        );
        coords[1][0] += 0.1;
        let nonstationary = model
            .minimize_lbfgs_bounded_analytic(&mut coords, 0)
            .unwrap();
        assert!(!nonstationary.converged);
        assert_eq!(nonstationary.iterations, 0);
    }
}
