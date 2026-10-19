//! Eager E/Z reference kept independent of the lazy ranker.
use super::*;
const CORPORA: [&str; 2] = [
    include_str!("../../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"),
    include_str!("../../../validation/cosmolkit_comparison/corpora/stress-v1.smi"),
];
fn eager_assign_ez(
    mol: &Molecule,
    bond_idx: BondIdx,
    ranker: &chematic_cip::SubstituentRanker,
) -> Option<(AtomIdx, CipCode)> {
    let bond = mol.bond(bond_idx);
    // A ring of fewer than eight atoms holds its double bond cis: no E/Z
    // (RDKit drops such stereo too).
    if bond.order != BondOrder::Double
        || !end_has_marked_substituent(mol, bond.atom1, bond.atom2)
        || !end_has_marked_substituent(mol, bond.atom2, bond.atom1)
        || double_bond_in_small_ring(mol, bond_idx)
    {
        return None;
    }
    let a1 = bond.atom1;
    let a2 = bond.atom2;
    let subs = |end: AtomIdx, other: AtomIdx| -> Vec<AtomIdx> {
        mol.neighbors(end)
            .filter(|&(nb, bidx)| nb != other && mol.bond(bidx).order != BondOrder::Double)
            .map(|(nb, _)| nb)
            .collect()
    };
    let subs_a1 = subs(a1, a2);
    let subs_a2 = subs(a2, a1);
    if subs_a1.is_empty() || subs_a2.is_empty() {
        return None;
    }
    let (_, up_a1) = eager_highest_sub(mol, a1, &subs_a1, ranker)?;
    let (_, up_a2) = eager_highest_sub(mol, a2, &subs_a2, ranker)?;
    let code = if up_a1 == up_a2 {
        CipCode::Z
    } else {
        CipCode::E
    };
    Some((a1, code))
}

/// [`highest_stereo_sub`] with the accurate comparator (at most two
/// substituents on a double-bond end).
fn eager_highest_sub(
    mol: &Molecule,
    alkene_end: AtomIdx,
    subs: &[AtomIdx],
    ranker: &chematic_cip::SubstituentRanker,
) -> Option<(AtomIdx, bool)> {
    // No marker on either substituent: nothing to read (see
    // `highest_stereo_sub`).
    if subs
        .iter()
        .all(|&sub| substituent_is_up(mol, alkene_end, sub).is_none())
    {
        return None;
    }
    let (top, other) = match subs {
        [only] => (*only, None),
        [a, b] => match ranker.compare(mol, alkene_end, *a, *b).ok()?? {
            std::cmp::Ordering::Greater => (*a, Some(*b)),
            std::cmp::Ordering::Less => (*b, Some(*a)),
            std::cmp::Ordering::Equal => return None,
        },
        _ => return None,
    };
    if let Some(up) = substituent_is_up(mol, alkene_end, top) {
        return Some((top, up));
    }
    let other_up = substituent_is_up(mol, alkene_end, other?)?;
    Some((top, !other_up))
}

#[test]
fn lazy_ez_matches_eager_on_raw_and_kekule_corpora() {
    let mut count = 0;
    for text in CORPORA {
        for line in text.lines().filter(|l| !l.trim().is_empty()) {
            let smiles = line.split_whitespace().next().unwrap();
            let mol = chematic_smiles::parse(smiles).unwrap();
            let mut forms = vec![mol.clone()];
            if let Ok(kekule) = chematic_core::kekulization::kekulize(&mol) {
                forms.push(chematic_core::kekulization::apply_kekule(&mol, &kekule));
            }
            for m in forms {
                let eager = chematic_cip::SubstituentRanker::new(
                    &m,
                    chematic_cip::CipBudget::default_budget(),
                );
                let expected: Vec<_> = (0..m.bond_count())
                    .filter_map(|i| {
                        let b = BondIdx(i as u32);
                        eager_assign_ez(&m, b, &eager).map(|(_, c)| (b, c))
                    })
                    .collect();
                assert_eq!(
                    assign_ez_bonds_with_mode(&m, CipMode::Accurate),
                    expected,
                    "{smiles}"
                );
            }
            count += 1;
        }
    }
    assert_eq!(count, 14_072);
}
#[test]
fn ranker_initializes_only_when_two_explicit_substituents_need_comparison() {
    for (smiles, initialized) in [
        ("c1ccccc1", false),
        ("F/C=C/F", false),
        ("F/C(Cl)=C(/Br)I", true),
    ] {
        let m = chematic_smiles::parse(smiles).unwrap();
        let ranker = LazyEzRanker::new(&m, chematic_cip::CipBudget::default_budget());
        for i in 0..m.bond_count() {
            assign_ez_accurate(&m, BondIdx(i as u32), &ranker);
        }
        assert_eq!(ranker.value.get().is_some(), initialized, "{smiles}");
    }
}
