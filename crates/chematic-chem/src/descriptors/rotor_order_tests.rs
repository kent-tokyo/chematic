//! Freeze the previous strict SMARTS and pair filtering independently of the
//! reordered predicates. Compare endpoints, not only their count.
use super::*;

const FROZEN_STRICT_SMARTS: &str = "[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)\
&!$(C([CH3])([CH3])[CH3])&!$([CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])\
&!$([#7,O,S!D1]-!@[CD3]=[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]\
-,:;!@[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])\
&!$([CH3])]";

fn frozen_pairs(mol: &Molecule) -> Vec<(AtomIdx, AtomIdx)> {
    static QUERY: OnceLock<chematic_smarts::QueryMolecule> = OnceLock::new();
    let query = QUERY.get_or_init(|| parse_smarts(FROZEN_STRICT_SMARTS).unwrap());
    let mut pairs: Vec<_> = chematic_smarts::find_match_atom_sets_perceived(
        query,
        mol,
        &chematic_smarts::MatchConfig::default(),
    )
    .into_iter()
    .filter_map(|set| match set.as_slice() {
        &[a, b] => Some((AtomIdx(a.min(b) as u32), AtomIdx(a.max(b) as u32))),
        _ => None,
    })
    .filter(|&(a, b)| {
        [a, b].into_iter().all(|end| {
            mol.neighbors(end)
                .filter(|&(nb, _)| !crate::hydrogen::is_removable_explicit_h(mol, nb))
                .count()
                > 1
        })
    })
    .collect();
    pairs.sort_unstable_by_key(|&(a, b)| (a.0, b.0));
    pairs.dedup();
    pairs
}

fn compare(mol: &Molecule) {
    let expected = frozen_pairs(mol);
    assert_eq!(rotatable_bond_atom_pairs(mol), expected);
    assert_eq!(rotatable_bond_count(mol), expected.len());
}

#[test]
fn reordered_rotor_predicates_keep_every_pinned_pair() {
    let corpora = [
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"
        )),
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/cosmolkit_comparison/corpora/stress-v1.smi"
        )),
    ];
    let mut rows = 0;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            let mol =
                chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
            compare(&mol);
            let kekule = chematic_core::kekulize(&mol)
                .map(|k| chematic_core::apply_kekule(&mol, &k))
                .unwrap_or_else(|_| mol.clone());
            compare(&kekule);
            rows += 1;
        }
    }
    assert_eq!(rows, 14072);
}

#[test]
fn rotor_exclusions_keep_explicit_hydrogens_and_functional_groups() {
    compare(&chematic_core::MoleculeBuilder::new().build());
    for text in [
        "[H]NCC",
        "[2H]NCC",
        "[2H]C([2H])([2H])NC=O",
        "N#CCCCC",
        "CC(F)(F)F",
        "CCC(Cl)(Cl)Cl",
        "CCC(Br)(Br)Br",
        "CCC(C)(C)C",
        "CC(=O)NCC",
        "CCC(=N)OCC",
        "CCC(=[NH2+])NCC",
        "CCON=C(C)C",
        "CSSC",
        "c1ccccc1-c1ccccc1",
    ] {
        compare(&chematic_smiles::parse_template(text).unwrap());
    }
}
