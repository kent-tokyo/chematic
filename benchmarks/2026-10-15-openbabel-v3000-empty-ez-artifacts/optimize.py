from pathlib import Path
s=Path('/private/tmp/chematic-openbabel-ez-empty-source-20261015');r=Path('/private/tmp/chematic-openbabel-ez-empty-20261015');p=s/'crates/chematic-perception/src/stereo2d_ez_direction.rs';old=p.read_text();a=old.index('pub fn apply_ez_directions_from_2d_ex(');b=old.index('\n// ---------------------------------------------------------------------------\n// Classification',a);reference=old[a:b].replace('pub fn apply_ez_directions_from_2d_ex(', 'fn frozen_reference(',1)
start=old.index('    // One-time, non-mutating aromaticity query',a);end=old.index('    // Phase 2:',start)
replacement='''    // With no literal double bond, every later phase has no outcomes and
    // cannot write a direction or diagnostic. Avoid the unused ring work.
    if double_bonds.is_empty() {
        return Vec::new();
    }
    // Exclude Kekule aromatic ring bonds with the same non-mutating query.
    let aromaticity = assign_aromaticity(mol);
    // Classify independently; no double bond reads another bond's outcome.
    let mut outcomes = classify_double_bonds(
        mol,
        coords,
        &double_bonds,
        explicitly_unspecified,
        &aromaticity,
    );

'''
t=old[:start]+replacement+old[end:]
helper='''#[inline]
fn classify_double_bonds(
    mol: &Molecule,
    coords: &[(f64, f64)],
    double_bonds: &[BondIdx],
    explicitly_unspecified: &HashSet<BondIdx>,
    aromaticity: &AromaticityModel,
) -> HashMap<BondIdx, EzOutcome> {
    let mut outcomes = HashMap::with_capacity(double_bonds.len());
    for &bidx in double_bonds {
        let outcome = classify_double_bond(mol, coords, bidx, explicitly_unspecified, aromaticity);
        outcomes.insert(bidx, outcome);
    }
    outcomes
}

'''
loc=t.index('/// Per-double-bond classification');t=t[:loc]+helper+t[loc:]
t+='\n#[cfg(test)]\n#[path = "stereo2d_ez_empty_tests.rs"]\nmod empty_double_tests;\n';p.write_text(t)
tests='''//! Compare the empty-double guard and extracted classification with the frozen API.
use super::*;

'''+reference+'''
fn compare(input: &Molecule, coords: &[(f64, f64)], unspecified: &HashSet<BondIdx>) {
    let mut expected = input.clone();
    let mut actual = input.clone();
    let old = frozen_reference(&mut expected, coords, unspecified);
    let new = apply_ez_directions_from_2d_ex(&mut actual, coords, unspecified);
    assert_eq!(old, new);
    assert_eq!(expected.atoms().collect::<Vec<_>>(), actual.atoms().collect::<Vec<_>>());
    assert_eq!(expected.bonds().collect::<Vec<_>>(), actual.bonds().collect::<Vec<_>>());
    assert_eq!(expected.stereo_groups(), actual.stereo_groups());
    assert_eq!(expected.rdkit_bond_order(), actual.rdkit_bond_order());
    for (i, _) in expected.atoms() {
        assert_eq!(expected.neighbor_slice(i), actual.neighbor_slice(i));
        assert_eq!(expected.stereo_neighbor_order(i), actual.stereo_neighbor_order(i));
        assert_eq!(expected.atom_tag(i), actual.atom_tag(i));
        assert_eq!(expected.r_group_label(i), actual.r_group_label(i));
    }
    for (i, _) in expected.bonds() {
        assert_eq!(expected.bond_direction(i), actual.bond_direction(i));
        assert_eq!(expected.bond_direction_anchor(i), actual.bond_direction_anchor(i));
        assert_eq!(expected.smiles_ring_closure_metadata(i), actual.smiles_ring_closure_metadata(i));
    }
    for i in 0..expected.stereo_groups().len() {
        assert_eq!(expected.stereo_group_bond_indices(i), actual.stereo_group_bond_indices(i));
    }
}

#[test]
fn no_double_bonds_preserve_existing_annotations_and_ignore_unused_coordinates() {
    let empty = chematic_core::MoleculeBuilder::new().build();
    compare(&empty, &[], &HashSet::from([BondIdx(u32::MAX)]));
    for smiles in ["CCCC", "C1CCCCC1", "C12C3C4C1C5C2C3C45", "c1ccccc1", "[Na+].[Cl-]", "F[C@](Cl)(Br)I", "[*:1]CC", "CCN->[Cu]"] {
        let mut mol = chematic_smiles::parse_template(smiles).unwrap();
        assert!(!mol.bonds().any(|(_, b)| b.order == BondOrder::Double));
        mol.set_tag(AtomIdx(0), Some(91));
        if mol.bond_count() > 0 { mol.set_bond_direction(BondIdx(0), BondOrder::Up); }
        for coords in [Vec::new(), vec![(0.0, 0.0); mol.atom_count()], vec![(f64::NAN, f64::INFINITY); mol.atom_count()]] {
            compare(&mol, &coords, &HashSet::from([BondIdx(u32::MAX)]));
        }
    }
}

#[test]
fn fixed_ten_thousand_molecules_keep_all_ez_results() {
    let corpus = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"));
    let mut rows = 0;
    for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
        let mol = chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
        let coords = (0..mol.atom_count()).map(|i| (i as f64, (i * i % 7) as f64)).collect::<Vec<_>>();
        compare(&mol, &coords, &HashSet::new());
        rows += 1;
    }
    assert_eq!(rows, 10000);
}

#[test]
fn double_query_dative_and_unspecified_cases_match_reference() {
    for smiles in ["CC=CC", "C/C=C/C=C/C", "CC=C=C(C)C", "C1=CC=CC=C1", "CC1=CCCCC1", "N->[Cu]", "CC"] {
        let mol = chematic_smiles::parse_template(smiles).unwrap();
        let unspecified = mol.bonds().filter_map(|(i,b)| (b.order == BondOrder::Double).then_some(i)).collect::<HashSet<_>>();
        for coords in [Vec::new(), vec![(0.0,0.0);mol.atom_count()], (0..mol.atom_count()).map(|i| (i as f64,(i%2) as f64)).collect()] {
            compare(&mol,&coords,&HashSet::new());
            compare(&mol,&coords,&unspecified);
        }
    }
    for order in [BondOrder::QueryAny, BondOrder::QueryDoubleOrAromatic, BondOrder::Aromatic] {
        let mut mol = chematic_smiles::parse_template("CC").unwrap();
        mol.set_bond_order(BondIdx(0),order);
        compare(&mol,&[],&HashSet::new());
    }
}

#[test]
fn graph_edits_recheck_the_double_bond_gate() {
    let mut mol = chematic_smiles::parse_template("CC=CC").unwrap();
    let coords = [(0.0,1.0),(1.0,0.0),(2.0,0.0),(3.0,1.0)];
    for order in [BondOrder::Single, BondOrder::Double, BondOrder::Single, BondOrder::Double] {
        mol.set_bond_order(BondIdx(1),order);
        compare(&mol,&coords,&HashSet::new());
    }
}
'''
(s/'crates/chematic-perception/src/stereo2d_ez_empty_tests.rs').write_text(tests)
p=s/'CHANGELOG.md';t=p.read_text().replace('### Performance\n','### Performance\n\n- Skip unused E/Z aromaticity work for molecules without literal double\n  bonds, retaining existing stereo annotations and all nonempty-bond\n  classification and diagnostic behavior.\n',1);p.write_text(t)
print('Prepared no-double E/Z guard and frozen full-API differential tests')
