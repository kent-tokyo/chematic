from pathlib import Path
r=Path('/private/tmp/chematic-openbabel-ez-gate-20261015');s=Path('/private/tmp/chematic-openbabel-ez-gate-source-20261015');p=s/'crates/chematic-perception/src/stereo2d_ez_direction.rs';old=p.read_text();start=old.index('pub fn apply_ez_directions_from_2d_ex(');end=old.index('// ---------------------------------------------------------------------------\n// Classification',start)
reference=old[start:end].replace('pub fn apply_ez_directions_from_2d_ex(', 'fn reference_ez(',1).replace('classify_double_bonds(', 'reference_classify_double_bonds(')
a=old.index('#[inline]\nfn classify_double_bonds(');b=old.index('/// Per-double-bond classification',a);ref_classifier=old[a:b].replace('fn classify_double_bonds(', 'fn reference_classify_double_bonds(',1)
x=old.index('    // With no literal double bond',start);y=old.index('    // Phase 2:',x)
t=old[:x]+'''    // Classification returns no outcomes when every double bond is already
    // known to be NotRequested from the existing topology rules.
    let Some(mut outcomes) =
        classify_double_bonds(mol, coords, &double_bonds, explicitly_unspecified)
    else {
        return Vec::new();
    };

'''+old[y:]
a=t.index('#[inline]\nfn classify_double_bonds(');b=t.index('/// Per-double-bond classification',a)
classifier='''#[inline]
fn classify_double_bonds(
    mol: &Molecule,
    coords: &[(f64, f64)],
    double_bonds: &[BondIdx],
    explicitly_unspecified: &HashSet<BondIdx>,
) -> Option<HashMap<BondIdx, EzOutcome>> {
    // Both tests imply NotRequested regardless of aromaticity or coordinates.
    // Terminal ends qualify only after preserving the earlier cumulene and
    // explicit-unspecified diagnostics. Small rings are unconditionally out.
    if double_bonds.iter().all(|&bidx| {
        let bond = mol.bond(bidx);
        let terminal = !explicitly_unspecified.contains(&bidx)
            && !has_other_double_bond(mol, bond.atom1, bidx)
            && !has_other_double_bond(mol, bond.atom2, bidx)
            && [bond.atom1, bond.atom2].into_iter().any(|end| {
                let other = if end == bond.atom1 { bond.atom2 } else { bond.atom1 };
                mol.neighbors(end)
                    .all(|(nb, bi)| nb == other || mol.bond(bi).order == BondOrder::Double)
            });
        terminal || in_ring_smaller_than(mol, bidx, 8)
    }) {
        return None;
    }
    // Any potentially requested bond retains the original complete query,
    // classification order, diagnostics and carrier-conflict resolution.
    let aromaticity = assign_aromaticity(mol);
    let mut outcomes = HashMap::with_capacity(double_bonds.len());
    for &bidx in double_bonds {
        let outcome = classify_double_bond(mol, coords, bidx, explicitly_unspecified, &aromaticity);
        outcomes.insert(bidx, outcome);
    }
    Some(outcomes)
}

'''
t=t[:a]+classifier+t[b:];t+='\n#[cfg(test)]\n#[path = "stereo2d_ez_inactive_tests.rs"]\nmod inactive_double_tests;\n';p.write_text(t)
previous=(s/'crates/chematic-perception/src/stereo2d_ez_empty_tests.rs').read_text();x=previous.index('fn compare(');y=previous.index('#[test]',x);compare=previous[x:y].replace('frozen_reference(', 'reference_ez(')
# Compare the same independently cloned molecules again after the first call.
compare=compare.replace('    let old = reference_ez', '    for _ in 0..2 {\n    let old = reference_ez',1);closing=compare.rfind('}');compare=compare[:closing]+'    }\n'+compare[closing:]
test='//! Existing pinned inputs and exact pre-change API/classifier oracle.\nuse super::*;\n\n'+reference+ref_classifier+compare+'''
#[test]
fn inactive_double_gate_preserves_full_api_on_pinned_inputs() {
    let corpora = [include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi")),include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/cosmolkit_comparison/corpora/stress-v1.smi"))];
    let mut rows=0;
    let mut comparisons=0;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            let parsed=chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
            let kekule=chematic_core::kekulize(&parsed).map(|k| chematic_core::apply_kekule(&parsed,&k)).unwrap_or_else(|_|parsed.clone());
            for mut mol in [parsed,kekule] {
                // Preserve annotations too; gate eligibility does not read them.
                if rows % 2 == 1 {
                    let ids=mol.bonds().filter_map(|(bi,b)| (bi.0 % 3 == 0 && b.order != BondOrder::Double).then_some(bi)).collect::<Vec<_>>();
                    for bi in ids { mol.set_bond_direction(bi,BondOrder::Down); }
                    if mol.atom_count()>0 { mol.set_tag(AtomIdx(0),Some(91)); }
                }
                let all=mol.bonds().filter_map(|(bi,b)| (b.order==BondOrder::Double).then_some(bi)).collect::<HashSet<_>>();
                let even=all.iter().copied().filter(|bi|bi.0 % 2 == 0).collect::<HashSet<_>>();
                let normal=(0..mol.atom_count()).map(|i|(i as f64,(i*i % 7) as f64)).collect::<Vec<_>>();
                for coords in [normal,Vec::new(),vec![(f64::NAN,f64::INFINITY);mol.atom_count()]] {
                    for unspecified in [&HashSet::new(),&all,&even] {
                        compare(&mol,&coords,unspecified);
                        comparisons+=2;
                    }
                }
            }
            rows+=1;
        }
    }
    assert_eq!(rows,14072);
    assert_eq!(comparisons,506592);
}
''';(s/'crates/chematic-perception/src/stereo2d_ez_inactive_tests.rs').write_text(test)
p=s/'CHANGELOG.md';p.write_text(p.read_text().replace('### Performance\n','### Performance\n\n- Skip E/Z aromaticity and carrier work when all literal double bonds are\n  already excluded by the existing small-ring or terminal-end rules, while\n  preserving cumulene and explicitly unspecified stereo diagnostics.\n',1));print('Inactive E/Z topology candidate and 506592-case frozen oracle written')
