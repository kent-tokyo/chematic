//! Tests extracted from `transform.rs`; kept as a child module so private helpers remain testable.

use super::*;
use chematic_smiles::parse;

fn canonical(mol: &Molecule) -> String {
    chematic_smiles::canonical_smiles(mol)
}

fn canonical_set(set: Vec<Molecule>) -> Vec<String> {
    set.into_iter().map(|mol| canonical(&mol)).collect()
}

#[test]
fn centre_gaining_a_neighbour_in_place_of_its_h_keeps_stereo() {
    // Native: the new bond takes the implicit H's place (geometric
    // retention). RDKit 2026.03.6 copies the raw tag onto a bond list
    // with the new bond appended, which is the opposite arrangement; the
    // RDKit profile reproduces that (its `RunReactants`, sanitized).
    for (smirks, smiles, native, rdkit) in [
        (
            "[C;$(CO):1]>>[C:1]C(=O)C",
            "C[C@H](N)O",
            "C[C@](C(C)=O)(N)O",
            "CC(=O)[C@](C)(N)O",
        ),
        (
            "[CH1:1]>>[C:1]F",
            "C[C@@H](Cl)CC",
            "C[C@@](F)(Cl)CC",
            "CC[C@@](C)(F)Cl",
        ),
        (
            "[C:1]>>[C:1]C",
            "F[C@H](Cl)Br",
            "F[C@](C)(Cl)Br",
            "C[C@](F)(Cl)Br",
        ),
    ] {
        let mol = parse(smiles).unwrap();
        let got: Vec<String> = run_reactants(smirks, &[&mol])
            .unwrap()
            .into_iter()
            .flatten()
            .map(|m| canonical(&m))
            .collect();
        let want = canonical(&parse(native).unwrap());
        assert!(
            !got.is_empty() && got.iter().all(|g| *g == want),
            "native {smirks} on {smiles}: {got:?}"
        );
        let RdkitProfileOutcome::Report(report) = run_reactants_traced_rdkit_2026_03_6(
            smirks,
            &[&mol],
            &ReactionTransformLimits::default(),
        )
        .unwrap() else {
            panic!("RDKit profile refused {smirks}");
        };
        let got: Vec<String> = report
            .products
            .iter()
            .flatten()
            .map(|p| canonical(&p.molecule))
            .collect();
        let want = canonical(&parse(rdkit).unwrap());
        assert!(
            !got.is_empty() && got.iter().all(|g| *g == want),
            "RDKit profile {smirks} on {smiles}: {got:?}"
        );
    }
}

#[test]
fn identity_single_atom() {
    let mol = parse("C").unwrap();
    let results = run_reactants("[C:1]>>[C:1]", &[&mol]).unwrap();
    assert_eq!(results.len(), 1);
    assert_eq!(results[0].len(), 1);
    assert_eq!(results[0][0].atom_count(), 1);
}

#[test]
fn carried_bonds_follow_product_atom_order() {
    // Carried atoms are added breadth-first from the template atoms in
    // reactant order and their bonds in product atom order, not in hash
    // map order (BTMR1032's Kekulé product changed with a map's
    // capacity, #734).
    let mol = parse("OC(=O)c1ccc2ccccc2c1CCN(C)C").unwrap();
    let results = run_reactants("[C:1](=[O:2])[OH:3]>>[C:1](=[O:2])Cl", &[&mol]).unwrap();
    let product = &results[0][0];
    let carried: Vec<u32> = product
        .bonds()
        .skip(3)
        .map(|(_, b)| b.atom1.0.min(b.atom2.0))
        .collect();
    assert!(carried.windows(2).all(|w| w[0] <= w[1]), "{carried:?}");
    // The same product however often the reaction runs.
    for _ in 0..3 {
        let again = run_reactants("[C:1](=[O:2])[OH:3]>>[C:1](=[O:2])Cl", &[&mol]).unwrap();
        let bonds = |m: &Molecule| {
            m.bonds()
                .map(|(_, b)| (b.atom1.0, b.atom2.0, b.order))
                .collect::<Vec<_>>()
        };
        assert_eq!(bonds(&again[0][0]), bonds(product));
    }
}

#[test]
fn mapped_asymmetric_ether_keeps_both_embeddings_on_one_atom_set() {
    // The C-O-C query covers the same three target atoms in both
    // orientations, but the mapped cleavage yields two different product
    // tuples. Target-set uniquification must not discard either one.
    let mol = parse("COCC").unwrap();
    let results = run_reactants("[C:1][O:2][C:3]>>[C:1][O:2].[C:3]", &[&mol]).unwrap();
    let outcomes: std::collections::BTreeSet<Vec<String>> = results
        .into_iter()
        .map(|set| {
            let mut products = canonical_set(set);
            products.sort();
            products
        })
        .collect();
    let expected: std::collections::BTreeSet<Vec<String>> = [vec!["C", "CCO"], vec!["CC", "CO"]]
        .into_iter()
        .map(|set| {
            let mut products: Vec<String> = set
                .into_iter()
                .map(|smiles| canonical(&parse(smiles).unwrap()))
                .collect();
            products.sort();
            products
        })
        .collect();
    assert_eq!(outcomes, expected);
}

#[test]
fn applies_atomic_number_smirks_and_keeps_atom_maps() {
    let reactant = parse("NC=O").unwrap();
    let smirks = "[#7:1][C:2](=[O:3])>>[#7:1][C:2](=[O:3])";
    let normalized = crate::reaction::expand_atomic_number_primitives(smirks)
        .unwrap()
        .into_iter()
        .next()
        .unwrap();
    let matches = find_reaction_matches(&normalized, &[&reactant]).unwrap();
    assert_eq!(matches.len(), 1);
    assert_eq!(
        matches[0].atom_map_positions(&normalized).unwrap().len(),
        3,
        "map identity must survive query matching"
    );
    let results = run_reactants(smirks, &[&reactant]).unwrap();
    assert!(!results.is_empty(), "atomic-number SMIRKS must apply");
    let product = &results[0][0];
    assert_eq!(product.atom_count(), 3);
}

#[test]
fn styrene_hydrogenation_keeps_both_template_map_orientations() {
    let reactant = parse("C=Cc1ccccc1").unwrap();
    let smirks = "[C:1]=[C:2]>>[C:1][C:2]";
    let prepared = PreparedReaction::new(smirks).unwrap();
    let matches = prepared.find_matches(&[&reactant]).unwrap();
    let maps: Vec<_> = matches
        .iter()
        .map(|m| m.atom_map_positions(smirks).unwrap())
        .collect();
    assert_eq!(
        maps.len(),
        2,
        "each atom-map orientation is observable: {maps:?}"
    );
    let orientations: std::collections::BTreeSet<_> = maps
        .iter()
        .map(|positions| (positions[&1].1.0, positions[&2].1.0))
        .collect();
    assert_eq!(orientations, [(0, 1), (1, 0)].into());
}

#[test]
fn reactant_side_smarts_queries_apply(/* issue #734 */) {
    // The SMIRKS reactant side now parses as full SMARTS: query-only
    // features (D, X, ;, ',', !, $()) match and the identity transform
    // returns the input, where before every one of these failed to parse.
    let canon = |m: &Molecule| chematic_smiles::canonical_smiles(m);
    let ethane = parse("CC").unwrap();
    let ethanol = parse("CCO").unwrap();
    for (smirks, reactant) in [
        ("[CD1:1]>>[C:1]", &ethane),
        ("[CX4:1]>>[C:1]", &ethane),
        ("[C;H3:1]>>[C:1]", &ethane),
        ("[C,N:1]>>[C:1]", &ethane),
        ("[!O:1]>>[C:1]", &ethane),
        ("[$([OH]):1]>>[O:1]", &ethanol),
    ] {
        let products = run_reactants(smirks, &[reactant]).unwrap();
        assert!(!products.is_empty(), "{smirks}: expected a product");
        assert!(
            products.iter().all(|set| canon(&set[0]) == canon(reactant)),
            "{smirks}: identity transform must return the input"
        );
    }
    // Product-side query features are specifications in RDKit's sense:
    // query-only parts are dropped (`[OX2H1]` applies H1), a mapped atom
    // without a single element keeps the reactant's element, and an
    // unmapped one cannot be built.
    let anisole = parse("COc1ccccc1").unwrap();
    let out = run_reactants("[c:1][OX2:2][CH3:3]>>[c:1][OX2H1:2]", &[&anisole]).unwrap();
    assert_eq!(out.len(), 1);
    assert_eq!(canon(&out[0][0]), canon(&parse("Oc1ccccc1").unwrap()));
    let out = run_reactants("[C:1]>>[C,N:1]", &[&ethane]).unwrap();
    assert!(out.iter().all(|set| canon(&set[0]) == canon(&ethane)));
    assert!(matches!(
        run_reactants("[C:1]>>[C:1][C,N]", &[&ethane]),
        Err(TransformError::SmirksParse(
            crate::reaction::RxnError::UnsupportedProductPrimitive { .. }
        ))
    ));
    // `;H<n>` is the SMARTS spelling of an H count.
    let out = run_reactants("[C:1]>>[C;H3:1]", &[&ethane]).unwrap();
    assert!(!out.is_empty());
}

#[test]
fn mapped_atom_takes_a_changed_product_element(/* issue #734 item 2 */) {
    // RDKit 2026.03.6 reference products for each case.
    let canon = |m: &Molecule| chematic_smiles::canonical_smiles(m);
    for (smirks, reactant, expected) in [
        ("[C:1]>>[N:1]", "CC", vec!["CN"]),
        ("[#6:1]>>[#7:1]", "CC", vec!["CN"]),
        ("[CH3:1]>>[N:1]", "CC", vec!["CN"]),
        ("[C:1][O:2]>>[C:1][S:2]", "CO", vec!["CS"]),
        ("[c:1]>>[n:1]", "c1ccccc1", vec!["c1ccncc1"]),
        ("[*:1]>>[N:1]", "C", vec!["N"]),
        // No change: wildcard product, or the product element equals the
        // reactant query's (first) element.
        ("[C:1]>>[*:1]", "CO", vec!["CO"]),
        ("[C,N:1]>>[C:1]", "CN", vec!["CC", "CN"]),
        ("[C,N:1]>>[N:1]", "CC", vec!["CN"]),
    ] {
        let mol = parse(reactant).unwrap();
        let mut got: Vec<String> = run_reactants(smirks, &[&mol])
            .unwrap()
            .iter()
            .map(|set| canon(&set[0]))
            .collect();
        got.sort();
        got.dedup();
        let mut want: Vec<String> = expected.iter().map(|s| canon(&parse(s).unwrap())).collect();
        want.sort();
        assert_eq!(got, want, "{smirks} on {reactant}");
    }
}

#[test]
fn charge_change_keeps_rdkit_hydrogen_count(/* issue #734 item 4 */) {
    // RDKit 2026.03.6 reference products.
    let canon = |m: &Molecule| chematic_smiles::canonical_smiles(m);
    for (smirks, reactant, expected) in [
        ("[N:1]>>[N+2:1]", "CN", "C[NH2+2]"),
        ("[N:1]>>[N++:1]", "CN", "C[NH2+2]"),
        ("[N:1]>>[N+:1]", "CN", "C[NH3+]"),
        ("[N:1]>>[N-:1]", "CN", "C[NH-]"),
        ("[O:1]>>[O-:1]", "CO", "C[O-]"),
        ("[C:1]>>[C-:1]", "CC", "[CH2-]C"),
        ("[C:1]>>[C+:1]", "CC", "[CH2+]C"),
    ] {
        let mol = parse(reactant).unwrap();
        let out = run_reactants(smirks, &[&mol]).unwrap();
        assert!(!out.is_empty(), "{smirks}");
        assert!(
            out.iter()
                .all(|set| canon(&set[0]) == canon(&parse(expected).unwrap())),
            "{smirks} on {reactant}: {:?}",
            out.iter().map(|s| canon(&s[0])).collect::<Vec<_>>()
        );
    }
}

#[test]
fn explicit_hydrogen_reactants_give_the_implicit_products(/* issue #734 item 5 */) {
    // Every case the issue lists: the explicit-H reactant must give the
    // same products (compared without explicit H) as its implicit form.
    let canon = |m: &Molecule| chematic_smiles::canonical_smiles(m);
    let strip = |m: &Molecule| canon(&chematic_chem::remove_hydrogens(m));
    for (smirks, reactant) in [
        ("[C:1]-[C:2]>>[C:1]=[C:2]", "CC"),
        ("[O:1]>>[O:1]C(=O)C", "CCO"),
        ("[C:1]>>[C:1]Cl", "CC"),
        ("[O:1]>>[O:1]C", "CCO"),
        ("[C:1][O:2]>>[C:1]=[O:2]", "CO"),
        ("[C:1][O:2]>>[C:1]=[O:2]", "CCO"),
        ("[C:1][N:2]>>[C:1]=[N+:2]", "CN(C)C"),
        ("[CH3:1]>>[CH3:1]", "CC"),
        ("[C:1]>>[C:1]", "CC"),
        ("[N:1]>>[N+:1]", "CN"),
        ("[C:1][O:2]>>[C:1].[O:2]", "CCO"),
        (
            "[#6:1](-[OX2:2])-;@[CX4&+0:3]>>[*:1]-[*:3]",
            "CC(=O)C1=C(O)CN(C2OCCc3ccccc32)C1=O",
        ),
    ] {
        let implicit = parse(reactant).unwrap();
        let explicit = chematic_chem::add_hydrogens(&implicit);
        let mut want: Vec<String> = run_reactants(smirks, &[&implicit])
            .unwrap()
            .iter()
            .map(|set| set.iter().map(canon).collect::<Vec<_>>().join("."))
            .collect();
        let mut got: Vec<String> = run_reactants(smirks, &[&explicit])
            .unwrap()
            .iter()
            .map(|set| set.iter().map(&strip).collect::<Vec<_>>().join("."))
            .collect();
        want.sort();
        want.dedup();
        got.sort();
        got.dedup();
        assert!(!want.is_empty(), "{smirks}: implicit form must react");
        assert_eq!(got, want, "{smirks} on explicit-H {reactant}");
    }
    // The edited atoms keep the input's explicit-H mode: ethane with
    // explicit H gives ethylene as C2 plus four H atoms.
    let ethane = chematic_chem::add_hydrogens(&parse("CC").unwrap());
    let out = run_reactants("[C:1]-[C:2]>>[C:1]=[C:2]", &[&ethane]).unwrap();
    assert!(!out.is_empty());
    for set in out {
        assert_eq!(set[0].atom_count(), 6);
        assert!(
            set[0]
                .atoms()
                .all(|(i, _)| chematic_core::implicit_hcount(&set[0], i) == 0)
        );
    }
    // A template that matches explicit hydrogens still sees them.
    let methane = chematic_chem::add_hydrogens(&parse("C").unwrap());
    assert!(
        !run_reactants("[C:1][H]>>[C:1]O", &[&methane])
            .unwrap()
            .is_empty()
    );
}

#[test]
fn reactants_follow_rdkit_sanitized_input(/* issues #734, #754 sweep */) {
    let canon = |m: &Molecule| chematic_smiles::canonical_smiles(m);
    // Kekulé-written aromatic reactants are matched as their aromatic
    // form, like RDKit's sanitized input: an aliphatic C template does
    // not match benzene, pyranone or tropone.
    for reactant in ["C1=CC=CC=C1", "O=C1C=COC=C1", "O=C1C=CC=CC=C1"] {
        let mol = parse(reactant).unwrap();
        assert!(
            run_reactants("[C:1]>>[C:1]", &[&mol]).unwrap().is_empty(),
            "{reactant}"
        );
        assert!(
            !run_reactants("[c:1]>>[c:1]", &[&mol]).unwrap().is_empty(),
            "{reactant}"
        );
    }
    // Breaking a bond of an explicit-H reactant refills both ends: no
    // radical on the former spiro atom.
    let spiro = chematic_chem::add_hydrogens(&parse("C1CCC2(C1)CCCCC2").unwrap());
    let out = run_reactants("[C:1][C:2]>>[C:1].[C:2]", &[&spiro]).unwrap();
    assert!(!out.is_empty());
    for set in &out {
        for p in set {
            assert!(
                p.atoms().all(|(i, a)| a.element.atomic_number() != 6
                    || chematic_core::implicit_hcount(p, i) + p.neighbors(i).count() as u8 == 4),
                "radical in {}",
                canon(p)
            );
        }
    }
    // A cationic aromatic carbon without H kekulizes, so RDKit's
    // sanitize keeps `[c:1]>>[c+:1]` products.
    let benzene = parse("c1ccccc1").unwrap();
    let out = run_reactants("[c:1]>>[c+:1]", &[&benzene]).unwrap();
    assert_eq!(out.len(), 6);
    for set in &out {
        assert_eq!(set[0].atoms().filter(|(_, a)| a.charge == 1).count(), 1);
        assert_eq!(
            set[0]
                .atoms()
                .map(|(i, _)| u32::from(chematic_core::implicit_hcount(&set[0], i)))
                .sum::<u32>(),
            5
        );
    }
    let naphthalene = parse("c1ccc2ccccc2c1").unwrap();
    assert_eq!(
        run_reactants("[c;H1:1]>>[c+:1]", &[&naphthalene])
            .unwrap()
            .len(),
        8
    );
    assert_eq!(
        run_reactants("[c:1]>>[c-:1]", &[&benzene]).unwrap().len(),
        6
    );
    // Lowering the bond order of H-free atoms of an explicit-H reactant
    // refills them too: RDKit gives `CC(O)O`.
    let acid = chematic_chem::add_hydrogens(&parse("CC(=O)O").unwrap());
    let out = run_reactants("[C:1]=[O:2]>>[C:1]-[O:2]", &[&acid]).unwrap();
    assert_eq!(out.len(), 1);
    assert_eq!(
        canon(&chematic_chem::remove_hydrogens(&out[0][0])),
        canon(&parse("CC(O)O").unwrap())
    );
    // RDKit keeps `[n:1]>>[n+2:1]` on pyrrole as `C1=C[NH+2]C=C1`, but its
    // own SMILES parser rejects the aromatic spelling `c1cc[nH+2]c1`;
    // chematic gives no product (a typed refusal in the checked API).
    let pyrrole = parse("c1cc[nH]c1").unwrap();
    assert!(
        run_reactants("[n:1]>>[n+2:1]", &[&pyrrole])
            .unwrap()
            .is_empty()
    );
    // An explicit-H pyrrole N keeps its H through an identity edit.
    let pyrrole_h = chematic_chem::add_hydrogens(&parse("c1cc[nH]c1").unwrap());
    let out = run_reactants("[n:1]>>[n:1]", &[&pyrrole_h]).unwrap();
    assert_eq!(out.len(), 1);
    assert_eq!(
        canon(&chematic_chem::remove_hydrogens(&out[0][0])),
        canon(&parse("c1cc[nH]c1").unwrap())
    );
    // A charged H-free atom takes RDKit's H count: `C[OH+]C`.
    let ether = chematic_chem::add_hydrogens(&parse("COC").unwrap());
    let out = run_reactants("[O:1]>>[O+:1]", &[&ether]).unwrap();
    assert_eq!(
        canon(&chematic_chem::remove_hydrogens(&out[0][0])),
        canon(&parse("C[OH+]C").unwrap())
    );
}

/// Product sets under the RDKit 2026.03.6 profile, as canonical SMILES
/// (explicit-H inputs are spelled with `[H]` atoms).
fn rdkit_profile_sets(smirks: &str, reactant: &str) -> Vec<Vec<String>> {
    let mol = parse(reactant).unwrap();
    let canon =
        |m: &Molecule| chematic_smiles::canonical_smiles(&chematic_chem::remove_hydrogens(m));
    match run_reactants_traced_rdkit_2026_03_6(smirks, &[&mol], &ReactionTransformLimits::default())
        .unwrap()
    {
        RdkitProfileOutcome::Report(report) => {
            let mut sets: Vec<Vec<String>> = report
                .products
                .iter()
                .map(|set| set.iter().map(|p| canon(&p.molecule)).collect())
                .collect();
            sets.sort();
            sets.dedup();
            sets
        }
        RdkitProfileOutcome::Unsupported(u) => panic!("{smirks}: unsupported {u:?}"),
    }
}

fn canon_of(smiles: &str) -> String {
    chematic_smiles::canonical_smiles(&parse(smiles).unwrap())
}

#[test]
fn biotransformer_corpus_follow_ups_follow_rdkit(/* #734 follow-ups */) {
    // Each row: SMIRKS, reactant (explicit-H rows spelled with `[H]`
    // atoms, as RDKit's AddHs writes them) and RDKit 2026.03.6's
    // sanitized product sets. Covers the ring model for `R<n>`
    // (quinuclidine N is in three rings for RDKit), dummy product atoms,
    // E/Z of created and carried double bonds, an aromatic product atom
    // losing its tag, the H count of a ring CH beside a rewritten C=C,
    // an indole [nH] beside a zero-order bond, and multi-digit charges.
    let cases: &[(&str, &str, &[&[&str]])] = &[
        (
            r#"[#6:1][#7;A;X3;R1,R2;+0:2]>>[#6:1][#7;A;X4+:2][#8-]"#,
            r#"C1CN2CCC1CC2"#,
            &[],
        ),
        (
            r#"[#6:1][#7;A;X3;R3;+0:2]>>[#6:1][#7;A;X4+:2][#8-]"#,
            r#"C1CN2CCC1CC2"#,
            &[&[r#"[O-][N+]12CCC(CC1)CC2"#]],
        ),
        (
            r#"[C:1][O:2]>>[C:1][O:2]P([!#1!#6;O,$([O-])])([!#1!#6;O,$([O-])])=O"#,
            r#"CO"#,
            &[&[r#"*P(*)(=O)OC"#]],
        ),
        (r#"[C:1]>>[C:1]*"#, r#"C"#, &[&[r#"*C"#]]),
        (
            r#"[C:1][O:2]>>[C:1][O:2][!#1;-]"#,
            r#"CO"#,
            &[&[r#"[*-]OC"#]],
        ),
        (
            r#"[C:2]#[C:1]>>[H]/[#6:2]=[#6:1]/[H]"#,
            r#"C/C=C/C#CC"#,
            &[&[r#"C/C=C/C=C/C"#]],
        ),
        (
            r#"[C:1][CH2:3][OH:4]>>[C:1]"#,
            r#"CC/C(CO)=C/C"#,
            &[&[r#"C/C=C/CC"#]],
        ),
        (
            r#"[C:1][CH2:3][OH:4]>>[C:1]O[CH2:3][OH:4]"#,
            r#"CC/C(CO)=C/C"#,
            &[&[r#"CC=C(CC)OCO"#]],
        ),
        (
            r#"[C:1][CH2:3][OH:4]>>[C:1].[CH2:3][OH:4]"#,
            r#"CC/C(CO)=C/C"#,
            &[&[r#"C/C=C/CC"#, r#"[CH2]O"#]],
        ),
        (
            r#"[C:1][CH2:3][CH3:4]>>[C:1]O[CH2:3][CH3:4]"#,
            r#"CC/C(CO)=C/C"#,
            &[&[r#"C/C=C(\CO)OCC"#]],
        ),
        (
            r#"[C:4][C:1]=[C:2][CH2:3]>>[C:4][C:1]=[C:2][CH:3]O"#,
            r#"CC/C=C/CC"#,
            &[&[r#"CCC=CC(C)O"#]],
        ),
        (
            r#"[C:1]=[C:2]/[CH2:3]>>[C:1]=[C:2]/[CH:3]O"#,
            r#"CC/C=C/CC"#,
            &[&[r#"CC/C=C/C(C)O"#]],
        ),
        (
            r#"[#6:1][#16;D2:2][#6,#7:3]>>[#6:1]-,:[S:2](-,:[#6,#7:3])=O"#,
            r#"Cn1c2c(s/c1=N\C(=O)C)CCC2"#,
            &[&[r#"CC(=O)/N=C1/N(C)C2=C(CCC2)S1=O"#]],
        ),
        (
            r#"[O;X1:1]=[#6:2]-1[#6;A;X4;H1,H2:3][#6;A;X4;H1,H2:4][#6:5]-2-[#6:6]-3-[#6:7]-[#6:8]-[#6:9]-4-[#6:10]-[#6:11]-[#6:12]-[#6:13]-4-[#6:14]-3-[#6:15]-[#6:16]-[#6:17]-2=[#6:18]-1>>[#8;H1X2:1]-[c:2]1[c:3][c:4][c:5][c:17](-[#6:16]-[#6:15]-[#6:14]-2-[#6:13]-3-[#6:12]-[#6:11]-[#6:10]-[#6:9]-3-[#6:8]-[#6:7]-[#6:6]-2=O)[c:18]1"#,
            r#"C[C@]12CCC3C(CCC4=CC(=O)CC[C@]43C3CS3)C1CCC2=O"#,
            &[&[r#"C[C@]12CCC(=O)C(CCc3cc(O)ccc3C3CS3)C1CCC2=O"#]],
        ),
        (
            r#"[#6:7][#6;X3:5](=,:[#6:4]!@-[#6;X3:2](-[#6:1])=[O;X1:3])-[#8;X2:6][H:8]>>[H][#8]-[#6;X3:2](-[#6:1])=[O;X1:3].[#6:7][#6;X3:5](=,:[#6:4][H])-[#8;X2:6][H:8]"#,
            r#"[H]Oc1c([H])c(C([H])([H])[H])oc(=O)c1C(=O)C([H])([H])[H]"#,
            &[&[r#"CC(=O)O"#, r#"Cc1cc(O)cc(=O)o1"#]],
        ),
        (
            r#"[#6:1]@[#6;R2:2]1=,:[#6;R2:3](@[#6,#8,#7,#16:4])[#6;R1:5]([H:6])=,:[#6;R1:7](-[!#8;*,#1:11])[#6;R1:8]=,:[#6;R1:9]1[H:10]>>[#6:1]@[#6;R2:2]1=,:[#6;R2:3](@[#6,#8,#7,#16:4])[#6;R1:5](-[#8][H:6])=,:[#6;R1:7](-[!#8;*,#1:11])[#6;R1:8]=,:[#6;R1:9]1[H:10]"#,
            r#"[H]c1c([H])c([H])c2c(c1[H])c(C([H])([H])[H])c([H])n2[H]"#,
            &[],
        ),
        (r#"[C:1]>>[C;+10:1]"#, r#"CC"#, &[&[r#"[C+10]C"#]]),
        (r#"[C:1]>>[C;H15:1]"#, r#"CC"#, &[]),
    ];
    for (smirks, reactant, want) in cases {
        let mut want: Vec<Vec<String>> = want
            .iter()
            .map(|set| set.iter().map(|p| canon_of(p)).collect())
            .collect();
        want.sort();
        assert_eq!(
            rdkit_profile_sets(smirks, reactant),
            want,
            "{smirks} on {reactant}"
        );
    }
    // Surrounding whitespace is not part of the template.
    assert_eq!(
        rdkit_profile_sets(" [C:1][CH2:3][OH:4]>>[C:1] ", r"CC/C(CO)=C/C"),
        vec![vec![canon_of(r"C/C=C/CC")]]
    );
}

#[test]
fn smirks_fuzz_probes_follow_rdkit(/* #754: scripts/validation/smirks_property_fuzz.py */) {
    // An element change takes the template atom's charge, not the
    // matched atom's: RDKit gives a neutral sulfur here, not `[SH+]`.
    assert_eq!(
        rdkit_profile_sets("[N:1]>>[S:1]", "C[N+](C)(C)C"),
        vec![vec![canon_of("CS(C)(C)C")]]
    );
    assert_eq!(
        rdkit_profile_sets("[N:1]>>[S:1]", "CN"),
        vec![vec![canon_of("CS")]]
    );
    // A ring carbon rewritten to sulfur keeps the reactant's aromatic
    // flag in RDKit; three-connected, it takes a double bond in the
    // kekulization (RDKit writes the product `ClS1=CC=CC=N1`) instead
    // of failing it and the product being refused.
    assert_eq!(
        rdkit_profile_sets("[#6&D3+0:1]>>[S:1]", "Clc1ccccn1"),
        vec![vec![canon_of("Cls1ncccc1")]]
    );
    // An explicit-H reactant gives the implicit-H products when the
    // edit changes an element next to an aromatic bond (was `[SH]`, a
    // radical).
    let smirks = "[#6:1]([#6:2]):[cH1:3]>>[*:1](-[S:2])-[*:3]";
    let mut implicit = rdkit_profile_sets(smirks, "Clc1ccccn1");
    let mut explicit = rdkit_profile_sets(smirks, "[H]c1nc(Cl)c([H])c([H])c1[H]");
    implicit.sort();
    explicit.sort();
    implicit.dedup();
    explicit.dedup();
    assert!(!implicit.is_empty());
    assert_eq!(explicit, implicit);

    // Deleting an explicitly hydrogenated atom can lower the valence of
    // a mapped core reached through a ring. The RDKit profile must refill
    // that core instead of emitting a carbon radical (v1.0.42 residual).
    let explicit_reactant =
        chematic_chem::add_hydrogens(&parse("CC(=O)C1=C(O)CN(C2OCCc3ccccc32)C1=O").unwrap());
    let explicit_reactant = chematic_smiles::write(&explicit_reactant);
    let mut products = rdkit_profile_sets(
        "[#6:1](-[OX2:2])-;@[CX4&+0:3]>>[*:1]-[*:3]",
        &explicit_reactant,
    );
    products.sort();
    let mut expected = vec![
        vec![canon_of("CC(=O)C1=CCN(C2OCCc3ccccc32)C1=O")],
        vec![canon_of("CCc1ccccc1CN1CC(O)=C(C(C)=O)C1=O")],
    ];
    expected.sort();
    assert_eq!(products, expected);
}

#[test]
fn biotransformer_rule_shapes_follow_rdkit(/* #734: BioTransformer rule corpus */) {
    // More `[#6:n]` atoms than the aromatic/aliphatic expansion can
    // enumerate (2^9 variants): one SMARTS reading instead of a refusal.
    let many = "[#6:1]1[#6:2][#6:3][#6:4][#6:5][#6:6]1[#6:7][#6:8][#6:9]>>\
                    [#6:1]1[#6:2][#6:3][#6:4][#6:5][#6:6]1[#6:7][#6:8][#6:9]O";
    // RDKit's product atoms take the template's (aliphatic) aromaticity
    // and its implicit product bonds are single, so the ring the template
    // spells is saturated; a ring bond the template leaves alone stays
    // aromatic after RDKit's sanitize.
    assert_eq!(
        rdkit_profile_sets(many, "CCCc1ccccc1"),
        vec![vec![canon_of("OCCCC1CCCCC1")]]
    );
    assert_eq!(
        rdkit_profile_sets(
            "[#6:1]1[#6:2][#6:3][#6:4][#6:5][#6:6]1>>[#6:1]1[#6:2][#6:3][#6:4][#6:5][#6:6]1O",
            "c1ccccc1"
        ),
        vec![vec![canon_of("OC1CCCCC1")]]
    );
    assert_eq!(
        rdkit_profile_sets("[#6:1][#6:2]>>[#6:1][#6:2]O", "c1ccccc1"),
        vec![vec![canon_of("Oc1ccccc1")]]
    );
    // Product bond expressions take RDKit's bond type: the first order
    // primitive, or no order (a zero-order `~` bond) after `@`/`!@`; a
    // lone `~` keeps the reactant bond.
    for (smirks, reactant, want) in [
        ("[C:1][O:2]>>[C:1]-;!@[O:2]", "CO", "CO"),
        ("[C:1][O:2]>>[C:1]!=[O:2]", "CO", "C=O"),
        ("[C:1][O:2]>>[C:1]~[O:2]", "CO", "CO"),
    ] {
        assert_eq!(
            rdkit_profile_sets(smirks, reactant),
            vec![vec![canon_of(want)]],
            "{smirks}"
        );
    }
    for smirks in ["[C:1][O:2]>>[C:1]!@-[O:2]", "[C:1][O:2]>>[C:1]@[O:2]"] {
        let mol = parse("CO").unwrap();
        let RdkitProfileOutcome::Report(report) = run_reactants_traced_rdkit_2026_03_6(
            smirks,
            &[&mol],
            &ReactionTransformLimits::default(),
        )
        .unwrap() else {
            panic!("{smirks}");
        };
        let product = &report.products[0][0].molecule;
        // C~O with CH4 and OH2: the bond adds nothing to valence.
        assert_eq!(product.bond(BondIdx(0)).order, BondOrder::Zero, "{smirks}");
        assert_eq!(
            chematic_core::implicit_hcount(product, AtomIdx(0)),
            4,
            "{smirks}"
        );
    }
    let new_bond = rdkit_profile_sets("[N+:1]>>[N+:1]~C", "C[NH3+]");
    assert_eq!(new_bond.len(), 1);
    assert!(new_bond[0][0].contains("NH3+"), "{new_bond:?}");
    // `/` `\` in a SMARTS-only reactant template match like plain SMARTS.
    assert_eq!(
        rdkit_profile_sets(
            "[#6;X3:1]=[#6:2](/[#6:3])>>[#6:1]1[#6:2]([#6:3])O1",
            "CC=CC"
        )
        .len(),
        1
    );
    // A product atom naming no single element keeps the matched one;
    // H count and charge spelled in a list take the first alternative.
    for (smirks, reactant, want) in [
        (
            "[C:1][Cl:2]>>[C:1].[F,Cl,Br,I;-:2]",
            "CCl",
            vec!["C", "[Cl-]"],
        ),
        ("[N:1]>>[N;H1,H2:1]", "CN", vec!["C[NH]"]),
        ("[N:1]>>[N;X4H2+,X3H1+0:1]", "CN", vec!["C[NH2+]"]),
    ] {
        let mut want: Vec<String> = want.into_iter().map(canon_of).collect();
        want.sort();
        let mut got = rdkit_profile_sets(smirks, reactant);
        got.iter_mut().for_each(|s| s.sort());
        assert_eq!(got, vec![want], "{smirks}");
    }
    // E/Z next to the reaction survives explicit H atoms: on the
    // template's own single bond, and when the H that carried it folds.
    for (smirks, reactant, want) in [
        (
            "[#6:1]-[#6;X3:2](-[#6:4])=[O;X1:3]>>[H][#8;X2:3][C;X4:2]([H])([#6:1])[#6:4]",
            "O=C(C)/C(C)=C/C",
            "CC(O)/C(C)=C/C",
        ),
        ("[C:1]=[C:2]>>[C:1]=[C:2]", "[H]/C(C)=C(/[H])CC", "C/C=C/CC"),
    ] {
        let mol = chematic_chem::add_hydrogens(&parse(reactant).unwrap());
        let input = chematic_smiles::write(&mol);
        assert_eq!(
            rdkit_profile_sets(smirks, &input),
            vec![vec![canon_of(want)]],
            "{smirks}"
        );
    }
    // A template that writes Kekulé bonds into an aromatic ring is judged
    // after RDKit's sanitize: kekulize, perceive aromaticity again.
    let demethylation = "[H:1][*:2][#6:3](=,:[*:4][H:5])-[#8;X2R0:6][C:7]([H:8])([H])[H:9]>>\
                             [H:8][#6:7]([H:9])=O.[H:1][*:2][#6:3](=,:[*:4][H:5])-[#8;X2R0:6][H]";
    let anisole =
        chematic_smiles::write(&chematic_chem::add_hydrogens(&parse("COc1ccccc1").unwrap()));
    let mut got = rdkit_profile_sets(demethylation, &anisole);
    got.iter_mut().for_each(|s| s.sort());
    let mut want = vec![canon_of("C=O"), canon_of("Oc1ccccc1")];
    want.sort();
    assert_eq!(got, vec![want]);
    // `$` is a quadruple bond on the reactant side too (it matched any
    // bond when the SMARTS reader did not know it).
    assert!(rdkit_profile_sets("[C:1]$[O:2]>>[C:1].[O:2]", "CCO").is_empty());
    assert!(rdkit_profile_sets("[C:1]$[C:2]>>[C:1]-[C:2]", "N#CC").is_empty());
    // A double bond between two ring atoms that stay aromatic satisfies
    // neither: RDKit's cumulene for benzene.
    assert_eq!(
        rdkit_profile_sets("[c:1]:[c:2]>>[c:1]=[c:2]", "c1ccccc1"),
        vec![vec![canon_of("C1=C=CC=CC=1")]]
    );
    // A Kekulé ring spelled over one ring of naphthalene keeps it
    // naphthalene, its fusion bond single or double.
    for product in [
        "[#6:1]1=[#6:2][#6:3]=[#6:4][#6:5]=[#6:6]1",
        "[#6:1]1-[#6:2]=[#6:3]-[#6:4]=[#6:5]-[#6:6]=1",
    ] {
        let smirks = format!("[#6:1]1:[#6:2]:[#6:3]:[#6:4]:[#6:5]:[#6:6]:1>>{product}");
        let sets = rdkit_profile_sets(&smirks, "c1ccc2ccccc2c1");
        assert!(!sets.is_empty(), "{smirks}");
        assert!(
            sets.iter()
                .all(|set| set == &vec![canon_of("c1ccc2ccccc2c1")]),
            "{smirks}: {sets:?}"
        );
    }
    // E/Z of an exocyclic imine whose ring-side direction sits on an
    // aromatic ring bond survives a template elsewhere (BTMR0884), in
    // both profiles.
    let imine = "ClCC/N=c1\\c(O)c(O)\\c1=N/C";
    assert_eq!(
        rdkit_profile_sets("[Cl:1]>>[Cl:1]", imine),
        vec![vec![canon_of(imine)]]
    );
    let mol = parse(imine).unwrap();
    let native = run_reactants("[Cl:1]>>[Cl:1]", &[&mol]).unwrap();
    assert_eq!(
        chematic_smiles::canonical_smiles(&native[0][0]),
        canon_of(imine)
    );
    // An explicit-H ring the template dearomatizes takes RDKit's H
    // counts: the carbon that loses its `[H]` and the ring CH beside it
    // are refilled (BTMR0655 on toluene).
    let ring_oh = "[H][#6;a:5]1[#6,#7;a:4][c:3][c:2][c:1][c:6]1!@-[#6:7]>>\
                       [H][#8]-[#6;a:5]1[#6,#7;a:4][c:3][c:2][c:1][c:6]1!@-[#6:7]";
    let sets = rdkit_profile_sets(ring_oh, "[H]c1c([H])c([H])c(C([H])([H])[H])c([H])c1[H]");
    let reread: Vec<Vec<String>> = sets
        .iter()
        .map(|set| set.iter().map(|smi| canon_of(smi)).collect())
        .collect();
    assert_eq!(reread, vec![vec![canon_of("C~C1=CC=CCC1O")]]);
    // So is a product `[#6;a:1]`: RDKit's dioxygenation BTMR0842 leaves
    // one aliphatic atom in a ring of aromatic ones, which RDKit's
    // sanitize cannot kekulize; the native reading gives the catechol.
    let dioxygenation = "[#8:7]([c:6]1[c:1][c:2][c:3][c:4][c:5]1)[c:8]1[c:10][c:11][c:12]\
                             [c:13][#6;a;H1X3:9]1>>[#8:7]-[c:6]1[c:1][c:2][c:3][c:4][c:5]1.\
                             [#8]-[c:8]1[c:10][c:11][c:12][c:13][#6;a:9]1-[#8]";
    assert!(rdkit_profile_sets(dioxygenation, "c1ccccc1Oc1ccccc1").is_empty());
    let mol = parse("c1ccccc1Oc1ccccc1").unwrap();
    assert!(!run_reactants(dioxygenation, &[&mol]).unwrap().is_empty());
    // A product `[#6:1]` mapped to a reactant `[c:1]` is aliphatic, as
    // RDKit builds it (BTMR0826 ring cleavage leaves no aromatic atom).
    let cleavage = "[#6:8]-[c:2]1[#6;a;H1X3:3][c:4]([#8;A;H1X2:9])[#6;a;H1X3:5]\
                        [#6;a;H1X3:6][c:1]1[#8;A;H1X2:7]>>[#6:8]-[#6:2](=O)[#6;A;H1X3:3]\
                        [#6:4](=[#8;A;X1:9])-[#6;A;H1X3:5]=[#6;A;H1X3:6]-[#6:1]([#8;A;X1-:7])=O";
    assert_eq!(
        rdkit_profile_sets(cleavage, "CC(C)(C)c1cc(O)ccc1O"),
        vec![vec![canon_of("CC(C)(C)C(=O)[CH]C(=O)C=CC(=O)[O-]")]]
    );
    // A ring-fusion N next to a `c(=O)` carbon takes no double bond
    // (BTMR0678 methyl hydroxylation of an imidazo[1,2-a]pyrimidinone).
    let methyl = "[#6;A;H3X4:1][c:2]>>[#8][#6;A;H2X4:1][c:2]";
    assert_eq!(
        rdkit_profile_sets(methyl, "Cc1cn2c(=O)ccnc2n1C"),
        vec![vec![canon_of("OCc1cn2c(=O)ccnc2n1C")]]
    );
    // S-oxidation that writes a Kekulé thiophene into benzothiophene: the
    // `[s+]` with three bonds needs no double bond (BTMR0102).
    let s_oxidation = "[#6:1]=,:1[#6:2]=,:[#6:3][#16;X2:4][#6:5]=,:1>>\
                           [#8-][S;X3+:4]-,:1-,:[#6:5]=,:[#6:1][#6:2]=,:[#6:3]-,:1";
    let sets = rdkit_profile_sets(s_oxidation, "c1ccc2sccc2c1");
    assert!(!sets.is_empty());
    assert!(
        sets.iter()
            .all(|set| set == &vec![canon_of("[O-][s+]1ccc2ccccc21")]),
        "{sets:?}"
    );
}

#[test]
fn product_template_policies_follow_rdkit(/* issue #754 policy flags */) {
    let canon = |m: &Molecule| chematic_smiles::canonical_smiles(m);
    let sets = |smirks: &str, reactant: &str| -> Vec<Vec<String>> {
        let mol = parse(reactant).unwrap();
        run_reactants(smirks, &[&mol])
            .unwrap()
            .iter()
            .map(|set| set.iter().map(canon).collect())
            .collect()
    };
    let c2h4 = canon(&parse("C=C").unwrap());
    // A product bond list takes its first alternative.
    for smirks in [
        "[C:1]-[C:2]>>[C:1]=,:[C:2]",
        "[#6:1]-,:[#6:2]>>[#6:1]=,:[#6:2]",
    ] {
        let out = sets(smirks, "CC");
        assert!(!out.is_empty(), "{smirks}");
        assert!(out.iter().all(|set| set == &vec![c2h4.clone()]), "{smirks}");
    }
    // A grouped product component is one product object.
    let out = sets("[C:1][O:2]>>([C:1].[O:2]).[N]", "CCO");
    assert_eq!(out.len(), 1);
    assert_eq!(out[0].len(), 2);
    assert_eq!(out[0][0], canon(&parse("CC.O").unwrap()));
    assert!(run_reactants("[C:1][O:2]>>([C:1].[O:2]", &[&parse("CCO").unwrap()]).is_err());
    // An element change re-derives the H count, as RDKit does.
    let mut out = sets("[C:1]>>[N:1]", "[13CH3]C");
    out.sort();
    assert_eq!(
        out,
        vec![
            vec![canon(&parse("CN").unwrap())],
            vec![canon(&parse("[13CH3]N").unwrap())]
        ]
    );
    assert_eq!(
        sets("[C:1]>>[N:1]", "[CH2]C")[0],
        vec![canon(&parse("CN").unwrap())]
    );
    // Two-digit product charges and H counts parse, as in RDKit: `+10`
    // is kept, `H15` fails the product filter.
    assert_eq!(
        sets("[C:1]>>[C;+10:1]", "CC")[0],
        vec![canon(&parse("C[C+10]").unwrap())]
    );
    assert!(sets("[C:1]>>[C;H15:1]", "CC").is_empty());
    // An aromatic ring bond rewritten as `=` between aromatic atoms is
    // dropped, as RDKit's sanitize drops it for pyridine.
    assert!(sets("[#6:1]~[#7:2]>>[#6:1]=[#7:2]", "c1ccncc1").is_empty());
    // A quaternized pyrrole N kekulizes (RDKit: `C[N+]1(C)C=CC=C1`).
    assert_eq!(sets("[#7;a:1]>>[#7+:1]C", "Cn1cccc1").len(), 1);
    // `:` flags the bond's atoms aromatic: outside a ring that fails
    // RDKit's sanitize, so no product survives.
    assert!(sets("[#6:1]-[#6:2]>>[#6:1]:[#6:2]", "CC").is_empty());
}

#[test]
fn bug_check_products_match_rdkit(/* issue #734 bug check */) {
    // Each row: SMIRKS, reactant, RDKit 2026.03.6's sanitized products
    // (empty: no product survives sanitize).
    let canon = |s: &str| chematic_smiles::canonical_smiles(&parse(s).unwrap());
    for (smirks, reactant, expected) in [
        // A dative `->` in a template is a bond, not a reaction arrow.
        ("[N:1]->[Pt:2]>>[N:1].[Pt:2]", "CN->[Pt]", vec!["CN.[Pt]"]),
        // RDKit's clean-up spellings are valid products.
        (
            "[C:1][Br:2]>>[C:1][Cl:2]",
            "BrCCN(=O)=O",
            // The RDKit-parity view of the reactant carries RDKit's
            // nitro clean-up, as RDKit's reactant does.
            vec!["ClCC[N+](=O)[O-]"],
        ),
        (
            "[c:1][Cl:2]>>[c:1]N(=O)=O",
            "Clc1ccccc1",
            vec!["O=N(=O)c1ccccc1"],
        ),
        ("[O:1]>>[O:1]Cl(=O)(=O)=O", "CO", vec!["COCl(=O)(=O)=O"]),
        ("[C:1][Br:2]>>[C:1]N=N#N", "CBr", vec!["CN=N#N"]),
        // RDKit's valence model: Mg is unrestricted, Cl+ is S-like.
        ("[C:1]>>[Mg:1](C)(C)C", "CC", vec!["C[Mg](C)(C)C"]),
        ("[Cl:1]>>[Cl+:1](C)C", "CCl", vec!["C[ClH+](C)C"]),
        ("[S:1]>>[S+:1](C)(C)(C)(C)C", "CS", vec![]),
        // Unchanged template degree keeps the reactant's explicit H.
        ("[n:1]>>[n:1]", "c1cc[nH]c1", vec!["c1cc[nH]c1"]),
        ("[n:1]>>[n:1]C", "c1cc[nH]c1", vec!["Cn1cccc1"]),
        // RDKit implicit hydrogens beyond the organic subset.
        ("[C:1]>>[Se:1]", "CC", vec!["C[SeH]"]),
        ("[C:1]>>[Si:1]", "CC", vec!["C[SiH3]"]),
        // SMARTS reactant semantics.
        ("[*:1]>>[*:1]O", "CN", vec!["CNO", "NCO"]),
        ("[#7;H0:1]>>[#7:1]C", "CNC", vec![]),
        ("[NH0:1]>>[N:1]C", "CNC", vec![]),
        ("[N+0:1]>>[N:1]C", "C[NH3+]", vec![]),
        // Aromatic atoms outside rings are rejected.
        ("[C:1]>>[c:1]", "CC", vec![]),
        // Unspecified product charge and isotope keep the reactant's;
        // an explicit H0 applies.
        ("[O-:1]>>[O:1]", "C[O-]", vec!["C[O-]"]),
        ("[N:1]>>[N:1]C", "C[NH3+]", vec!["C[NH2+]C"]),
        ("[C:1]>>[13C:1]", "CC", vec!["C[13CH3]"]),
        ("[C:1]>>[CH0:1]", "CC", vec!["[C]C"]),
        // RDKit 2026.03.6 treats H0 on a newly created, unmapped atom
        // as unspecified, unlike H0 on a mapped product atom above.
        ("[C:1]>>[C:1][CH0]", "C", vec!["CC"]),
        ("[C:1]>>[C:1][CH0:2]", "C", vec!["CC"]),
        ("[C:1]>>[C:1][OH0]", "C", vec!["CO"]),
        ("[C:1]>>[C:1][NH0]", "C", vec!["CN"]),
    ] {
        let mol = parse(reactant).unwrap();
        let mut got: Vec<String> = run_reactants(smirks, &[&mol])
            .unwrap_or_else(|e| panic!("{smirks}: {e}"))
            .iter()
            .map(|set| {
                let mut parts: Vec<String> = set.iter().map(canonical).collect();
                parts.sort();
                parts.join(".")
            })
            .collect();
        got.sort();
        got.dedup();
        let mut want: Vec<String> = expected
            .iter()
            .map(|e| {
                let mut parts: Vec<String> = e.split('.').map(canon).collect();
                parts.sort();
                parts.join(".")
            })
            .collect();
        want.sort();
        assert_eq!(got, want, "{smirks} on {reactant}");
    }
    // Atom-map positions read only the reactant side, so a product
    // atomic-number primitive is fine.
    let mol = parse("CN").unwrap();
    let smirks = "[N:1]>>[#7:1]C";
    let matches = find_reaction_matches(smirks, &[&mol]).unwrap();
    assert_eq!(matches[0].atom_map_positions(smirks).unwrap().len(), 1);
}

#[test]
fn products_rdkit_sanitize_rejects_are_dropped(/* issue #734 item 3 */) {
    // RDKit 2026.03.6 rejects each of these raw products at sanitize; the
    // native model kept the last three (N valence 3 or 5, aromatic bond
    // counted as one).
    for (smirks, reactant) in [
        ("[C:1][O:2]>>[C:1]=[O:2]", "COC"),
        ("[C:1][N:2]>>[C:1]=[N:2]", "CN(C)C"),
        ("[N:1]>>[N:1]=O", "CNC"),
    ] {
        let mol = parse(reactant).unwrap();
        assert!(
            run_reactants(smirks, &[&mol]).unwrap().is_empty(),
            "{smirks} on {reactant}"
        );
    }
    // Ipso attack on toluene's methyl-bearing carbon is dropped; the
    // ortho/meta/para products remain (RDKit: the same three).
    let canon = |m: &Molecule| chematic_smiles::canonical_smiles(m);
    let toluene = parse("Cc1ccccc1").unwrap();
    let mut got: Vec<String> = run_reactants("[c:1]>>[c:1]O", &[&toluene])
        .unwrap()
        .iter()
        .map(|set| canon(&set[0]))
        .collect();
    got.sort();
    got.dedup();
    let mut want: Vec<String> = ["Cc1ccccc1O", "Cc1cccc(O)c1", "Cc1ccc(O)cc1"]
        .iter()
        .map(|s| canon(&parse(s).unwrap()))
        .collect();
    want.sort();
    assert_eq!(got, want);
    // The charged spellings are valid and kept.
    let tma = parse("CN(C)C").unwrap();
    assert!(
        !run_reactants("[C:1][N:2]>>[C:1]=[N+:2]", &[&tma])
            .unwrap()
            .is_empty()
    );
    assert!(
        !run_reactants("[N:1]>>[N+:1][O-]", &[&tma])
            .unwrap()
            .is_empty()
    );
}

#[test]
fn reactant_smarts_query_selects_the_matching_atom(/* issue #734 */) {
    // A reactant query that only some atoms satisfy must transform only
    // those atoms, proving the query is actually evaluated (not ignored).
    // Primary alcohol carbon oxidation: only the CH2 next to OH matches.
    let canon = |m: &Molecule| chematic_smiles::canonical_smiles(m);
    let ethanol = parse("CCO").unwrap();
    let products = run_reactants("[CX4;H2:1][OX2H:2]>>[C:1]=[O:2]", &[&ethanol]).unwrap();
    assert_eq!(products.len(), 1);
    assert_eq!(canon(&products[0][0]), canon(&parse("CC=O").unwrap()));
    // The isopropanol carbon is CH, so the H2 query must not match it.
    let isopropanol = parse("CC(C)O").unwrap();
    assert!(
        run_reactants("[CX4;H2:1][OX2H:2]>>[C:1]=[O:2]", &[&isopropanol])
            .unwrap()
            .is_empty()
    );
}

#[test]
fn reactant_smarts_only_stereo_is_a_typed_error(/* issue #734 */) {
    // A SMARTS-only reactant template (one that is not valid SMILES) that
    // carries tetrahedral or E/Z stereo is refused with a typed error
    // rather than silently dropping the constraint.
    let subject = parse("F[C@H](Cl)Br").unwrap();
    assert!(matches!(
        run_reactants("[C@;X4:1](F)(Cl)Br>>[C:1](F)(Cl)I", &[&subject]),
        Err(TransformError::SmirksParse(
            crate::reaction::RxnError::UnsupportedReactantTemplate { .. }
        ))
    ));
    // A purely SMILES reactant template with the same stereo still works
    // (it goes through the SMILES parser and the existing post-check).
    assert!(run_reactants("[C@H:1](F)(Cl)Br>>[C:1](F)(Cl)I", &[&subject]).is_ok());
}

#[test]
fn bracket_hydrogen_atom_parses_as_hydrogen(/* issue #734 */) {
    // `[H]` inside a SMARTS reactant template is a hydrogen atom, not an
    // H-count primitive: it matches an explicit hydrogen neighbour.
    let mol = parse("C").unwrap();
    let mol = chematic_chem::add_hydrogens(&mol);
    let products = run_reactants("[C:1][H]>>[C:1]", &[&mol]).unwrap();
    assert!(!products.is_empty(), "[H] must match an explicit hydrogen");
}

#[test]
fn product_atomic_number_primitives_are_literal_atoms_with_implicit_h() {
    // Issue #679: an unmapped product-side `[#6](=[#8])[#6]` adds an acetyl
    // group exactly like `C(=O)C`; no aromatic-flag variants, no `[C]`
    // radicals, one product.
    let ethanol = parse("CCO").unwrap();
    let canon = |m: &Molecule| chematic_smiles::canonical_smiles(m);
    let organic = run_reactants("[O:1]>>[*:1]C(=O)C", &[&ethanol]).unwrap();
    let atomic = run_reactants("[O:1]>>[*:1][#6](=[#8])[#6]", &[&ethanol]).unwrap();
    assert_eq!(organic.len(), 1);
    assert_eq!(atomic.len(), 1);
    assert_eq!(canon(&atomic[0][0]), canon(&organic[0][0]));
    assert_eq!(canon(&atomic[0][0]), "CC(=O)OCC");
    // Bracket product atoms behave the same: `[C]` is an unspecified
    // hydrogen count in a product template, as for mapped atoms (#18).
    let bracket = run_reactants("[O:1]>>[*:1][C](=[O])[C]", &[&ethanol]).unwrap();
    assert_eq!(canon(&bracket[0][0]), "CC(=O)OCC");
    // `;H1` pins the count; non-organic-subset elements stay bracketed.
    let vinyl = run_reactants("[O:1]>>[*:1][#6;H1]=[#6]", &[&ethanol]).unwrap();
    assert_eq!(canon(&vinyl[0][0]), canon(&parse("C=COCC").unwrap()));
    let silyl = run_reactants("[O:1]>>[*:1][#14]", &[&ethanol]).unwrap();
    assert_eq!(silyl.len(), 1);
    assert_eq!(
        silyl[0][0]
            .atoms()
            .filter(|(_, a)| a.element.atomic_number() == 14)
            .count(),
        1
    );
}

#[test]
fn mapped_product_atomic_number_primitive_keeps_reactant_aromaticity() {
    // `[#7:1]` on both sides shares one spelling per variant, so an
    // aliphatic amine gives one aliphatic product and pyridine one aromatic
    // product (previously each gave an extra product with the flag flipped).
    let canon = |m: &Molecule| chematic_smiles::canonical_smiles(m);
    let amine = parse("CCN").unwrap();
    let out = run_reactants("[#7:1]>>[#7:1]C", &[&amine]).unwrap();
    assert_eq!(out.len(), 1);
    assert_eq!(canon(&out[0][0]), canon(&parse("CCNC").unwrap()));
    // Pyridine needs the charge (a neutral N-methylpyridine has a
    // four-valent N, rejected like RDKit's sanitize does, #734).
    let pyridine = parse("c1ccncc1").unwrap();
    assert!(
        run_reactants("[#7:1]>>[#7:1]C", &[&pyridine])
            .unwrap()
            .is_empty()
    );
    let out = run_reactants("[#7:1]>>[#7+:1]C", &[&pyridine]).unwrap();
    assert_eq!(out.len(), 1);
    assert!(
        out[0][0]
            .atoms()
            .any(|(_, a)| a.aromatic && a.element.atomic_number() == 7)
    );
    // A product map spelled with a plain symbol on the reactant side
    // follows that spelling.
    let variants = crate::reaction::expand_atomic_number_primitives("[n:1]>>[#7:1]C").unwrap();
    assert_eq!(variants, vec!["[n:1]>>[n:1]C".to_string()]);
    let variants = crate::reaction::expand_atomic_number_primitives("[N:1]>>[#7:1]C").unwrap();
    assert_eq!(variants, vec!["[N:1]>>[N:1]C".to_string()]);
}

#[test]
fn prepared_reaction_applies_atomic_number_smirks_variants() {
    let reactant = parse("NC=O").unwrap();
    let prepared = PreparedReaction::new("[#7:1][C:2](=[O:3])>>[#7:1][C:2](=[O:3])").unwrap();
    let products = prepared.run_reactants(&[&reactant]).unwrap();
    // Two variants (`[N:1]`/`[n:1]` on both sides together); only the
    // aliphatic one matches formamide, and the product keeps the reactant's
    // aliphatic nitrogen (issue #679: no aromatic-flag duplicate).
    assert_eq!(products.len(), 1);
    assert!(products.iter().all(|set| set.len() == 1));
    assert!(!products[0][0].atoms().any(|(_, a)| a.aromatic));
}

#[test]
fn prepared_reaction_reports_each_atomic_number_variant() {
    let reactant = parse("NC=O").unwrap();
    let prepared = PreparedReaction::new("[#7:1][C:2](=[O:3])>>[#7:1][C:2](=[O:3])").unwrap();
    let reports = prepared
        .run_reactants_with_variant_diagnostics(&[&reactant], &ReactionTransformLimits::default())
        .unwrap();
    assert_eq!(reports.len(), 2);
    assert_eq!(reports[0].variant_index, 0);
    assert_eq!(reports[1].variant_index, 1);
    assert_eq!(
        reports[0].normalized_smirks,
        "[N:1][C:2](=[O:3])>>[N:1][C:2](=[O:3])"
    );
    assert_eq!(
        reports[1].normalized_smirks,
        "[n:1][C:2](=[O:3])>>[n:1][C:2](=[O:3])"
    );
    assert!(reports[0].diagnostics.accepted_matches > 0);
    assert!(reports[0].diagnostics.applied_products > 0);
    assert!(reports[1].diagnostics.accepted_matches == 0);
}

#[test]
fn atomic_number_smirks_rejects_unsupported_compound_product_primitive() {
    let reactant = parse("N").unwrap();
    // An unmapped product atom naming no single element cannot be built
    // and is rejected (#734); `;X3` alone would just be dropped.
    let err = run_reactants("[#7:1]>>[#7:1][#6,#7]", &[&reactant]);
    assert!(matches!(
        err,
        Err(TransformError::SmirksParse(
            crate::reaction::RxnError::UnsupportedProductPrimitive { .. }
        ))
    ));
    assert!(run_reactants("[#7:1]>>[#7;X3:1]", &[&reactant]).is_ok());
}

#[test]
fn prepared_reaction_is_send_sync_and_matches_legacy_output() {
    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<PreparedReaction>();
    assert_send_sync::<ReactionMatchContext>();

    let mol = parse("CCOC(=O)C").unwrap();
    let smirks = "[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])O.[O:3][C:4]";
    let prepared = PreparedReaction::new(smirks).unwrap();
    let legacy = run_reactants(smirks, &[&mol]).unwrap();
    let reused = prepared.run_reactants(&[&mol]).unwrap();
    let rings = chematic_perception::find_sssr(&mol);
    let reused_with_rings = prepared
        .run_reactants_with_rings(&[&mol], &[&rings])
        .unwrap();
    let canonical_sets = |sets: Vec<Vec<Molecule>>| {
        sets.into_iter()
            .map(|set| set.into_iter().map(|m| canonical(&m)).collect::<Vec<_>>())
            .collect::<Vec<_>>()
    };
    let legacy = canonical_sets(legacy);
    assert_eq!(canonical_sets(reused), legacy);
    assert_eq!(canonical_sets(reused_with_rings), legacy);

    let strict = prepared
        .run_reactants_strict_with_rings(&[&mol], &[&rings])
        .unwrap();
    let strict_limited = prepared
        .run_reactants_strict_with_rings_and_limits(
            &[&mol],
            &[&rings],
            &ReactionTransformLimits::default(),
        )
        .unwrap();
    assert_eq!(canonical_sets(strict), canonical_sets(strict_limited));
    assert_eq!(
        prepared
            .find_matches_with_rings_and_limits(
                &[&mol],
                &[&rings],
                &ReactionTransformLimits::default(),
            )
            .unwrap(),
        prepared
            .find_matches_with_rings(&[&mol], &[&rings])
            .unwrap()
    );
    let context = ReactionMatchContext::new(&mol, Some(&rings));
    assert_eq!(
        prepared.find_matches_with_context(&context).unwrap(),
        prepared
            .find_matches_with_rings(&[context.target()], &[context.rings()])
            .unwrap()
    );
}

#[test]
fn prepared_reaction_validates_once_and_preserves_limit_errors() {
    assert!(PreparedReaction::new("not-smirks").is_err());
    let mol = parse("CC").unwrap();
    let prepared = PreparedReaction::new("[C:1]>>[C:1]").unwrap();
    let err = match prepared
        .run_reactants_with_limits(&[&mol], &ReactionTransformLimits { max_matches: 0 })
    {
        Ok(_) => panic!("zero match limit must reject the first match set"),
        Err(err) => err,
    };
    assert!(matches!(
        err,
        TransformError::ResourceLimit {
            resource: "reaction matches",
            ..
        }
    ));
}

#[test]
fn prepared_find_and_apply_match_legacy_split_api() {
    let mol = parse("NCCN").unwrap();
    let smirks = "[N:1]>>[N:1]";
    let prepared = PreparedReaction::new(smirks).unwrap();

    let legacy_matches = find_reaction_matches(smirks, &[&mol]).unwrap();
    let prepared_matches = prepared.find_matches(&[&mol]).unwrap();
    assert_eq!(prepared_matches, legacy_matches);
    for m in &prepared_matches {
        let legacy = apply_reaction_match(smirks, &[&mol], m, true).unwrap();
        let reused = prepared.apply_match(&[&mol], m, true).unwrap();
        assert_eq!(
            reused.map(canonical_set),
            legacy.map(canonical_set),
            "prepared match application must preserve legacy output"
        );
    }
}

#[test]
fn reaction_transform_limits_reject_matches_before_cartesian_product() {
    let mol = parse("CC").unwrap();
    let err = find_reaction_matches_with_limits(
        "[C:1]>>[C:1]",
        &[&mol],
        &ReactionTransformLimits { max_matches: 0 },
    );
    assert!(matches!(
        err,
        Err(TransformError::ResourceLimit {
            resource: "reaction matches",
            ..
        })
    ));
}

#[test]
fn no_match_returns_empty() {
    let mol = parse("C").unwrap();
    let results = run_reactants("[N:1]>>[N:1]", &[&mol]).unwrap();
    assert!(
        results.is_empty(),
        "nitrogen template must not match methane"
    );
}

#[test]
fn multiple_matches_in_single_mol() {
    let mol = parse("NCCN").unwrap();
    let results = run_reactants("[N:1]>>[N:1]", &[&mol]).unwrap();
    assert_eq!(results.len(), 2, "two N atoms in NCCN → two product sets");
}

#[test]
fn bond_formation_two_mols() {
    let n_mol = parse("N").unwrap();
    let c_mol = parse("C").unwrap();
    let results = run_reactants("[N:1].[C:2]>>[N:1][C:2]", &[&n_mol, &c_mol]).unwrap();
    assert!(!results.is_empty());
    let prod = &results[0][0];
    assert_eq!(prod.atom_count(), 2, "product must have 2 atoms");
    assert_eq!(prod.bonds().count(), 1, "product must have 1 bond");
}

#[test]
fn bond_cleavage_two_products() {
    let mol = parse("CC").unwrap();
    let results = run_reactants("[C:1][C:2]>>[C:1].[C:2]", &[&mol]).unwrap();
    assert!(!results.is_empty());
    let products = &results[0];
    assert_eq!(products.len(), 2, "two product templates → two products");
    assert_eq!(products[0].atom_count(), 1);
    assert_eq!(products[1].atom_count(), 1);
}

#[test]
fn reactant_count_mismatch_error() {
    let mol = parse("C").unwrap();
    let err = run_reactants("[N:1].[C:2]>>[N:1][C:2]", &[&mol]);
    assert!(
        matches!(
            err,
            Err(TransformError::ReactantCountMismatch {
                expected: 2,
                got: 1
            })
        ),
        "two-template SMIRKS with one reactant must error"
    );
}

#[test]
fn invalid_smirks_error() {
    let mol = parse("C").unwrap();
    let err = run_reactants("[X]>>[X]", &[&mol]);
    assert!(
        matches!(err, Err(TransformError::SmirksParse(_))),
        "unknown element must yield SmirksParse error"
    );
}

#[test]
fn overvalent_product_filtered_oxygen() {
    // O normally has max valence 2.
    // SMIRKS adds two carbons to an oxygen that already has one bond → 3 bonds on O → invalid.
    // CCO: the O is bonded to 1 C (bond_sum=1). Template [O:1]>>[O:1](C)C adds 2 more.
    let ethanol = parse("CCO").unwrap();
    let results = run_reactants("[O:1]>>[O:1](C)C", &[&ethanol]).unwrap();
    // The O that already had 1 bond would get 3 → over-valenced → filtered out.
    // The only match is the terminal O (1 bond → +2 = 3 bonds, invalid).
    assert!(
        results.is_empty(),
        "product with O having 3 bonds must be filtered out, got {} sets",
        results.len()
    );
}

#[test]
fn valid_charged_product_kept() {
    // N with charge +1 can have up to 4 bonds (normal valences [3,5], +1 allows 4).
    // trimethylamine N(C)(C)C has N with bond_sum=3, charge=0.
    // Template [N:1]>>[N+:1] just changes charge, keeps 3 bonds → valid.
    let tma = parse("N(C)(C)C").unwrap();
    let results = run_reactants("[N:1]>>[N+:1]", &[&tma]).unwrap();
    assert!(
        !results.is_empty(),
        "N+ with 3 bonds must be valid and kept"
    );
}

#[test]
fn new_atom_in_product() {
    let mol = parse("C").unwrap();
    let results = run_reactants("[C:1]>>[C:1]=O", &[&mol]).unwrap();
    assert!(!results.is_empty());
    let prod = &results[0][0];
    assert_eq!(prod.atom_count(), 2, "C + new O = 2 atoms");
}

#[test]
fn amide_bond_formation() {
    // NH3 + H-C(=O)-Cl → H-C(=O)-NH2 (formamide)
    let nh3 = parse("N").unwrap();
    let hcocl = parse("C(=O)Cl").unwrap();
    let results = run_reactants("[N:1].[C:2](=O)Cl>>[C:2](=O)[N:1]", &[&nh3, &hcocl]).unwrap();
    assert!(!results.is_empty());
    let prod = &results[0][0];
    assert_eq!(prod.atom_count(), 3, "C + O(new) + N = 3 atoms");
}

#[test]
fn double_bond_product() {
    let mol = parse("CC").unwrap();
    let results = run_reactants("[C:1][C:2]>>[C:1]=[C:2]", &[&mol]).unwrap();
    assert!(!results.is_empty());
    let prod = &results[0][0];
    assert_eq!(prod.atom_count(), 2);
    let bond_orders: Vec<BondOrder> = prod.bonds().map(|(_, b)| b.order).collect();
    assert!(
        bond_orders.contains(&BondOrder::Double),
        "product must contain a double bond"
    );
}

#[test]
fn substituent_carry_through() {
    // Methylamine + acetyl chloride → N-methylacetamide (5 heavy atoms)
    // CH3-NH2 + CH3-C(=O)-Cl → CH3-C(=O)-NH-CH3
    let methylamine = parse("NC").unwrap();
    let acetyl_cl = parse("CC(=O)Cl").unwrap();
    let results = run_reactants(
        "[N:1].[C:2](=O)Cl>>[C:2](=O)[N:1]",
        &[&methylamine, &acetyl_cl],
    )
    .unwrap();
    assert!(!results.is_empty(), "must produce at least one product set");
    let prod = &results[0][0];
    assert_eq!(
        prod.atom_count(),
        5,
        "N-methylacetamide has 5 heavy atoms, got {}",
        prod.atom_count()
    );
}

#[test]
fn bfs_no_leakage_into_other_product_template_atoms() {
    // Issue #13: in diethylamine (CCNCC), the SMIRKS [N:1][C:2]>>[N:1].[C:2]
    // should cleave the N-C bond and produce:
    //   product1 [N:1] = N + right ethyl chain  (3 atoms: N, C, C)
    //   product2 [C:2] = left ethyl fragment     (2 atoms: C, C)
    //
    // Before the #13 fix, the BFS for product2 was seeded from BOTH N and C:2,
    // causing the right ethyl chain (atoms beyond N) to leak into product2
    // → product2 would have 4 atoms instead of 2.
    let diethylamine = parse("CCNCC").unwrap(); // C-C-N-C-C, 5 heavy atoms
    let results = run_reactants("[N:1][C:2]>>[N:1].[C:2]", &[&diethylamine]).unwrap();
    assert!(
        !results.is_empty(),
        "should find at least one N-C bond match"
    );

    // Find a result where product2 ([C:2]) has exactly 2 atoms (ethyl fragment)
    // — this is only possible when BFS does NOT leak the other ethyl chain.
    let clean_cleavage = results.iter().find(|ps| {
        ps.len() == 2
            && ((ps[0].atom_count() == 3 && ps[1].atom_count() == 2)
                || (ps[0].atom_count() == 2 && ps[1].atom_count() == 3))
    });
    assert!(
        clean_cleavage.is_some(),
        "expected at least one product set with sizes {{3, 2}} (N+ethyl, ethyl); \
             all sets: {:?}",
        results
            .iter()
            .map(|ps| ps.iter().map(|p| p.atom_count()).collect::<Vec<_>>())
            .collect::<Vec<_>>()
    );
}

#[test]
fn single_product_no_leakage_from_other_template_core() {
    // Ethane cleavage: each product should be a single carbon atom.
    let ethane = parse("CC").unwrap();
    let results = run_reactants("[C:1][C:2]>>[C:1].[C:2]", &[&ethane]).unwrap();
    assert!(!results.is_empty());
    for ps in &results {
        assert_eq!(ps.len(), 2, "two product templates → two products");
        assert_eq!(ps[0].atom_count(), 1, "each product is a single carbon");
        assert_eq!(ps[1].atom_count(), 1, "each product is a single carbon");
    }
}

// ── Stereo SMIRKS tests ───────────────────────────────────────────────────

#[test]
fn stereo_cleared_when_all_neighbors_are_unmapped_template_literals() {
    // Product template [C:1] has no chirality, and all three substituents
    // (F, Cl, Br) are unmapped literal atoms in BOTH templates -- SMIRKS
    // gives no per-atom reactant<->product correspondence for an unmapped
    // atom (only `:n` atom-map numbers establish identity), so there is no
    // derivable stereo_neighbor_order to carry the inherited flag with.
    // Fail-closed: clear rather than keep a flag with an undefined order
    // (chematic-cip already treats order-less chirality as unassigned, so
    // keeping the raw flag here was never actually meaningful downstream).
    let mol = parse("[C@@H](F)(Cl)Br").unwrap();
    let results = run_reactants("[C@@H:1](F)(Cl)Br>>[C:1](F)(Cl)Br", &[&mol]).unwrap();
    assert!(!results.is_empty(), "should match and produce a product");
    let prod = &results[0][0];
    // The core C atom is first in the builder (index 0).
    let core_chirality = prod.atom(AtomIdx(0)).chirality;
    assert_eq!(
        core_chirality,
        Chirality::None,
        "no derivable stereo_neighbor_order for an all-unmapped-substituent \
             inherited flag must clear chirality, not keep an order-less flag"
    );
}

#[test]
fn stereo_preserved_when_all_neighbors_are_mapped_and_remappable() {
    // Every substituent around the stereocenter is atom-mapped, so each
    // has a real, unique reactant->product correspondence via src_to_new.
    // This is the case the remap mechanism *can* and must handle: keep
    // both the flag and a validated, remapped stereo_neighbor_order.
    let mol = parse("[C@@H](F)(Cl)Br").unwrap();
    let results = run_reactants(
        "[C@@H:1]([F:2])([Cl:3])[Br:4]>>[C:1]([F:2])([Cl:3])[Br:4]",
        &[&mol],
    )
    .unwrap();
    assert!(!results.is_empty(), "should match and produce a product");
    let prod = &results[0][0];
    let core_chirality = prod.atom(AtomIdx(0)).chirality;
    assert_eq!(
        core_chirality,
        Chirality::Clockwise,
        "source @@ chirality must be preserved when every neighbor is \
             mapped and remappable, with an identity-order template"
    );
    assert!(
        prod.stereo_neighbor_order(AtomIdx(0)).is_some(),
        "a kept inherited flag must always carry a validated stereo_neighbor_order"
    );
}

#[test]
fn stereo_unmapped_leaving_groups_fully_removed_clears_chirality() {
    // Blocker #1 from review: the reactant's literal unmapped
    // substituents (F, Cl, Br) are removed entirely by the product
    // template (bare [C:1], zero neighbours) -- bug_b_topology_unchanged
    // must not special-case an empty product-side unmapped set as
    // "unchanged" when the reactant side was non-empty. In the current
    // architecture this is enforced twice over: the topology gate itself
    // (fixed here) AND independently by remap_reactant_stereo_order,
    // since a literal/unmapped reactant-template neighbour never has a
    // src_to_new entry regardless of what the product does with it. The
    // observable behavior this test pins is the one the reviewer asked
    // for either way: chirality must clear, not survive.
    let mol = parse("[C@@H](F)(Cl)Br").unwrap();
    let results = run_reactants("[C@@H:1](F)(Cl)Br>>[C:1]", &[&mol]).unwrap();
    assert!(!results.is_empty(), "should match and produce a product");
    let prod = &results[0][0];
    assert_eq!(
        prod.atom(AtomIdx(0)).chirality,
        Chirality::None,
        "deleting all of the stereocenter's unmapped substituents must clear chirality"
    );
}

#[test]
fn stereo_remote_reaction_preserves_stereocenter() {
    // The bond change (Br leaving, N arriving) happens two bonds away
    // from the stereocenter, through a mapped, unchanged carrier chain
    // (F, Cl, C4 all keep the same atom-map numbers and the same
    // relative template order on both sides) -- the stereocenter's own
    // immediate neighbor set and order are completely untouched by the
    // remote reaction. Configuration must survive.
    let mol = parse("[C@@H](F)(Cl)CCBr").unwrap();
    let amine = parse("N").unwrap();
    let results = run_reactants(
        "[C@@H:1]([F:2])([Cl:3])[C:4][C:5][Br:6].[N:7]\
             >>[C:1]([F:2])([Cl:3])[C:4][C:5][N:7]",
        &[&mol, &amine],
    )
    .unwrap();
    assert!(!results.is_empty(), "should match and produce a product");
    let prod = &results[0][0];
    assert_eq!(
        prod.atom(AtomIdx(0)).chirality,
        Chirality::Clockwise,
        "a remote bond change 2 bonds away must not disturb the \
             stereocenter's own unchanged, fully-mapped neighbor set"
    );
    assert!(
        prod.stereo_neighbor_order(AtomIdx(0)).is_some(),
        "preserved chirality must carry a validated stereo_neighbor_order"
    );
}

#[test]
fn stereo_survives_16_plus_atom_map_relabelings_smiles_and_cip_invariant() {
    // The chemistry must not depend on which integers a SMIRKS author
    // picks for `:n` map numbers -- only the correspondence they encode.
    // Re-run the same reaction with 20 distinct map-number assignments
    // (same template text shape and order each time, only the four
    // integers change) and assert canonical SMILES and accurate CIP
    // agree with the first run every time.
    let mol = parse("[C@@H](F)(Cl)Br").unwrap();
    // Central atom is fixed at map `:1`; offset the other three ranges
    // so none of them ever collides with it or each other.
    let labelings: Vec<[u32; 3]> = (2..=21).map(|i| [i, i + 100, i + 200]).collect();
    assert!(labelings.len() >= 16, "need at least 16 relabelings");

    let mut baseline_smiles: Option<String> = None;
    let mut baseline_cip: Option<Option<chematic_core::CipCode>> = None;
    for [b, c, d] in labelings {
        let smirks =
            format!("[C@@H:1]([F:{b}])([Cl:{c}])[Br:{d}]>>[C:1]([F:{b}])([Cl:{c}])[Br:{d}]");
        let results = run_reactants(&smirks, &[&mol]).unwrap();
        assert!(!results.is_empty(), "labeling {b},{c},{d} must still match");
        let prod = &results[0][0];
        let smi = canonical(prod);
        let cip = chematic_chem::assign_cip_with_mode(prod, chematic_chem::CipMode::Accurate)
            .unwrap()
            .get(AtomIdx(0));

        match (&baseline_smiles, &baseline_cip) {
            (None, None) => {
                baseline_smiles = Some(smi);
                baseline_cip = Some(cip);
            }
            (Some(base_smi), Some(base_cip)) => {
                assert_eq!(
                    &smi, base_smi,
                    "canonical SMILES must be invariant to atom-map relabeling \
                         (labeling {b},{c},{d})"
                );
                assert_eq!(
                    &cip, base_cip,
                    "accurate CIP code must be invariant to atom-map relabeling \
                         (labeling {b},{c},{d})"
                );
            }
            _ => unreachable!(),
        }
    }
}

#[test]
fn stereo_chirality_implies_valid_stereo_neighbor_order_invariant() {
    // Cross-cutting invariant the reviewer required: whenever an atom
    // carries a non-None chirality flag, it must always have a valid
    // stereo_neighbor_order alongside it -- never a raw flag with no
    // order (which chematic-cip and the SMILES writer cannot interpret
    // meaningfully). Checked across every atom of a representative
    // spread of already-covered product-generating scenarios.
    let assert_invariant_holds = |mol: &Molecule, label: &str| {
        for i in 0..mol.atom_count() {
            let idx = AtomIdx(i as u32);
            if mol.atom(idx).chirality != Chirality::None {
                assert!(
                    mol.stereo_neighbor_order(idx).is_some(),
                    "{label}: atom {i} has chirality {:?} but no stereo_neighbor_order",
                    mol.atom(idx).chirality
                );
            }
        }
    };

    let mol = parse("[C@@H](F)(Cl)Br").unwrap();
    let chain_mol = parse("[C@@H](F)(Cl)CCBr").unwrap();
    let amine = parse("N").unwrap();

    let scenarios: Vec<(&str, Vec<Molecule>)> = vec![
        (
            "identity template order",
            vec![
                run_reactants("[C@@H:1](F)(Cl)Br>>[C@@H:1](F)(Cl)Br", &[&mol]).unwrap()[0][0]
                    .clone(),
            ],
        ),
        (
            "reordered template inverts",
            vec![
                run_reactants("[C@@H:1](F)(Cl)Br>>[C@@H:1](Cl)(F)Br", &[&mol]).unwrap()[0][0]
                    .clone(),
            ],
        ),
        (
            "all neighbors mapped and remappable",
            vec![
                run_reactants(
                    "[C@@H:1]([F:2])([Cl:3])[Br:4]>>[C:1]([F:2])([Cl:3])[Br:4]",
                    &[&mol],
                )
                .unwrap()[0][0]
                    .clone(),
            ],
        ),
        (
            "all unmapped substituents deleted",
            vec![run_reactants("[C@@H:1](F)(Cl)Br>>[C:1]", &[&mol]).unwrap()[0][0].clone()],
        ),
        (
            "substituent replacement",
            vec![run_reactants("[C@@H:1](F)(Cl)Br>>[C:1](F)(Cl)I", &[&mol]).unwrap()[0][0].clone()],
        ),
        (
            "remote reaction",
            vec![
                run_reactants(
                    "[C@@H:1]([F:2])([Cl:3])[C:4][C:5][Br:6].[N:7]\
                         >>[C:1]([F:2])([Cl:3])[C:4][C:5][N:7]",
                    &[&chain_mol, &amine],
                )
                .unwrap()[0][0]
                    .clone(),
            ],
        ),
    ];

    for (label, products) in &scenarios {
        for prod in products {
            assert_invariant_holds(prod, label);
        }
    }
}

#[test]
fn stereo_inverted_by_template() {
    // Product template [C@H:1] has @ (CounterClockwise) → overrides source @@ (Clockwise).
    let mol = parse("[C@@H](F)(Cl)Br").unwrap();
    let results = run_reactants("[C@@H:1](F)(Cl)Br>>[C@H:1](F)(Cl)Br", &[&mol]).unwrap();
    assert!(!results.is_empty(), "should match and produce a product");
    let prod = &results[0][0];
    let core_chirality = prod.atom(AtomIdx(0)).chirality;
    assert_eq!(
        core_chirality,
        Chirality::CounterClockwise,
        "product template @ must override source @@ → CounterClockwise"
    );
}

#[test]
fn stereo_identity_template_order_preserves_configuration() {
    // Product template repeats the exact same neighbour order as the
    // reactant template (F, Cl, Br) with an explicit @@ -- the simplest
    // "nothing changed" case for the parity-aware correction.
    let mol = parse("[C@@H](F)(Cl)Br").unwrap();
    let results = run_reactants("[C@@H:1](F)(Cl)Br>>[C@@H:1](F)(Cl)Br", &[&mol]).unwrap();
    let prod = &results[0][0];
    let smi = canonical(prod);
    let input_smi = canonical(&mol);
    assert_eq!(
        smi, input_smi,
        "identity template order should reproduce the input unchanged"
    );
}

#[test]
fn stereo_reordered_template_inverts_absolute_configuration() {
    // Issue found while surveying RDKit's open issues (analogous to
    // RDKit #9257): reordering two substituents in the product template
    // while keeping the SAME @@ symbol must invert the absolute
    // configuration -- the symbol describes an order-relative sense,
    // not an absolute one. Cross-checked against a live RDKit oracle
    // (`rdkit.Chem.rdChemReactions`): RDKit's `reordered` output equals
    // its own `inverted` output, both differing from `identity`.
    let mol = parse("[C@@H](F)(Cl)Br").unwrap();
    let identity =
        canonical(&run_reactants("[C@@H:1](F)(Cl)Br>>[C@@H:1](F)(Cl)Br", &[&mol]).unwrap()[0][0]);
    let reordered =
        canonical(&run_reactants("[C@@H:1](F)(Cl)Br>>[C@@H:1](Cl)(F)Br", &[&mol]).unwrap()[0][0]);
    let explicit_invert =
        canonical(&run_reactants("[C@@H:1](F)(Cl)Br>>[C@H:1](F)(Cl)Br", &[&mol]).unwrap()[0][0]);
    assert_ne!(
        reordered, identity,
        "reordering two substituents under the same @@ symbol must change \
             the absolute configuration"
    );
    assert_eq!(
        reordered, explicit_invert,
        "reordering two substituents under @@ must match the explicit-@ \
             (same order, opposite symbol) result -- both express the same \
             inverted configuration"
    );
}

#[test]
fn stereo_substituent_replacement_clears_chirality() {
    // Neighbour SET changes (Br -> I): the old flag can no longer mean
    // anything -- must clear to Chirality::None, not carry a stale @@
    // onto a differently-substituted center. Pre-existing behavior,
    // now driven by correct_product_stereo's Bug-B gate instead of the
    // deleted ad hoc heuristic; must not regress.
    let mol = parse("[C@@H](F)(Cl)Br").unwrap();
    let results = run_reactants("[C@@H:1](F)(Cl)Br>>[C:1](F)(Cl)I", &[&mol]).unwrap();
    let prod = &results[0][0];
    let core_chirality = prod.atom(AtomIdx(0)).chirality;
    assert_eq!(
        core_chirality,
        Chirality::None,
        "substituent replacement (Br->I) must clear chirality, not preserve a stale flag"
    );
}

#[test]
fn stereo_lost_mapped_neighbor_clears_spurious_chirality() {
    // Issue found while surveying RDKit's open issues: a mapped/core
    // neighbour (N) is dropped by the product template, while the
    // template's own C:1 atom carries no explicit chirality. The old
    // invalidation heuristic only ever compared *unmapped* substituent
    // element sets (both empty here, since every neighbour of C:1 is
    // mapped) and missed this case entirely, producing a spurious
    // chirality tag on a carbon with 3 identical implicit hydrogens.
    let mol = parse("[C@H](N)(O)Cl").unwrap();
    let results = run_reactants("[C@H:1]([N:2])([O:3])Cl>>[C:1][O:3]", &[&mol]).unwrap();
    let prod = &results[0][0];
    let smi = canonical(prod);
    assert!(
        !smi.contains('@'),
        "carbon losing its N neighbour is no longer a stereocenter (3 identical \
             implicit H's) -- product must carry no chirality symbol, got {smi}"
    );
}

#[test]
fn stereo_polyol_ring_contraction_no_spurious_ch2oh_chirality() {
    // The originally reported repro (analogous to RDKit #9257), minimized
    // versions of which are the two tests directly above. A pyranose ->
    // furanose ring contraction where one ring carbon (mapped :2) loses
    // its own explicit chirality in the product template and becomes a
    // plain CH2OH (degree 2, bonded to a degree-1 O) -- pre-fix this
    // atom kept a spurious inherited chirality tag ([C@H2], a carbon
    // with 3 identical implicit hydrogens). Does NOT assert a specific
    // sign at the two atoms whose absolute configuration this fix does
    // not by itself resolve either way (mapped :11, unchanged reactant->
    // product template and therefore architecturally untouched by this
    // fix; and :15, the debated center in the original RDKit report) --
    // only that the reaction produces exactly 4 real stereocenters, not
    // 5, and that the correct atom (the one that structurally lost its
    // stereocenter) is the one that lost its tag.
    let smirks = "[O:1][C@H:2]1[O:3][C@H:4]([C:5][C:6])[C@@H:11]([O:12])[C@H:13]([O:14])\
                       [C@H:15]1[O:16]>>[O:1][C:2][C@:15]([O:16])1[O:3][C@H:4]([C:5][C:6])\
                       [C@@H:11]([O:12])[C@H:13]1[O:14]";
    let mol = parse("CC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O").unwrap();
    let results = run_reactants(smirks, &[&mol]).unwrap();
    let prod = &results[0][0];

    let n = prod.atom_count();
    let chiral_atoms: Vec<AtomIdx> = (0..n)
        .map(|i| AtomIdx(i as u32))
        .filter(|&a| prod.atom(a).chirality != Chirality::None)
        .collect();
    assert_eq!(
        chiral_atoms.len(),
        4,
        "expected exactly 4 real stereocenters in the product, got {}: {:?}",
        chiral_atoms.len(),
        chiral_atoms
    );

    // Structurally identify the CH2OH carbon: degree 2, bonded to a
    // degree-1 oxygen -- not by a hardcoded index, since atom order
    // depends on build_product's own construction order.
    let ch2oh = (0..n).map(|i| AtomIdx(i as u32)).find(|&a| {
        let atom = prod.atom(a);
        atom.element == chematic_core::Element::C
            && prod.degree(a) == 2
            && prod.neighbors(a).any(|(nb, _)| {
                prod.atom(nb).element == chematic_core::Element::O && prod.degree(nb) == 1
            })
    });
    let ch2oh = ch2oh.expect("product must contain a CH2OH carbon");
    assert_eq!(
        prod.atom(ch2oh).chirality,
        Chirality::None,
        "the CH2OH carbon (lost its own chirality in the product template, degree 2, \
             no longer has 4 distinct substituents) must not carry a stereo tag"
    );
}

// ── run_reactants_strict tests ────────────────────────────────────────────

#[test]
fn strict_mode_excludes_substituents() {
    // Methylamine (NC): in normal mode [N:1]>>[N:1] carries C through as substituent.
    // In strict mode only N is returned (no C).
    let mol = parse("NC").unwrap();
    let normal = run_reactants("[N:1]>>[N:1]", &[&mol]).unwrap();
    let strict = run_reactants_strict("[N:1]>>[N:1]", &[&mol]).unwrap();
    assert!(!normal.is_empty());
    assert!(!strict.is_empty());
    let normal_atoms = normal[0][0].atom_count();
    let strict_atoms = strict[0][0].atom_count();
    assert!(
        normal_atoms > strict_atoms,
        "normal mode carries substituent C (got {normal_atoms}), \
             strict mode only mapped N (got {strict_atoms})"
    );
    assert_eq!(strict_atoms, 1, "strict mode: only the mapped N atom");
}

#[test]
fn strict_mode_bond_cleavage() {
    // Ethane cleavage: strict mode gives 1-atom products, same as normal here
    // (no unmapped substituents on either C).
    let ethane = parse("CC").unwrap();
    let results = run_reactants_strict("[C:1][C:2]>>[C:1].[C:2]", &[&ethane]).unwrap();
    assert!(!results.is_empty());
    for ps in &results {
        assert_eq!(ps[0].atom_count(), 1);
        assert_eq!(ps[1].atom_count(), 1);
    }
}

// ── Issue #18: product bracket notation cleanup ───────────────────────────

#[test]
fn product_removes_bracket_from_bare_bracket_atoms() {
    // Issue #18: [O:1] in product template (hydrogen_count=Some(0)) must produce
    // clean `O` SMILES, not `[O]`.
    use chematic_smiles::canonical_smiles;
    let mol = parse("OCC").unwrap();
    let results = run_reactants("[OH:1]>>[O:1]", &[&mol]).unwrap();
    assert!(!results.is_empty(), "should match hydroxyl");
    let prod_smi = canonical_smiles(&results[0][0]);
    assert!(
        !prod_smi.contains("[O]"),
        "bare [O:1] product must write as O, not [O], got: {prod_smi}"
    );
}

#[test]
fn product_preserves_explicit_h_from_template() {
    // [NH2:1] in product template specifies 2H explicitly — the built
    // atom's `hydrogen_count` must be `Some(2)`, not silently dropped or
    // recomputed.
    //
    // Checks the data-level field directly rather than the canonical
    // SMILES string's bracket notation (issue #205): after
    // `initial_invariant`/`emit_atom`'s explicit/implicit-H-count
    // unification fix, `canonical_smiles` correctly stops forcing
    // brackets for an atom whose explicit H count merely repeats what
    // organic-subset valence inference would already give -- this atom
    // (N bonded to one carbon) infers to 2 implicit H anyway, so its
    // canonical form is now the fully standard, bracket-free "CN"
    // (methylamine), not "C[NH2]". The template's explicit
    // specification is still honored -- it is what the "2" came from --
    // just no longer forced into visible bracket notation once it's
    // redundant with inference.
    let mol = parse("NC").unwrap();
    let results = run_reactants("[N:1]>>[NH2:1]", &[&mol]).unwrap();
    assert!(!results.is_empty(), "should match amine N");
    let product = &results[0][0];
    let n_atom = product
        .atoms()
        .find(|(_, a)| a.element == chematic_core::Element::N)
        .map(|(_, a)| a)
        .expect("product must contain the templated N atom");
    assert_eq!(
        n_atom.hydrogen_count,
        Some(2),
        "explicit [NH2:1] in product must set hydrogen_count = Some(2)"
    );
}

#[test]
fn reaction_derived_matches_direct_parse_chlorobenzene() {
    // The exact fixture from kent-tokyo/renkin PR #65: a reaction-
    // derived molecule (built via `run_reactants`, whose product-
    // template atom comes from bracket notation `[Cl]` in the SMIRKS)
    // must canonicalize identically to a directly-parsed organic-subset
    // "Clc1ccccc1" of the same compound (issue #205).
    use chematic_smiles::canonical_smiles;
    let fwd = "[c:1][Br]>>[c:1][Cl]";
    let known = parse("Brc1ccccc1").unwrap();
    let results = run_reactants(fwd, &[&known]).unwrap();
    assert_eq!(results.len(), 1);
    assert_eq!(results[0].len(), 1);
    let reaction_derived = &results[0][0];
    let direct = parse("Clc1ccccc1").unwrap();

    // Structural sanity first: same atom count and element multiset,
    // so canonical-string equality below is a real invariance proof,
    // not an assumption that the two are the same molecule.
    assert_eq!(reaction_derived.atom_count(), direct.atom_count());
    let mut els_a: Vec<_> = reaction_derived.atoms().map(|(_, a)| a.element).collect();
    let mut els_b: Vec<_> = direct.atoms().map(|(_, a)| a.element).collect();
    els_a.sort();
    els_b.sort();
    assert_eq!(els_a, els_b, "same element multiset required");

    assert_eq!(
        canonical_smiles(reaction_derived),
        canonical_smiles(&direct),
        "reaction-derived and directly-parsed chlorobenzene must canonicalize identically"
    );
}

// ── Issue #20: SMIRKS stereo filtering ───────────────────────────────────

#[test]
fn stereo_filter_rejects_wrong_enantiomer() {
    // Issue #20: SMIRKS with @@ reactant template must match @@ but NOT @ reactant.
    let l_ala = parse("N[C@@H](C)C(=O)O").unwrap(); // L-alanine (@@)
    let d_ala = parse("N[C@H](C)C(=O)O").unwrap(); // D-alanine (@)

    let smirks = "[N:1][C@@H:2](C)C(=O)O>>[N:1][C@@H:2](C)C(=O)O";
    let results_l = run_reactants(smirks, &[&l_ala]).unwrap();
    let results_d = run_reactants(smirks, &[&d_ala]).unwrap();

    assert!(
        !results_l.is_empty(),
        "L-alanine (@@) must match @@ template"
    );
    assert!(
        results_d.is_empty(),
        "D-alanine (@) must NOT match @@ template (stereo filter, issue #20)"
    );
}

#[test]
fn incomplete_tetrahedral_template_cannot_verify_handedness() {
    // The reactant's Br substituent is absent from the chiral template.
    // Its local neighbour orders cannot be aligned, so accepting either
    // enantiomer would silently discard the template's stereo constraint.
    let template = parse("[C@H](F)Cl").unwrap();
    assert!(template.stereo_neighbor_order(AtomIdx(0)).is_some());
    for input in ["Br[C@H](F)Cl", "Br[C@@H](F)Cl"] {
        let reactant = parse(input).unwrap();
        let mut mapping = FxHashMap::default();
        mapping.insert(0, AtomIdx(1));
        mapping.insert(1, AtomIdx(2));
        mapping.insert(2, AtomIdx(3));
        assert!(reactant.stereo_neighbor_order(AtomIdx(1)).is_some());
        let products = run_reactants("[C@H:1](F)Cl>>[C@H:1](F)Cl", &[&reactant]).unwrap();
        assert!(
            products.is_empty(),
            "must not invent handedness for {input}"
        );
        assert!(!smirks_chirality_ok(&template, &reactant, &mapping));
        assert!(
            !rdkit_profile("[C@H:1](F)Cl>>[C@H:1](F)Cl", &[input])
                .unwrap()
                .is_empty(),
            "the RDKit compatibility profile stays permissive for {input}"
        );
    }
}

#[test]
fn stereo_neutral_smirks_matches_both_enantiomers() {
    // SMIRKS without @/@@ must still match both enantiomers (backward compat).
    let l_ala = parse("N[C@@H](C)C(=O)O").unwrap();
    let d_ala = parse("N[C@H](C)C(=O)O").unwrap();
    let smirks = "[N:1][CH:2](C)C(=O)O>>[N:1][CH:2](C)C(=O)O";
    let r_l = run_reactants(smirks, &[&l_ala]).unwrap();
    let r_d = run_reactants(smirks, &[&d_ala]).unwrap();
    assert!(!r_l.is_empty(), "L-alanine must match non-stereo template");
    assert!(!r_d.is_empty(), "D-alanine must match non-stereo template");
}

/// Products under the RDKit profile, as canonical SMILES per set.
fn rdkit_profile(smirks: &str, reactants: &[&str]) -> Result<Vec<String>, &'static str> {
    let mols: Vec<Molecule> = reactants.iter().map(|s| parse(s).unwrap()).collect();
    let refs: Vec<&Molecule> = mols.iter().collect();
    match run_reactants_traced_rdkit_2026_03_6(smirks, &refs, &Default::default()).unwrap() {
        RdkitProfileOutcome::Report(report) => Ok(report
            .products
            .iter()
            .map(|set| {
                set.iter()
                    .map(|p| canonical(&p.molecule))
                    .collect::<Vec<_>>()
                    .join(".")
            })
            .collect()),
        RdkitProfileOutcome::Unsupported(reason) => Err(reason.reason_code()),
    }
}

fn same_molecule(a: &str, b: &str) -> bool {
    canonical(&parse(a).unwrap()) == canonical(&parse(b).unwrap())
}

#[test]
fn rdkit_profile_reproduces_rdkit_alanine_identity_rows(/* issue #734 */) {
    // RDKit 2026.03.6: the reactant `@@` does not filter, and with all
    // other neighbours unmapped the product takes the reactant's raw tag
    // in the product's bond order — so the reordered L-alanine spelling
    // comes out as D-alanine, exactly as RDKit gives it.
    let smirks = "[N:1][C@@H:2](C)C(=O)O>>[N:1][C@@H:2](C)C(=O)O";
    for (reactant, rdkit) in [
        ("N[C@@H](C)C(=O)O", "C[C@H](N)C(=O)O"),
        ("C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O"),
        ("N[C@H](C)C(=O)O", "C[C@@H](N)C(=O)O"),
    ] {
        let got = rdkit_profile(smirks, &[reactant]).unwrap();
        assert_eq!(got.len(), 1, "{reactant}");
        assert!(same_molecule(&got[0], rdkit), "{reactant}: {got:?}");
    }
}

#[test]
fn rdkit_profile_inversion_flags_follow_rdkit() {
    // Invert (tags differ), with every neighbour traced: geometric.
    let got = rdkit_profile(
        "[C:1][C@H:2]([N:3])[C:4]>>[C:1][C@@H:2]([N:3])[C:4]",
        &["N[C@@H](C)CC"],
    )
    .unwrap();
    assert!(
        got.iter().all(|p| same_molecule(p, "CC[C@@H](C)N")),
        "{got:?}"
    );
    // Replacing a leaving group keeps the traced neighbours' frame.
    let got = rdkit_profile(
        "[C:1][C@H:2](Br)[C:3]>>[C:1][C@@H:2](O)[C:3]",
        &["C[C@@H](Br)CC"],
    )
    .unwrap();
    assert!(
        got.iter().all(|p| same_molecule(p, "CC[C@H](C)O")),
        "{got:?}"
    );
    // Remove: stereo in the reactant template only.
    let got = rdkit_profile(
        "[C:1][C@H:2]([O:3])[C:4]>>[C:1][CH:2]([O:3])[C:4]",
        &["C[C@@H](O)CC"],
    )
    .unwrap();
    assert!(got.iter().all(|p| same_molecule(p, "CCC(C)O")), "{got:?}");
    // Ring cis/trans: the mapped centre inverts, the carried one stays.
    let got = rdkit_profile(
        "[C:1][C@H:2]([N:3])[C:4]>>[C:1][C@@H:2]([N:3])[C:4]",
        &["N[C@H]1CC[C@H](CC1)O"],
    )
    .unwrap();
    assert!(
        got.iter().all(|p| same_molecule(p, "N[C@H]1CC[C@H](O)CC1")),
        "{got:?}"
    );
}

#[test]
fn rdkit_profile_declines_unrecorded_ring_closure_order() {
    // Raw tag copy on a ring-fusion atom with two ring closures: RDKit's
    // answer depends on the ring-closure numbers, which a parsed molecule
    // does not keep.
    assert_eq!(
        rdkit_profile(
            "[CH3:1][C@:2]([CH2:3])([CH2:4])[CH:5]>>[CH3:1][C@:2](O)(N)F",
            &["C[C@]12CCC[C@H]1CCC2"],
        ),
        Err("ambiguous_stereo_bond_order")
    );
    // When the atom is no stereocentre in the product, any order does.
    let got = rdkit_profile(
        "[CH3:1][C@@:2](C)(C)C>>[CH3:1][C@@:2](C)(C)C",
        &["C[C@]12CCC[C@H]1CCC2"],
    )
    .unwrap();
    assert_eq!(got.len(), 6, "{got:?}");
}

#[test]
fn rdkit_profile_reads_smarts_only_stereo_and_ez_matching() {
    // SMARTS-only reactant templates may carry `@` here (RDKit matches
    // without it); the native constructor still refuses them.
    let smirks = "[N;H2:1][C@@H:2]([CH3:3])[C:4]>>[N:1][C@@H:2]([CH3:3])[C:4]";
    assert!(PreparedReaction::new(smirks).is_err());
    let got = rdkit_profile(smirks, &["C[C@H](N)C(=O)O"]).unwrap();
    assert!(same_molecule(&got[0], "C[C@H](N)C(=O)O"), "{got:?}");
    // RDKit ignores reactant `/` `\` when matching; the native check
    // rejecting a match is reported rather than turned into no products.
    assert_eq!(
        rdkit_profile(
            "[C:1]/[C:2]=[C:3]/[C:4]>>[C:1][C:2][C:3][C:4]",
            &["C/C=C\\C"]
        ),
        Err("ez_reactant_template_semantics")
    );
    assert!(
        rdkit_profile(
            "[C:1]/[C:2]=[C:3]/[C:4]>>[C:1][C:2][C:3][C:4]",
            &["C/C=C/C"]
        )
        .is_ok()
    );
}

#[test]
fn rdkit_profile_drops_tags_rdkit_parsing_drops() {
    // `C[C@H](C)O` has no stereocentre; RDKit's parser drops the tag and
    // its bracket H, so the new centre comes out unspecified.
    let got = rdkit_profile("[CH3:1][C:2]>>[CH2:1](F)[C:2]", &["C[C@H](C)O"]).unwrap();
    assert!(got.iter().all(|p| same_molecule(p, "CC(O)CF")), "{got:?}");
    // A real centre keeps its tag.
    let got = rdkit_profile("[CH3:1][C:2]>>[CH2:1](F)[C:2]", &["CC[C@H](C)O"]).unwrap();
    assert!(
        got.iter().any(|p| same_molecule(p, "C[C@H](O)CCF")),
        "{got:?}"
    );
    // So does a 1,4-disubstituted ring pair (cis/trans).
    let got = rdkit_profile("[C:1]>>[C:1]", &["C[C@H]1CC[C@@H](C)CC1"]).unwrap();
    assert!(
        got.iter()
            .all(|p| same_molecule(p, "C[C@H]1CC[C@@H](C)CC1")),
        "{got:?}"
    );
}

#[test]
fn carried_stereocentres_keep_their_configuration() {
    // A stereocentre carried through unchanged keeps its configuration
    // in whatever bond order the product has (v1.0.33 inverted it when
    // the template atoms came first in the product).
    for reactant in ["N[C@@H](C)C(O)=O", "C[C@H](N)C(O)=O", "OC(=O)[C@@H](N)C"] {
        let mol = parse(reactant).unwrap();
        for set in run_reactants("[C:1](=[O:2])[OH:3]>>[C:1](=[O:2])OC", &[&mol]).unwrap() {
            assert!(
                same_molecule(&canonical(&set[0]), "COC(=O)[C@H](C)N"),
                "{reactant}: {}",
                canonical(&set[0])
            );
        }
    }
}

#[test]
fn bond_between_mapped_atoms_absent_from_template_is_kept() {
    // RDKit (github #1387): the reactant template does not bond the two
    // mapped ring atoms, so their ring bond is carried, not broken.
    let mol = parse("NC1CC1").unwrap();
    for set in run_reactants("[N:1][C:2][C:3]>>[N:1][C:2][C:3]", &[&mol]).unwrap() {
        assert!(same_molecule(&canonical(&set[0]), "NC1CC1"));
    }
    // A bond the reactant template spells and the product omits is broken.
    for set in run_reactants("[C:2]1[C:3][C:4]1>>[C:2][C:3].[C:4]", &[&mol]).unwrap() {
        assert_eq!(set.len(), 2);
    }
}

#[test]
fn stereo_filter_same_config_different_write_order() {
    // Parity-aware matching must accept the same absolute configuration
    // regardless of SMILES atom write order (the confirmed bug in raw flag
    // comparison). Both molecules ARE L-alanine.
    //   Form A: N[C@@H](C)C(=O)O  — N written first, stored as Clockwise
    //   Form B: C[C@H](N)C(=O)O   — C_methyl first, stored as CounterClockwise
    // A raw flag comparison would reject Form B against an @@ template.
    let l_form_a = parse("N[C@@H](C)C(=O)O").unwrap();
    let l_form_b = parse("C[C@H](N)C(=O)O").unwrap(); // same absolute config, diff write order
    let d_form = parse("N[C@H](C)C(=O)O").unwrap(); // D-alanine (opposite config)

    let smirks = "[N:1][C@@H:2](C)C(=O)O>>[N:1][C@@H:2](C)C(=O)O";

    let r_a = run_reactants(smirks, &[&l_form_a]).unwrap();
    let r_b = run_reactants(smirks, &[&l_form_b]).unwrap();
    let r_d = run_reactants(smirks, &[&d_form]).unwrap();

    assert!(!r_a.is_empty(), "L-alanine form A (N-first @@) must match");
    assert!(
        !r_b.is_empty(),
        "L-alanine form B (C-first @, same absolute config) must also match \
             — parity-aware comparison required"
    );
    assert!(r_d.is_empty(), "D-alanine must still be rejected");
}

// ── RDKit #9339: orphaned stereo bonds cleared from products ─────────────

#[test]
fn smirks_reaction_clears_orphaned_stereo_bonds() {
    // (E)-2-butene C/C=C/C has Up/Down bonds adjacent to the C=C double bond.
    // After SMIRKS [C:1]=[C:2]>>[C:1][C:2] the double bond is reduced to a
    // single bond.  The Up/Down single bonds on C0-C1 and C2-C3 are no longer
    // adjacent to ANY double bond → orphaned → must be cleared (RDKit PR #9339).
    let mol = parse("C/C=C/C").unwrap(); // (E)-2-butene
    let results = run_reactants("[C:1]=[C:2]>>[C:1][C:2]", &[&mol]).unwrap();
    assert!(!results.is_empty(), "should produce at least one product");
    for prod_set in &results {
        for prod in prod_set {
            for (_, bond) in prod.bonds() {
                assert_ne!(
                    bond.order,
                    BondOrder::Up,
                    "stray Up bond in product after C=C→C-C (RDKit #9339)"
                );
                assert_ne!(
                    bond.order,
                    BondOrder::Down,
                    "stray Down bond in product after C=C→C-C (RDKit #9339)"
                );
            }
        }
    }
}

#[test]
fn smirks_preserves_stereo_bonds_adjacent_to_remaining_double() {
    // If the double bond is kept unchanged, the *exact* E/Z geometry must be
    // preserved — not merely "some Up/Down bond survives" (which passed even
    // while the geometry was flipping E↔Z, issue #50). Verify by comparing the
    // product's canonical SMILES to the canonical SMILES of the known input.
    use chematic_smiles::canonical_smiles;
    for input in ["C/C=C/C", "C/C=C\\C"] {
        let mol = parse(input).unwrap();
        let results = run_reactants("[C:1]=[C:2]>>[C:1]=[C:2]", &[&mol]).unwrap();
        assert!(!results.is_empty());
        let expected = canonical_smiles(&mol);
        let got = canonical_smiles(&results[0][0]);
        assert_eq!(
            got, expected,
            "identity SMIRKS must preserve exact E/Z geometry for {input}"
        );
    }
}

// -----------------------------------------------------------------------
// E/Z double-bond stereo filtering (issue #21)
// -----------------------------------------------------------------------

#[test]
fn ez_stereo_e_template_matches_e_alkene() {
    // Template specifies E: [C:1]/[C:2]=[C:3]/[C:4]
    // E-2-butene (C/C=C/C) should produce 1 result.
    let e_alkene = parse("C/C=C/C").unwrap();
    let smirks = "[C:1]/[C:2]=[C:3]/[C:4]>>[C:1][C:2][C:3][C:4]";
    let results = run_reactants(smirks, &[&e_alkene]).unwrap();
    assert!(!results.is_empty(), "E-template must match E-alkene");
}

#[test]
fn ez_stereo_e_template_rejects_z_alkene() {
    // Template specifies E: [C:1]/[C:2]=[C:3]/[C:4]
    // Z-2-butene (C/C=C\C) should produce 0 results.
    let z_alkene = parse("C/C=C\\C").unwrap();
    let smirks = "[C:1]/[C:2]=[C:3]/[C:4]>>[C:1][C:2][C:3][C:4]";
    let results = run_reactants(smirks, &[&z_alkene]).unwrap();
    assert!(results.is_empty(), "E-template must reject Z-alkene");
}

#[test]
fn ez_stereo_e_template_rejects_unspecified_alkene() {
    let smirks = "[F:1]/[C:2]=[C:3]/[Cl:4]>>[F:1]/[C:2]=[C:3]/[Cl:4]";
    for input in ["FC=CCl", "F/C=CCl"] {
        let unspecified = parse(input).unwrap();
        assert!(
            run_reactants(smirks, &[&unspecified]).unwrap().is_empty(),
            "an E/Z-specific template must not invent geometry from {input}"
        );
    }
    let unspecified = parse("FC=CCl").unwrap();
    let specified = parse("F/C=C/Cl").unwrap();
    assert!(!run_reactants(smirks, &[&specified]).unwrap().is_empty());
    let RdkitProfileOutcome::Report(compatibility) = run_reactants_traced_rdkit_2026_03_6(
        smirks,
        &[&unspecified],
        &ReactionTransformLimits::default(),
    )
    .unwrap() else {
        panic!("the RDKit profile must retain its permissive matching contract");
    };
    assert!(!compatibility.products.is_empty());
}

#[test]
fn ez_stereo_neutral_template_matches_both_geometries() {
    // Template without stereo: [C:1][C:2]=[C:3][C:4]>>[C:1]
    // Both E and Z alkenes should match.
    let e_alkene = parse("C/C=C/C").unwrap();
    let z_alkene = parse("C/C=C\\C").unwrap();
    let smirks = "[C:1][C:2]=[C:3][C:4]>>[C:1]";
    assert!(
        !run_reactants(smirks, &[&e_alkene]).unwrap().is_empty(),
        "neutral template must match E-alkene"
    );
    assert!(
        !run_reactants(smirks, &[&z_alkene]).unwrap().is_empty(),
        "neutral template must match Z-alkene"
    );
}

#[test]
fn ez_stereo_one_sided_template_matches_both_geometries() {
    // Single-sided stereo bond in template: [C:1]/[C:2]=[C:3][C:4]
    // Without both sides specified, E/Z is ambiguous → no filtering.
    let e_alkene = parse("C/C=C/C").unwrap();
    let z_alkene = parse("C/C=C\\C").unwrap();
    let smirks = "[C:1]/[C:2]=[C:3][C:4]>>[C:1]";
    assert!(
        !run_reactants(smirks, &[&e_alkene]).unwrap().is_empty(),
        "one-sided template must match E-alkene"
    );
    assert!(
        !run_reactants(smirks, &[&z_alkene]).unwrap().is_empty(),
        "one-sided template must match Z-alkene"
    );
    let unspecified = parse("CC=CC").unwrap();
    assert!(
        !run_reactants(smirks, &[&unspecified]).unwrap().is_empty(),
        "one-sided template has no E/Z constraint"
    );
}

#[test]
fn ez_stereo_retro_wittig_z_matches_z_hexene() {
    // Retro-Wittig (Z-alkene → two carbonyls).
    // SMIRKS: [C:1]/[C:2]=[C:3]\[C:4]>>[C:1][C:2]=O.[O:3]=[C:4]
    //   reads: C:2 and C:3 on OPPOSITE sides (E/trans for those two)
    //   but the substituents C:1 and C:4 are on the SAME side (Z-selectivity)
    //
    // Z-3-hexene (CC/C=C\CC) should match; E-3-hexene (CC/C=C/CC) should not.
    let z_hexene = parse("CC/C=C\\CC").unwrap();
    let e_hexene = parse("CC/C=C/CC").unwrap();
    let smirks = "[C:1]/[C:2]=[C:3]\\[C:4]>>[C:1][C:2]=O.[O:3]=[C:4]";
    assert!(
        !run_reactants(smirks, &[&z_hexene]).unwrap().is_empty(),
        "Z-template must match Z-3-hexene"
    );
    assert!(
        run_reactants(smirks, &[&e_hexene]).unwrap().is_empty(),
        "Z-template must reject E-3-hexene"
    );
}

#[test]
fn ez_stereo_z_template_matches_z_alkene() {
    // Template specifies Z: [C:1]/[C:2]=[C:3]\[C:4]
    // Z-2-butene (C/C=C\C) should match.
    let z_alkene = parse("C/C=C\\C").unwrap();
    let e_alkene = parse("C/C=C/C").unwrap();
    let smirks = "[C:1]/[C:2]=[C:3]\\[C:4]>>[C:1][C:2][C:3][C:4]";
    assert!(
        !run_reactants(smirks, &[&z_alkene]).unwrap().is_empty(),
        "Z-template must match Z-alkene"
    );
    assert!(
        run_reactants(smirks, &[&e_alkene]).unwrap().is_empty(),
        "Z-template must reject E-alkene"
    );
}

// -----------------------------------------------------------------------
// E/Z stereo transfer & creation in products (issue #50)
//
// Geometry is verified by comparing the product's canonical SMILES to the
// canonical SMILES of a reference molecule of known E/Z — exact, not "some
// Up/Down survives". (CipCode-based verification lives in the Python tests;
// chematic-rxn cannot depend on chematic-chem without a dependency cycle.)
// -----------------------------------------------------------------------

/// Canonical SMILES of the single product of `smirks` applied to `inputs`.
fn product_canon(smirks: &str, inputs: &[&str]) -> String {
    use chematic_smiles::canonical_smiles;
    let mols: Vec<Molecule> = inputs.iter().map(|s| parse(s).unwrap()).collect();
    let refs: Vec<&Molecule> = mols.iter().collect();
    let results = run_reactants(smirks, &refs).unwrap();
    assert!(!results.is_empty(), "no product for {smirks} on {inputs:?}");
    canonical_smiles(&results[0][0])
}

fn canon(smiles: &str) -> String {
    chematic_smiles::canonical_smiles(&parse(smiles).unwrap())
}

/// Writer-invariant E/Z of the first C=C double bond: `Some(true)` = E,
/// `Some(false)` = Z, `None` = no specified geometry. Reuses the crate's
/// own `ez_stereo_outward` (same convention as the #21 filter): equal
/// outward directions → Z, opposite → E.
fn double_bond_is_e(smiles: &str) -> Option<bool> {
    let mol = parse(smiles).unwrap();
    let (a1, a2) = mol
        .bonds()
        .find(|(_, b)| b.order == BondOrder::Double)
        .map(|(_, b)| (b.atom1, b.atom2))?;
    let sa = ez_stereo_outward(&mol, a1, a2)?;
    let sb = ez_stereo_outward(&mol, a2, a1)?;
    Some(sa != sb)
}

#[test]
fn issue50_transfer_identity_preserves_e() {
    // Identity SMIRKS on E-2-butene must yield an E product (was Z before Fix A).
    assert_eq!(
        product_canon("[C:1]=[C:2]>>[C:1]=[C:2]", &["C/C=C/C"]),
        canon("C/C=C/C"),
    );
}

#[test]
fn issue50_transfer_identity_preserves_z() {
    assert_eq!(
        product_canon("[C:1]=[C:2]>>[C:1]=[C:2]", &["C/C=C\\C"]),
        canon("C/C=C\\C"),
    );
}

#[test]
fn issue50_create_e_from_template() {
    // Product template introduces an E double bond from a saturated chain.
    assert_eq!(
        product_canon("[C:1][C:2][C:3][C:4]>>[C:1]/[C:2]=[C:3]/[C:4]", &["CCCC"]),
        canon("C/C=C/C"),
    );
}

#[test]
fn issue50_create_z_from_template() {
    assert_eq!(
        product_canon("[C:1][C:2][C:3][C:4]>>[C:1]/[C:2]=[C:3]\\[C:4]", &["CCCC"]),
        canon("C/C=C\\C"),
    );
}

#[test]
fn issue50_transfer_remote_reaction_keeps_e() {
    // Reaction at a remote site (aldehyde→alcohol) must not disturb the
    // E geometry of a carried-through alkene. The canonical writer may pick
    // `/C=C/` or `\C=C\` (both E) depending on traversal, so assert geometry
    // directly rather than the exact string.
    let got = product_canon("[CH:1]=O>>[C:1]O", &["CC/C=C/CC=O"]);
    assert_eq!(
        double_bond_is_e(&got),
        Some(true),
        "E geometry must survive a remote edit"
    );
    // And a Z input stays Z.
    let got_z = product_canon("[CH:1]=O>>[C:1]O", &["CC/C=C\\CC=O"]);
    assert_eq!(
        double_bond_is_e(&got_z),
        Some(false),
        "Z geometry must survive a remote edit"
    );
}

#[test]
fn issue50_geometry_is_deterministic() {
    // The pre-fix bug was nondeterministic (FxHashMap iteration order).
    // The same transform must give the same geometry on every run.
    let first = product_canon("[C:1]=[C:2]>>[C:1]=[C:2]", &["CC/C=C/CC"]);
    for _ in 0..6 {
        assert_eq!(
            product_canon("[C:1]=[C:2]>>[C:1]=[C:2]", &["CC/C=C/CC"]),
            first,
            "product geometry must be deterministic across runs"
        );
    }
}

#[test]
fn rdkit_9691_unspecified_product_ez_is_reactant_order_invariant() {
    let amine = "NC";
    for acid in ["C/C=C\\C(=O)O", "C/C=C/C(=O)O"] {
        let acid_first = product_canon(
            "[C:6][C:1]=[C:2][C:3](=O)[OH].[NH2:4][#6:5]>>[C:6][C:1]=[C:2][C:3](=O)[N:4][#6:5]",
            &[acid, amine],
        );
        let acid_second = product_canon(
            "[NH2:4][#6:5].[C:6][C:1]=[C:2][C:3](=O)[OH]>>[C:6][C:1]=[C:2][C:3](=O)[N:4][#6:5]",
            &[amine, acid],
        );
        let single_reactant = product_canon(
            "[C:6][C:1]=[C:2][C:3](=O)[OH]>>[C:6][C:1]=[C:2][C:3](=O)NC",
            &[acid],
        );

        assert_eq!(acid_first, canon("CC=CC(=O)NC"));
        assert_eq!(acid_second, acid_first);
        assert_eq!(single_reactant, acid_first);
        assert_eq!(double_bond_is_e(&acid_first), None);
    }
}

// -----------------------------------------------------------------------
// Reaction-transform performance regression witnesses (see
// docs/rfcs/reaction_transform_perf.md). Root cause: `chematic-smiles`'s
// `canonical_smiles()` wrote the winning individualize-refine branch's
// string, threw it away, and had `winning_individualized_ranks`'s caller
// write it a *second* time -- one fully redundant DFS-and-format pass on
// every single call, tied or not. Fixed by returning the already-written
// string instead of recomputing it. These three cases mirror the ones
// used to characterize and fix the regression:
// (a) a highly symmetric molecule (many individualize-refine branches --
//     this is the case that actually reproduces a large, measured slowdown
//     between chematic 0.4.25 and 0.4.30, NOT run_reactants match/product
//     volume, which stayed flat across versions);
// (b) an E/Z stereo control, since the fix touches the same
//     canonical-writer code path `resolve_ez_markers` depends on;
// (c) a negative control (asymmetric, no ties) that should show only the
//     universal (small, single-redundant-write) improvement, not the
//     symmetric-molecule-specific one.
// -----------------------------------------------------------------------

#[test]
fn perf_witness_a_symmetric_molecule_product_is_correct() {
    // Positive witness: adamantane (Td cage symmetry, 24
    // individualize-refine branches at time of writing) run through a
    // simple ring-opening SMIRKS. The fix must not change *which* string
    // wins -- only how many times it gets written -- so the product's
    // canonical SMILES must still round-trip to the same structure.
    let mol = parse("C1C2CC3CC1CC(C2)C3").unwrap(); // adamantane
    let results = run_reactants("[C:1][C:2]>>[C:1][C:2]", &[&mol]).unwrap();
    assert!(!results.is_empty(), "expected at least one C-C bond match");
    let canon = chematic_smiles::canonical_smiles(&results[0][0]);
    // Adamantane's own canonical form is a fixed point of this
    // identity-shaped SMIRKS: it must reparse to the exact same molecule
    // (same atom/bond count -- the transform doesn't add/remove atoms).
    let reparsed = parse(&canon).unwrap();
    assert_eq!(reparsed.atom_count(), mol.atom_count());
    assert_eq!(reparsed.bond_count(), mol.bond_count());
}

#[test]
fn perf_witness_b_ez_stereo_control_survives_symmetric_fix() {
    // Stereo control: identity/remote transforms on the E/Z pair used in
    // the perf investigation's own fixture
    // (crates/chematic-rxn/fixtures/witness_molecules.smi) must still
    // preserve exact geometry -- this is the non-negotiable issue #50
    // gate, re-run here against the specific molecules this perf fix
    // touched.
    assert_eq!(
        product_canon("[C:1]=[C:2]>>[C:1]=[C:2]", &["CC/C=C/CC(=O)O"]),
        canon("CC/C=C/CC(=O)O"),
        "(E)-hex-3-enoic acid must keep its E geometry"
    );
    assert_eq!(
        product_canon("[C:1]=[C:2]>>[C:1]=[C:2]", &["CC/C=C\\CC(=O)O"]),
        canon("CC/C=C\\CC(=O)O"),
        "(Z)-hex-3-enoic acid must keep its Z geometry"
    );
}

#[test]
fn perf_witness_c_negative_control_asymmetric_molecule() {
    // Negative control: aspirin has no automorphism ties (every ring
    // carbon has a distinct substitution environment), so
    // `winning_individualized_ranks` takes exactly one branch. This
    // exercises the same code path as (a) but should show only the
    // universal single-redundant-write saving, not a
    // symmetric-molecule-specific one -- included so a future reader can
    // tell the two effects apart empirically, not just by reasoning.
    let mol = parse("CC(=O)OC1=CC=CC=C1C(=O)O").unwrap(); // aspirin
    let results = run_reactants("[OH:1]-[C:2]=[O:3]>>C-[O:1]-[C:2]=[O:3]", &[&mol]).unwrap();
    assert!(!results.is_empty(), "expected the carboxylic acid to match");
    let canon = chematic_smiles::canonical_smiles(&results[0][0]);
    let reparsed = parse(&canon).unwrap();
    assert_eq!(reparsed.atom_count(), mol.atom_count() + 1); // +1 methyl carbon
}

// -------------------------------------------------------------------
// Match-level reaction application (issue #225)
// -------------------------------------------------------------------

/// `run_reactants(smirks, reactants)` must equal
/// `find_reaction_matches(...).filter_map(|m| apply_reaction_match(...))`
/// -- the exact equivalence issue #225's proposed API is built on.
/// Checked against several existing SMIRKS/reactant pairs already used
/// elsewhere in this file, not just one.
#[test]
fn find_and_apply_match_equals_run_reactants() {
    let mut cases: Vec<(&str, Vec<chematic_core::Molecule>)> = vec![
        ("[N:1]>>[N:1]", vec![parse("NCCN").unwrap()]),
        (
            "[N:1].[C:2]>>[N:1][C:2]",
            vec![parse("N").unwrap(), parse("C").unwrap()],
        ),
        ("[C:1][C:2]>>[C:1].[C:2]", vec![parse("CC").unwrap()]),
        (
            "[OH:1]-[C:2]=[O:3]>>C-[O:1]-[C:2]=[O:3]",
            vec![parse("CC(=O)OC1=CC=CC=C1C(=O)O").unwrap()],
        ),
    ];

    for (_, mols) in &mut cases {
        for mol in mols {
            for i in 0..mol.atom_count() {
                // Intentionally reused across reactants: tags are labels,
                // while atom_sources disambiguates reactant identity.
                mol.set_tag(AtomIdx(i as u32), Some(i as u16 + 1));
            }
        }
    }
    for (smirks, mols) in &cases {
        let reactants: Vec<&Molecule> = mols.iter().collect();
        let direct = run_reactants(smirks, &reactants).unwrap();

        let matches = find_reaction_matches(smirks, &reactants).unwrap();
        let via_matches: Vec<Vec<Molecule>> = matches
            .iter()
            .filter_map(|m| apply_reaction_match(smirks, &reactants, m, true).unwrap())
            .collect();

        assert_eq!(
            direct.len(),
            via_matches.len(),
            "{smirks}: product-set count must match"
        );
        for (d, v) in direct.iter().zip(via_matches.iter()) {
            assert_eq!(
                d.len(),
                v.len(),
                "{smirks}: product count per set must match"
            );
            for (dp, vp) in d.iter().zip(v.iter()) {
                for (idx, _) in dp.atoms() {
                    assert_eq!(dp.atom_tag(idx), vp.atom_tag(idx), "{smirks}: product tags");
                }
                assert_eq!(
                    chematic_smiles::canonical_smiles(dp),
                    chematic_smiles::canonical_smiles(vp),
                    "{smirks}: product molecule must match run_reactants exactly"
                );
            }
        }
    }
}

/// Issue #650: the traced application returns the same molecules as the
/// untraced one, and every traced product atom is a copy of the reactant
/// atom it names; mapped template atoms trace to their matched atoms and
/// template-created atoms trace to `None`.
#[test]
fn traced_apply_match_reports_atom_sources() {
    let mut cases: Vec<(&str, Vec<chematic_core::Molecule>)> = vec![
        ("[N:1]>>[N:1]", vec![parse("NCCN").unwrap()]),
        (
            "[N:1].[C:2]>>[N:1][C:2]",
            vec![parse("NCC").unwrap(), parse("CO").unwrap()],
        ),
        ("[C:1][C:2]>>[C:1].[C:2]", vec![parse("CCO").unwrap()]),
        // A product-only map labels a newly created atom; it has no
        // reactant origin, including when emitted as a second product.
        ("[C:1]>>[C:1][O:2]", vec![parse("C").unwrap()]),
        ("[C:1]>>[C:1].[O:2]", vec![parse("C").unwrap()]),
        (
            "[OH:1]-[C:2]=[O:3]>>C-[O:1]-[C:2]=[O:3]",
            vec![parse("CC(=O)OC1=CC=CC=C1C(=O)O").unwrap()],
        ),
        (
            "[C:1](=[O:2])[OH:3].[OH:4][C:5]>>[C:1](=[O:2])[O:4][C:5]",
            vec![parse("CC(=O)O").unwrap(), parse("OCC").unwrap()],
        ),
    ];
    for (_, mols) in &mut cases {
        for mol in mols {
            for i in 0..mol.atom_count() {
                mol.set_tag(AtomIdx(i as u32), Some(i as u16 + 1));
            }
        }
    }
    for (smirks, mols) in &cases {
        let reactants: Vec<&Molecule> = mols.iter().collect();
        let prepared = PreparedReaction::new(smirks).unwrap();
        let matches = find_reaction_matches(smirks, &reactants).unwrap();
        assert!(!matches.is_empty(), "{smirks}: expected a match");
        for m in &matches {
            for carry in [true, false] {
                let plain = prepared.apply_match(&reactants, m, carry).unwrap();
                let traced = prepared.apply_match_traced(&reactants, m, carry).unwrap();
                let free = apply_reaction_match_traced(smirks, &reactants, m, carry).unwrap();
                let (Some(plain), Some(traced), Some(free)) = (plain, traced, free) else {
                    continue;
                };
                assert_eq!(plain.len(), traced.len());
                assert_eq!(traced.len(), free.len());
                let global = m.atom_map_positions(smirks).unwrap();
                for (((p, t), f), template) in plain
                    .iter()
                    .zip(&traced)
                    .zip(&free)
                    .zip(&prepared.rxn.products)
                {
                    assert_eq!(
                        chematic_smiles::write(p),
                        chematic_smiles::write(&t.molecule),
                        "{smirks}: traced molecule differs"
                    );
                    assert_eq!(t.atom_sources, f.atom_sources);
                    assert_eq!(t.atom_sources.len(), t.molecule.atom_count());
                    for (i, src) in t.atom_sources.iter().enumerate() {
                        let expected_tag =
                            src.and_then(|src| reactants[src.reactant].atom_tag(src.atom));
                        for product in [p, &t.molecule, &f.molecule] {
                            assert_eq!(
                                product.atom_tag(AtomIdx(i as u32)),
                                expected_tag,
                                "{smirks}: tag at product atom {i}, carry={carry}"
                            );
                        }
                        if let Some(src) = src {
                            assert_eq!(
                                t.molecule.atom(AtomIdx(i as u32)).element,
                                reactants[src.reactant].atom(src.atom).element,
                                "{smirks}: atom {i} element"
                            );
                        }
                    }
                    // Sources are injective.
                    let mut seen = FxHashSet::default();
                    for src in t.atom_sources.iter().flatten() {
                        assert!(seen.insert(*src), "{smirks}: duplicate source {src:?}");
                    }
                    // Product template atom i is product atom i: a mapped
                    // one traces to where the match put that map number,
                    // an unmapped one is new.
                    for i in 0..template.atom_count() {
                        let expected = template
                            .atom(AtomIdx(i as u32))
                            .atom_map
                            .and_then(|am| global.get(&am))
                            .map(|&(reactant, atom)| ReactantAtom { reactant, atom });
                        assert_eq!(t.atom_sources[i], expected, "{smirks}: template atom {i}");
                    }
                }
            }
        }
    }
    // Template-created atom: the methyl C in the ester SMIRKS is new.
    let mol = parse("OC(=O)c1ccccc1").unwrap();
    let smirks = "[OH:1]-[C:2]=[O:3]>>C-[O:1]-[C:2]=[O:3]";
    let m = &find_reaction_matches(smirks, &[&mol]).unwrap()[0];
    let product = &apply_reaction_match_traced(smirks, &[&mol], m, true)
        .unwrap()
        .unwrap()[0];
    assert_eq!(
        product.atom_sources.iter().filter(|s| s.is_none()).count(),
        1
    );
    assert_eq!(
        product.atom_sources.iter().flatten().count(),
        mol.atom_count(),
        "every reactant atom is carried"
    );
}

#[test]
fn traced_report_preserves_product_order_and_diagnostics() {
    let cases = [
        ("[C:1]>>[C:1].[O:2]", vec![parse("C").unwrap()]),
        ("[#6:1]>>[#6:1]", vec![parse("CC").unwrap()]),
        (
            "[N:1].[C:2]>>[N:1][C:2]",
            vec![parse("N").unwrap(), parse("C").unwrap()],
        ),
        (
            "[13CH3:1][O:2]>>[13CH3:1].[O:2]",
            vec![parse("[13CH3]O").unwrap()],
        ),
    ];
    for (smirks, molecules) in cases {
        let refs = molecules.iter().collect::<Vec<_>>();
        let prepared = PreparedReaction::new(smirks).unwrap();
        let limits = ReactionTransformLimits::default();
        let plain = prepared
            .run_reactants_with_diagnostics(&refs, &limits)
            .unwrap();
        let traced = prepared
            .run_reactants_traced_with_diagnostics(&refs, &limits)
            .unwrap();
        assert_eq!(plain.diagnostics, traced.diagnostics, "{smirks}");
        assert_eq!(plain.products.len(), traced.products.len(), "{smirks}");
        for (plain_set, traced_set) in plain.products.iter().zip(&traced.products) {
            assert_eq!(plain_set.len(), traced_set.len(), "{smirks}");
            for (plain_product, traced_product) in plain_set.iter().zip(traced_set) {
                assert_eq!(
                    chematic_smiles::write(plain_product),
                    chematic_smiles::write(&traced_product.molecule),
                    "{smirks}"
                );
                assert_eq!(
                    traced_product.atom_sources.len(),
                    traced_product.molecule.atom_count(),
                    "{smirks}"
                );
                assert_eq!(
                    traced_product.template_maps.len(),
                    traced_product.molecule.atom_count(),
                    "{smirks}"
                );
            }
        }
        if smirks == "[C:1]>>[C:1].[O:2]" {
            assert_eq!(traced.products[0][0].template_maps, vec![Some(1)]);
            assert_eq!(traced.products[0][1].template_maps, vec![Some(2)]);
            assert_eq!(traced.products[0][1].atom_sources, vec![None]);
        }
    }
}

#[test]
fn symmetric_ether_cleavage_retains_both_source_assignments() {
    let mut ether = parse("COC").unwrap();
    for i in 0..ether.atom_count() {
        ether.set_tag(AtomIdx(i as u32), Some(i as u16 + 1));
    }
    let smirks = "[C:1][O:2][C:3]>>[C:1][O:2].[C:3]";
    let products = run_reactants(smirks, &[&ether]).unwrap();
    let mut assignments: Vec<_> = products
        .iter()
        .map(|set| {
            assert_eq!(set.len(), 2);
            (
                set[0].atom_tag(AtomIdx(0)).unwrap().get(),
                set[1].atom_tag(AtomIdx(0)).unwrap().get(),
            )
        })
        .collect();
    assignments.sort_unstable();
    assert_eq!(assignments, [(1, 3), (3, 1)]);
}

/// `run_reactants_strict` (carry_substituents=false) must also compose
/// the same way as the `carry_substituents=true` case above.
#[test]
fn find_and_apply_match_equals_run_reactants_strict() {
    let mol = parse("CC(=O)OC1=CC=CC=C1C(=O)O").unwrap();
    let reactants = [&mol];
    let smirks = "[OH:1]-[C:2]=[O:3]>>C-[O:1]-[C:2]=[O:3]";

    let direct = run_reactants_strict(smirks, &reactants).unwrap();
    let matches = find_reaction_matches(smirks, &reactants).unwrap();
    let via_matches: Vec<Vec<Molecule>> = matches
        .iter()
        .filter_map(|m| apply_reaction_match(smirks, &reactants, m, false).unwrap())
        .collect();

    assert_eq!(direct.len(), via_matches.len());
    for (d, v) in direct.iter().zip(via_matches.iter()) {
        for (dp, vp) in d.iter().zip(v.iter()) {
            assert_eq!(
                chematic_smiles::canonical_smiles(dp),
                chematic_smiles::canonical_smiles(vp)
            );
        }
    }
}

/// The core motivating use case from issue #225: enumerate matches
/// independently of applying them, reject some based on a property of
/// the match itself, and apply only the accepted ones -- without
/// discarding the legitimate matches along with the rejected one.
#[test]
fn selective_match_application() {
    let mol = parse("NCCN").unwrap();
    let reactants = [&mol];
    let smirks = "[N:1]>>[N:1]";

    let matches = find_reaction_matches(smirks, &reactants).unwrap();
    assert_eq!(matches.len(), 2, "two N atoms in NCCN → two matches");

    // Reject the match touching the higher-numbered atom (arbitrary
    // match-specific policy standing in for RENKIN's real ring-bond
    // rejection rule), keep the other.
    let positions: Vec<_> = matches
        .iter()
        .map(|m| m.atom_map_positions(smirks).unwrap()[&1])
        .collect();
    let keep_idx = if positions[0].1.0 < positions[1].1.0 {
        0
    } else {
        1
    };

    let applied = apply_reaction_match(smirks, &reactants, &matches[keep_idx], true).unwrap();
    assert!(applied.is_some(), "the accepted match must still apply");

    // Applying only one match must yield exactly one of the two
    // products `run_reactants` would have returned for the whole call,
    // not both and not neither.
    let full = run_reactants(smirks, &reactants).unwrap();
    assert_eq!(full.len(), 2);
    let applied_canon = chematic_smiles::canonical_smiles(&applied.unwrap()[0]);
    assert!(
        full.iter()
            .any(|ps| chematic_smiles::canonical_smiles(&ps[0]) == applied_canon),
        "the selectively-applied product must be one of run_reactants's own outputs"
    );
}

/// `ReactionMatch::atom_map_positions` must resolve atom_map:1 to the
/// actual matched N atom, for each of the two NCCN matches separately.
#[test]
fn atom_map_positions_resolves_matched_atom() {
    let mol = parse("NCCN").unwrap();
    let reactants = [&mol];
    let smirks = "[N:1]>>[N:1]";

    let matches = find_reaction_matches(smirks, &reactants).unwrap();
    let mut matched_atoms: Vec<AtomIdx> = matches
        .iter()
        .map(|m| {
            let positions = m.atom_map_positions(smirks).unwrap();
            let (reactant_slot, atom_idx) = positions[&1];
            assert_eq!(reactant_slot, 0, "single-reactant SMIRKS: slot must be 0");
            assert_eq!(mol.atom(atom_idx).element.symbol(), "N");
            atom_idx
        })
        .collect();
    matched_atoms.sort_by_key(|a| a.0);
    assert_eq!(
        matched_atoms,
        vec![AtomIdx(0), AtomIdx(3)],
        "NCCN's two N atoms are at index 0 and 3"
    );
}

/// `apply_reaction_match` must return `Ok(None)` -- not an error and not
/// a product -- for a match whose product set fails the existing
/// valence filter, matching the case [`run_reactants`] silently drops
/// (`overvalent_product_filtered_oxygen` above).
#[test]
fn apply_reaction_match_none_on_valence_violation() {
    let ethanol = parse("CCO").unwrap();
    let reactants = [&ethanol];
    let smirks = "[O:1]>>[O:1](C)C";

    let matches = find_reaction_matches(smirks, &reactants).unwrap();
    assert_eq!(matches.len(), 1, "exactly one O in ethanol");

    let applied = apply_reaction_match(smirks, &reactants, &matches[0], true).unwrap();
    assert!(
        applied.is_none(),
        "over-valenced product must come back as Ok(None), not Some(..) or Err(..)"
    );
}

#[test]
fn diagnostics_classify_valence_filtered_match_without_changing_products() {
    let ethanol = parse("CCO").unwrap();
    let smirks = "[O:1]>>[O:1](C)C";
    let report =
        run_reactants_with_diagnostics(smirks, &[&ethanol], &ReactionTransformLimits::default())
            .unwrap();
    assert_eq!(report.products.len(), 0);
    assert_eq!(report.diagnostics.accepted_matches, 1);
    assert_eq!(report.diagnostics.applied_products, 0);
    assert_eq!(report.diagnostics.valence_rejected_matches, 1);
    assert!(!report.diagnostics.truncated_matches);
    assert!(run_reactants(smirks, &[&ethanol]).unwrap().is_empty());
}

#[test]
fn detailed_diagnostics_report_typed_reason_and_mapped_atom_context() {
    let ethanol = parse("CCO").unwrap();
    let outcome = run_reactants_traced_rdkit_2026_03_6_detailed(
        "[O:1]>>[O:1](C)C",
        &[&ethanol],
        &ReactionTransformLimits::default(),
    )
    .unwrap();
    let RdkitDetailedProfileOutcome::Report(report) = outcome else {
        panic!("ordinary valence rejection must be diagnosable");
    };
    assert_eq!(report.rejected_products.len(), 1);
    let diagnostic = &report.rejection_diagnostics[0][0];
    assert!(!diagnostic.accepted);
    assert_eq!(diagnostic.reason, Some(ProductRejectionReason::Valence));
    assert_eq!(diagnostic.atoms.len(), 1);
    let atom = &diagnostic.atoms[0];
    assert_eq!(atom.template_map, Some(1));
    assert_eq!(atom.element, chematic_core::Element::O);
    assert_eq!(atom.observed_valence, Some(3));
    assert_eq!(atom.max_allowed_valence, Some(2));
}

#[test]
fn detailed_diagnostics_distinguish_aromaticity_and_kekulization() {
    let mut outside_ring = parse("c1ccccc1").unwrap();
    outside_ring.remove_bond(BondIdx(0));
    let outside_ring = TracedProduct {
        atom_sources: vec![None; outside_ring.atom_count()],
        template_maps: vec![None; outside_ring.atom_count()],
        molecule: outside_ring,
    };
    let aromaticity = product_sanitization_diagnostic(&outside_ring, 0);
    assert_eq!(
        aromaticity.reason,
        Some(ProductRejectionReason::Aromaticity)
    );
    assert!(!aromaticity.atoms.is_empty());

    let mut odd_aromatic_ring = parse("C1CCCC1").unwrap();
    for atom in 0..odd_aromatic_ring.atom_count() {
        odd_aromatic_ring.set_atom_aromatic(AtomIdx(atom as u32), true);
    }
    for bond in 0..odd_aromatic_ring.bond_count() {
        odd_aromatic_ring.set_bond_order(BondIdx(bond as u32), BondOrder::Aromatic);
    }
    let odd_aromatic_ring = TracedProduct {
        atom_sources: vec![None; odd_aromatic_ring.atom_count()],
        template_maps: vec![None; odd_aromatic_ring.atom_count()],
        molecule: odd_aromatic_ring,
    };
    let kekulization = product_sanitization_diagnostic(&odd_aromatic_ring, 0);
    assert_eq!(
        kekulization.reason,
        Some(ProductRejectionReason::Kekulization)
    );
    assert_eq!(kekulization.atoms.len(), 5);
}

#[test]
fn mapped_isotope_template_matches_only_the_specified_nuclide() {
    let labeled = parse("[13CH3]O").unwrap();
    let unlabeled = parse("CO").unwrap();
    let smirks = "[13CH3:1][O:2]>>[13CH3:1][O:2]";

    let products = run_reactants(smirks, &[&labeled]).unwrap();
    assert_eq!(products.len(), 1);
    assert_eq!(products[0].len(), 1);
    assert!(
        products[0]
            .iter()
            .any(|m| m.atoms().any(|(_, a)| a.isotope == Some(13)))
    );
    assert!(run_reactants(smirks, &[&unlabeled]).unwrap().is_empty());
}

#[test]
fn diagnostics_match_ring_aware_products_and_counts() {
    let mol = parse("C1CCCCC1").unwrap();
    let rings = chematic_perception::find_sssr(&mol);
    let prepared = PreparedReaction::new("[C:1]>>[C:1]").unwrap();
    let report = prepared
        .run_reactants_with_rings_and_limits_with_diagnostics(
            &[&mol],
            &[&rings],
            &ReactionTransformLimits::default(),
        )
        .unwrap();
    assert_eq!(
        report.diagnostics.accepted_matches,
        report.diagnostics.applied_products
    );
    assert_eq!(report.diagnostics.valence_rejected_matches, 0);
    assert_eq!(report.products.len(), report.diagnostics.applied_products);
}

#[test]
fn variant_diagnostics_match_ring_aware_application() {
    let mol = parse("C1CCCCC1").unwrap();
    let rings = chematic_perception::find_sssr(&mol);
    let prepared = PreparedReaction::new("[#6:1]>>[#6:1]").unwrap();
    let reports = prepared
        .run_reactants_with_rings_and_limits_with_variant_diagnostics(
            &[&mol],
            &[&rings],
            &ReactionTransformLimits::default(),
        )
        .unwrap();
    assert_eq!(reports.len(), 2);
    assert_eq!(reports[0].variant_index, 0);
    assert!(reports[0].diagnostics.accepted_matches > 0);
    assert_eq!(
        reports
            .iter()
            .map(|report| report.diagnostics.applied_products)
            .sum::<usize>(),
        prepared
            .run_reactants_with_rings_and_limits(
                &[&mol],
                &[&rings],
                &ReactionTransformLimits::default(),
            )
            .unwrap()
            .len()
    );
}

/// `find_reaction_matches` and `apply_reaction_match` must propagate
/// the same `ReactantCountMismatch` error `run_reactants` does when the
/// number of input molecules doesn't match the SMIRKS's reactant-slot
/// count.
#[test]
fn find_and_apply_match_reactant_count_mismatch_errors() {
    let mol = parse("C").unwrap();
    let smirks = "[N:1].[C:2]>>[N:1][C:2]";

    let find_err = find_reaction_matches(smirks, &[&mol]);
    assert!(matches!(
        find_err,
        Err(TransformError::ReactantCountMismatch {
            expected: 2,
            got: 1
        })
    ));

    let dummy_match = ReactionMatch {
        per_reactant: vec![FxHashMap::default()],
    };
    let apply_err = apply_reaction_match(smirks, &[&mol], &dummy_match, true);
    assert!(matches!(
        apply_err,
        Err(TransformError::ReactantCountMismatch {
            expected: 2,
            got: 1
        })
    ));
}

/// `apply_reaction_match` must also reject a `ReactionMatch` whose own
/// `per_reactant` shape doesn't match the SMIRKS being applied against
/// (e.g. a match obtained from a different SMIRKS), rather than
/// panicking on an out-of-bounds index into `template_atom_maps`.
#[test]
fn apply_reaction_match_rejects_mismatched_match_shape() {
    let n_mol = parse("N").unwrap();
    let c_mol = parse("C").unwrap();
    let reactants = [&n_mol, &c_mol];
    let smirks = "[N:1].[C:2]>>[N:1][C:2]";

    // A match shaped for a single-reactant SMIRKS, applied against a
    // two-reactant one.
    let mismatched_match = ReactionMatch {
        per_reactant: vec![FxHashMap::default()],
    };
    let err = apply_reaction_match(smirks, &reactants, &mismatched_match, true);
    assert!(matches!(
        err,
        Err(TransformError::ReactantCountMismatch {
            expected: 2,
            got: 1
        })
    ));
}

/// `ReactionMatch::atom_map_positions` must likewise reject a
/// shape mismatch rather than panicking.
#[test]
fn atom_map_positions_rejects_mismatched_match_shape() {
    let mismatched_match = ReactionMatch {
        per_reactant: vec![FxHashMap::default()],
    };
    let err = mismatched_match.atom_map_positions("[N:1].[C:2]>>[N:1][C:2]");
    assert!(matches!(
        err,
        Err(TransformError::ReactantCountMismatch {
            expected: 2,
            got: 1
        })
    ));
}

#[test]
fn atom_tag_survives_apply_atom_map_does_not() {
    let mut mol = parse("CCO").unwrap();
    for i in 0..mol.atom_count() {
        mol.set_tag(AtomIdx(i as u32), Some(200 + i as u16));
    }
    let matches = find_reaction_matches("[C:1]>>[C:1]O", &[&mol]).expect("match");
    assert!(!matches.is_empty());
    let products =
        apply_reaction_match("[C:1]>>[C:1]O", &[&mol], &matches[0], true).expect("apply");
    let product = &products.expect("valence")[0];
    let tags: Vec<_> = product
        .atoms()
        .map(|(i, _)| product.atom_tag(i).map(core::num::NonZeroU16::get))
        .collect();
    assert!(tags.contains(&Some(200)) || tags.contains(&Some(201)) || tags.contains(&Some(202)));
    assert!(product.atoms().all(|(_, a)| a.atom_map.is_none()));
    // Born oxygen has no tag.
    assert!(
        product
            .atoms()
            .any(|(i, a)| product.atom_tag(i).is_none() && a.element.atomic_number() == 8)
    );
}

#[test]
fn native_detailed_report_exposes_typed_rejection_without_returning_product() {
    let reactant = parse("C[C@H](O)c1ccccc1").unwrap();
    let report = run_reactants_traced_with_rejection_diagnostics(
        "[C:1]-[CH:2](-[OH:3])-[c:4]>>[C:1]-[C:2](=[O:3])-[c:4]",
        &[&reactant],
        &ReactionTransformLimits { max_matches: 100 },
    )
    .unwrap();

    assert_eq!(report.report.diagnostics.accepted_matches, 1);
    assert_eq!(report.report.diagnostics.applied_products, 0);
    assert_eq!(report.report.diagnostics.valence_rejected_matches, 1);
    assert!(report.report.products.is_empty());
    assert_eq!(report.rejection_diagnostics.len(), 1);
    let rejected = &report.rejection_diagnostics[0];
    assert_eq!(rejected.len(), 1);
    assert!(!rejected[0].accepted);
    assert_eq!(rejected[0].reason, Some(ProductRejectionReason::Valence));
    assert!(!rejected[0].atoms.is_empty());
    assert!(
        rejected[0]
            .atoms
            .iter()
            .any(|atom| atom.template_map == Some(2))
    );
}
