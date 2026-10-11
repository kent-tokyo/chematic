use super::*;
use chematic_core::{AtomIdx, BondOrder, Element};

// -----------------------------------------------------------------------
// Test fixtures
// -----------------------------------------------------------------------

const METHANE_V3K: &str = "\
methane
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
";

const ETHANOL_V3K: &str = "\
ethanol
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  V30 2 C 1.5 0.0 0.0 0
M  V30 3 O 3.0 0.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
";

// -----------------------------------------------------------------------
// Test 1: methane — 1 atom, 0 bonds
// -----------------------------------------------------------------------
#[test]
fn test_methane_counts() {
    let (mol, _) = parse_mol_v3000(METHANE_V3K).expect("parse methane");
    assert_eq!(mol.atom_count(), 1);
    assert_eq!(mol.bond_count(), 0);
}

#[test]
fn rejects_malformed_coordinates_and_records_after_terminator() {
    assert!(parse_mol_v3000(&METHANE_V3K.replace("0.0 0.0 0.0 0", "nope 0.0 0.0 0")).is_err());
    assert!(parse_mol_v3000(&(METHANE_V3K.to_owned() + "M  V30 BEGIN CTAB\n")).is_err());
    assert!(parse_mol_v3000(&METHANE_V3K.replace("M  V30 END ATOM", "M V30 END ATOM")).is_err());
}

// -----------------------------------------------------------------------
// Test 2: ethanol — 3 atoms, 2 bonds
// -----------------------------------------------------------------------
#[test]
fn test_ethanol_counts() {
    let (mol, _) = parse_mol_v3000(ETHANOL_V3K).expect("parse ethanol");
    assert_eq!(mol.atom_count(), 3);
    assert_eq!(mol.bond_count(), 2);
}

// -----------------------------------------------------------------------
// Test 3: ethanol — bond between atoms 0 and 1 is Single
// -----------------------------------------------------------------------
#[test]
fn test_ethanol_bond_0_1_single() {
    let (mol, _) = parse_mol_v3000(ETHANOL_V3K).expect("parse ethanol");
    let (_, bond) = mol
        .bond_between(AtomIdx(0), AtomIdx(1))
        .expect("bond 0-1 exists");
    assert_eq!(bond.order, BondOrder::Single);
}

// -----------------------------------------------------------------------
// Test 4: ethanol — bond between atoms 1 and 2 is Single
// -----------------------------------------------------------------------
#[test]
fn test_ethanol_bond_1_2_single() {
    let (mol, _) = parse_mol_v3000(ETHANOL_V3K).expect("parse ethanol");
    let (_, bond) = mol
        .bond_between(AtomIdx(1), AtomIdx(2))
        .expect("bond 1-2 exists");
    assert_eq!(bond.order, BondOrder::Single);
}

// -----------------------------------------------------------------------
// Test 5: CHG=1 → atom.charge == 1
// -----------------------------------------------------------------------
#[test]
fn test_positive_charge() {
    let mol_str = "\
charged_pos
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N 0.0 0.0 0.0 0 CHG=1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
";
    let (mol, _) = parse_mol_v3000(mol_str).expect("parse charged_pos");
    assert_eq!(mol.atom(AtomIdx(0)).charge, 1);
}

// -----------------------------------------------------------------------
// Test 6: CHG=-1 → atom.charge == -1
// -----------------------------------------------------------------------
#[test]
fn test_negative_charge() {
    let mol_str = "\
charged_neg
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O 0.0 0.0 0.0 0 CHG=-1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
";
    let (mol, _) = parse_mol_v3000(mol_str).expect("parse charged_neg");
    assert_eq!(mol.atom(AtomIdx(0)).charge, -1);
}

// -----------------------------------------------------------------------
// Test 7: MASS=13 → atom.isotope == Some(13)
// -----------------------------------------------------------------------
#[test]
fn test_isotope() {
    let mol_str = "\
isotope
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0 MASS=13
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
";
    let (mol, _) = parse_mol_v3000(mol_str).expect("parse isotope");
    assert_eq!(mol.atom(AtomIdx(0)).isotope, Some(13));
}

// -----------------------------------------------------------------------
// Test 8: bond type=4 → BondOrder::Aromatic
// -----------------------------------------------------------------------
#[test]
fn test_aromatic_bond() {
    let mol_str = "\
aromatic
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  V30 2 C 1.5 0.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 4 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
";
    let (mol, _) = parse_mol_v3000(mol_str).expect("parse aromatic");
    let (_, bond) = mol
        .bond_between(AtomIdx(0), AtomIdx(1))
        .expect("bond exists");
    assert_eq!(bond.order, BondOrder::Aromatic);
}

// -----------------------------------------------------------------------
// Test 9: bond type=2 → BondOrder::Double
// -----------------------------------------------------------------------
#[test]
fn test_double_bond() {
    let mol_str = "\
double_bond
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  V30 2 O 1.2 0.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
";
    let (mol, _) = parse_mol_v3000(mol_str).expect("parse double_bond");
    let (_, bond) = mol
        .bond_between(AtomIdx(0), AtomIdx(1))
        .expect("bond exists");
    assert_eq!(bond.order, BondOrder::Double);
}

// -----------------------------------------------------------------------
// Test 10: metadata — name and comment parsed correctly
// -----------------------------------------------------------------------
#[test]
fn test_metadata() {
    let mol_str = "\
my_molecule
  some_prog
my comment line
  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
";
    let (_, meta) = parse_mol_v3000(mol_str).expect("parse metadata");
    assert_eq!(meta.name, "my_molecule");
    assert_eq!(meta.comment, "my comment line");
}

// -----------------------------------------------------------------------
// Test 11: line continuation (atom line split with trailing `-`)
// -----------------------------------------------------------------------
#[test]
fn test_line_continuation() {
    // The atom line for atom 1 is split after the z-coordinate.
    // The HCOUNT keyword appears on the continuation line.
    let mol_str = "\
continuation
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0 MASS=12 -
M  V30 HCOUNT=3
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
";
    let (mol, _) = parse_mol_v3000(mol_str).expect("parse continuation");
    // atom 0 should have isotope 12 and hydrogen_count 3 (set via continuation)
    assert_eq!(mol.atom(AtomIdx(0)).element, Element::C);
    assert_eq!(mol.atom(AtomIdx(0)).isotope, Some(12));
    assert_eq!(mol.atom(AtomIdx(0)).hydrogen_count, Some(3));
}

// -----------------------------------------------------------------------
// Test 12: missing END ATOM → returns Err
// -----------------------------------------------------------------------
#[test]
fn test_missing_end_atom_is_error() {
    let mol_str = "\
bad
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  END
";
    let result = parse_mol_v3000(mol_str);
    assert!(
        matches!(result, Err(MolParseError::V3000ParseError { .. })),
        "expected V3000ParseError but got a different result"
    );
}

// -----------------------------------------------------------------------
// Test 13: elements parsed correctly for ethanol (C, C, O)
// -----------------------------------------------------------------------
#[test]
fn test_ethanol_elements() {
    let (mol, _) = parse_mol_v3000(ETHANOL_V3K).expect("parse ethanol");
    let atoms: Vec<_> = mol.atoms().collect();
    assert_eq!(atoms[0].1.element, Element::C);
    assert_eq!(atoms[1].1.element, Element::C);
    assert_eq!(atoms[2].1.element, Element::O);
}

// -----------------------------------------------------------------------
// Test 14: triple bond (type=3 → BondOrder::Triple)
// -----------------------------------------------------------------------
#[test]
fn test_triple_bond() {
    let mol_str = "\
triple_bond
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  V30 2 N 1.2 0.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 3 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
";
    let (mol, _) = parse_mol_v3000(mol_str).expect("parse triple_bond");
    let (_, bond) = mol
        .bond_between(AtomIdx(0), AtomIdx(1))
        .expect("bond exists");
    assert_eq!(bond.order, BondOrder::Triple);
}

#[test]
fn test_query_bond_types_preserved() {
    let mol_str = "\
query_bonds
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0
M  V30 2 C 1.0 0.0 0.0 0
M  V30 3 C 2.0 0.0 0.0 0
M  V30 4 C 3.0 0.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 5 1 2
M  V30 2 8 3 4
M  V30 END BOND
M  V30 END CTAB
M  END
";
    let (mol, meta) = parse_mol_v3000(mol_str).expect("parse query bonds");
    let bonds: Vec<_> = mol.bonds().collect();
    assert_eq!(bonds[0].1.order, BondOrder::QuerySingleOrDouble);
    assert_eq!(bonds[1].1.order, BondOrder::QueryAny);

    let written = write_mol_v3000(&mol, &meta, &[]);
    assert!(written.contains("M  V30 1 5 1 2"), "{written}");
    assert!(written.contains("M  V30 2 8 3 4"), "{written}");
}

// -----------------------------------------------------------------------
// Test: MDL bond type 9 (dative/coordinate) round-trips, not silently
// corrupted to Single -- regression test for the platinum coordination-
// chemistry benchmark (validation/platinum/FEASIBILITY.md). The molblock
// below is not hand-crafted: it is exactly what RDKit 2026.03.3 writes
// for `Chem.MolFromSmiles` + `AddBond(n, pt, Chem.BondType.DATIVE)` +
// `Chem.MolToMolBlock(mol, forceV3000=True)` -- RDKit auto-upgrades to
// V3000 for any dative bond (it cannot express one in V2000 either).
// Before this fix, bond type 9 fell through the reader's `_ =>
// BondOrder::Single` catch-all with no error, no warning: a
// structurally valid file silently produced a different molecule.
// -----------------------------------------------------------------------
#[test]
fn test_dative_bond_type_9_round_trips() {
    // NOT built with a `"\` line-continuation: the real molblock's second
    // header line is blank-then-indented (RDKit writes no molecule name),
    // and a leading `\`-continuation would silently swallow that blank
    // line along with its indentation (rustc even warns "multiple lines
    // skipped by escaped newline" if written that way) -- explicit `\n`
    // keeps this byte-for-byte identical to RDKit's own output.
    let mol_str = "\n     RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 3 2 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 N 0.000000 0.000000 0.000000 0\nM  V30 2 Pt 1.299038 0.750000 0.000000 0 VAL=2\nM  V30 3 Cl 2.598076 1.500000 0.000000 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 9 1 2\nM  V30 2 1 2 3\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";
    let (mol, meta) = parse_mol_v3000(mol_str).expect("parse RDKit dative molblock");
    let bonds: Vec<_> = mol.bonds().collect();
    assert_eq!(
        bonds[0].1.order,
        BondOrder::Dative,
        "bond type 9 must read as Dative, not silently fall through to Single"
    );
    // atom1/atom2 preserve the file's donor/acceptor order (N -> Pt).
    assert_eq!(mol.atom(bonds[0].1.atom1).element, Element::N);
    assert_eq!(mol.atom(bonds[0].1.atom2).element, Element::PT);
    assert_eq!(bonds[1].1.order, BondOrder::Single);

    let written = write_mol_v3000(&mol, &meta, &[]);
    assert!(
        written.contains("M  V30 1 9 1 2"),
        "Dative must write back out as bond type 9: {written}"
    );

    // Full round trip: re-parsing the freshly-written molblock still
    // gives a Dative bond with the same donor/acceptor order.
    let (rt_mol, _) = parse_mol_v3000(&written).expect("re-parse written molblock");
    let rt_bonds: Vec<_> = rt_mol.bonds().collect();
    assert_eq!(rt_bonds[0].1.order, BondOrder::Dative);
    assert_eq!(rt_mol.atom(rt_bonds[0].1.atom1).element, Element::N);
    assert_eq!(rt_mol.atom(rt_bonds[0].1.atom2).element, Element::PT);
}

#[test]
fn shared_mol_v3000_roundtrip_contract_matches() {
    let document: serde_json::Value = serde_json::from_str(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/cross_binding_contract.json"
    )))
    .expect("contract JSON");
    let contract = &document["mol_v3000_contract"];
    let (mol, meta) = parse_mol_v3000(contract["input"].as_str().unwrap()).unwrap();
    let serialized = write_mol_v3000(&mol, &meta, &[]);
    let (roundtripped, _) = parse_mol_v3000(&serialized).unwrap();
    assert_eq!(
        roundtripped.atom_count(),
        contract["expected"]["atom_count"]
    );
    assert_eq!(roundtripped.bond_count(), 1);
}

// -----------------------------------------------------------------------
// Test 15: atom-map number stored when nonzero
// -----------------------------------------------------------------------
#[test]
fn test_atom_map() {
    let mol_str = "\
atommapped
  test

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 3
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
";
    let (mol, _) = parse_mol_v3000(mol_str).expect("parse atommapped");
    assert_eq!(mol.atom(AtomIdx(0)).atom_map, Some(3));
}
