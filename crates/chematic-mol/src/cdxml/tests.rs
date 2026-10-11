use super::*;
use chematic_core::Chirality;

// Minimal hand-crafted CDXML for ethanol (C-C-O) with 2 bonds.
const ETHANOL_CDXML: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd">
<CDXML>
<fragment>
<n id="1" p="10.0 20.0" Element="6" NumHydrogens="3"/>
<n id="2" p="25.0 20.0" Element="6" NumHydrogens="2"/>
<n id="3" p="40.0 20.0" Element="8" NumHydrogens="1"/>
<b B="1" E="2" Order="1"/>
<b B="2" E="3" Order="1"/>
</fragment>
</CDXML>"#;

#[test]
fn parse_cdxml_ethanol_atom_count() {
    let (mol, coords) = parse_cdxml(ETHANOL_CDXML).unwrap();
    assert_eq!(mol.atom_count(), 3, "ethanol: 3 heavy atoms");
    assert_eq!(mol.bond_count(), 2, "ethanol: 2 bonds");
    assert_eq!(coords.len(), 3, "one coord per atom");
}

#[test]
fn parse_cdxml_ethanol_elements() {
    let (mol, _) = parse_cdxml(ETHANOL_CDXML).unwrap();
    let elems: Vec<&str> = mol.atoms().map(|(_, a)| a.element.symbol()).collect();
    assert!(elems.contains(&"C"), "should contain C");
    assert!(elems.contains(&"O"), "should contain O");
}

#[test]
fn parse_cdxml_ethanol_coords() {
    let (_, coords) = parse_cdxml(ETHANOL_CDXML).unwrap();
    // Atom 1 (first C): p="10.0 20.0"
    assert!(
        (coords[0].0 - 10.0).abs() < 0.01,
        "first atom x=10.0: {:?}",
        coords[0]
    );
    assert!(
        (coords[0].1 - 20.0).abs() < 0.01,
        "first atom y=20.0: {:?}",
        coords[0]
    );
}

#[test]
fn parse_cdxml_carbon_element_6() {
    let cdxml = r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
</fragment></CDXML>"#;
    let (mol, _) = parse_cdxml(cdxml).unwrap();
    assert_eq!(mol.atom_count(), 1);
    let atom = mol.atom(chematic_core::AtomIdx(0));
    assert_eq!(atom.element.symbol(), "C", "Element=6 → Carbon");
}

#[test]
fn parse_cdxml_double_bond() {
    let cdxml = r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="8" p="10 0"/>
<b B="1" E="2" Order="2"/>
</fragment></CDXML>"#;
    let (mol, _) = parse_cdxml(cdxml).unwrap();
    let bond = mol.bond(chematic_core::BondIdx(0));
    assert_eq!(bond.order, BondOrder::Double, "Order=2 → Double");
}

#[test]
fn parse_cdxml_charge() {
    let cdxml = r#"<CDXML><fragment>
<n id="1" Element="7" Charge="1" p="0 0"/>
</fragment></CDXML>"#;
    let (mol, _) = parse_cdxml(cdxml).unwrap();
    let atom = mol.atom(chematic_core::AtomIdx(0));
    assert_eq!(atom.charge, 1, "Charge=1 → N+");
}

#[test]
fn parse_cdxml_unknown_atomic_number_returns_err() {
    let cdxml = r#"<CDXML><fragment>
<n id="1" Element="999" p="0 0"/>
</fragment></CDXML>"#;
    let result = parse_cdxml(cdxml);
    assert!(
        matches!(result, Err(CdxmlError::UnknownAtomicNumber(_))),
        "unknown atomic number should return Err"
    );
}

#[test]
fn parse_cdxml_element_above_255_is_rejected() {
    // Element=300 would silently truncate to 44 (Ru) via u32 as u8.
    // Must be caught as UnknownAtomicNumber before the cast.
    let cdxml = r#"<CDXML><fragment>
<n id="1" Element="300" p="0 0"/>
</fragment></CDXML>"#;
    let result = parse_cdxml(cdxml);
    assert!(
        matches!(result, Err(CdxmlError::UnknownAtomicNumber(300))),
        "Element=300 must return UnknownAtomicNumber(300)"
    );
}

// B4: multi-fragment CDXML document tests

const TWO_FRAGMENT_CDXML: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<CDXML>
<fragment>
<n id="1" p="10.0 20.0" Element="6"/>
<n id="2" p="25.0 20.0" Element="8"/>
<b B="1" E="2" Order="1"/>
</fragment>
<fragment>
<n id="3" p="60.0 20.0" Element="7"/>
<n id="4" p="75.0 20.0" Element="6"/>
<n id="5" p="90.0 20.0" Element="6"/>
<b B="3" E="4" Order="1"/>
<b B="4" E="5" Order="1"/>
</fragment>
</CDXML>"#;

#[test]
fn parse_cdxml_all_two_fragments_count() {
    let mols = parse_cdxml_all(TWO_FRAGMENT_CDXML).unwrap();
    assert_eq!(mols.len(), 2, "two <fragment> elements → two molecules");
}

#[test]
fn parse_cdxml_all_first_fragment_co() {
    let mols = parse_cdxml_all(TWO_FRAGMENT_CDXML).unwrap();
    let (mol, _) = &mols[0];
    assert_eq!(mol.atom_count(), 2, "first fragment: C + O");
    assert_eq!(mol.bond_count(), 1);
}

#[test]
fn parse_cdxml_all_second_fragment_ncc() {
    let mols = parse_cdxml_all(TWO_FRAGMENT_CDXML).unwrap();
    let (mol, _) = &mols[1];
    assert_eq!(mol.atom_count(), 3, "second fragment: N + C + C");
    assert_eq!(mol.bond_count(), 2);
}

#[test]
fn parse_cdxml_all_coords_independent() {
    let mols = parse_cdxml_all(TWO_FRAGMENT_CDXML).unwrap();
    let (_, coords0) = &mols[0];
    let (_, coords1) = &mols[1];
    // First fragment atom 0 at x=10
    assert!((coords0[0].0 - 10.0).abs() < 0.01);
    // Second fragment atom 0 at x=60
    assert!((coords1[0].0 - 60.0).abs() < 0.01);
}

#[test]
fn parse_cdxml_empty_doc_returns_empty_vec() {
    let cdxml = r#"<?xml version="1.0"?><CDXML></CDXML>"#;
    let mols = parse_cdxml_all(cdxml).unwrap();
    assert!(mols.is_empty(), "empty CDXML → empty Vec");
}

// -----------------------------------------------------------------------
// Order="1.5" aromatic bond (OpenBabel / some CDXML writers)
// -----------------------------------------------------------------------

#[test]
fn parse_cdxml_aromatic_bond_order_1_5() {
    // Some CDXML producers (e.g. tools derived from OpenBabel) write
    // aromatic bonds as Order="1.5".  These must be stored as
    // BondOrder::Aromatic rather than falling through to Single.
    let cdxml = r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="6" p="10 0"/>
<b B="1" E="2" Order="1.5"/>
</fragment></CDXML>"#;
    let (mol, _) = parse_cdxml(cdxml).unwrap();
    let bond = mol.bond(chematic_core::BondIdx(0));
    assert_eq!(bond.order, BondOrder::Aromatic, "Order=1.5 → Aromatic");
}

#[test]
fn parse_cdxml_benzene_all_aromatic_bonds() {
    // Benzene written with all Order="1.5" bonds.
    let cdxml = r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="6" p="10 0"/>
<n id="3" Element="6" p="20 0"/>
<n id="4" Element="6" p="30 0"/>
<n id="5" Element="6" p="20 10"/>
<n id="6" Element="6" p="10 10"/>
<b B="1" E="2" Order="1.5"/>
<b B="2" E="3" Order="1.5"/>
<b B="3" E="4" Order="1.5"/>
<b B="4" E="5" Order="1.5"/>
<b B="5" E="6" Order="1.5"/>
<b B="6" E="1" Order="1.5"/>
</fragment></CDXML>"#;
    let (mol, _) = parse_cdxml(cdxml).unwrap();
    assert_eq!(mol.atom_count(), 6, "benzene has 6 atoms");
    assert_eq!(mol.bond_count(), 6, "benzene has 6 bonds");
    let all_aromatic = mol.bonds().all(|(_, b)| b.order == BondOrder::Aromatic);
    assert!(all_aromatic, "all bonds must be Aromatic");
}

// -----------------------------------------------------------------------
// write_cdxml round-trip
// -----------------------------------------------------------------------

#[test]
fn write_cdxml_roundtrip_ethanol() {
    let (mol, coords) = parse_cdxml(ETHANOL_CDXML).unwrap();
    let written = write_cdxml(&mol, &coords);
    let (mol2, coords2) = parse_cdxml(&written).unwrap();

    assert_eq!(mol2.atom_count(), mol.atom_count());
    assert_eq!(mol2.bond_count(), mol.bond_count());
    let elems: Vec<&str> = mol.atoms().map(|(_, a)| a.element.symbol()).collect();
    let elems2: Vec<&str> = mol2.atoms().map(|(_, a)| a.element.symbol()).collect();
    assert_eq!(elems, elems2);
    for (c1, c2) in coords.iter().zip(coords2.iter()) {
        assert!((c1.0 - c2.0).abs() < 0.01);
        assert!((c1.1 - c2.1).abs() < 0.01);
    }
}

#[test]
fn write_cdxml_roundtrip_double_bond() {
    let cdxml = r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="8" p="10 0"/>
<b B="1" E="2" Order="2"/>
</fragment></CDXML>"#;
    let (mol, coords) = parse_cdxml(cdxml).unwrap();
    let written = write_cdxml(&mol, &coords);
    let (mol2, _) = parse_cdxml(&written).unwrap();
    let bond2 = mol2.bond(chematic_core::BondIdx(0));
    assert_eq!(bond2.order, BondOrder::Double);
}

#[test]
fn write_cdxml_roundtrip_wedge_bond() {
    let cdxml = r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="6" p="10 0"/>
<b B="1" E="2" Order="1" Display="WedgeBegin"/>
</fragment></CDXML>"#;
    let (mol, coords) = parse_cdxml(cdxml).unwrap();
    let written = write_cdxml(&mol, &coords);
    assert!(written.contains("Display=\"WedgeBegin\""));
    let (mol2, _) = parse_cdxml(&written).unwrap();
    let bond2 = mol2.bond(chematic_core::BondIdx(0));
    assert_eq!(bond2.order, BondOrder::Up);
}

/// Issue found while surveying RDKit's open issues (analogous to RDKit
/// #9359, "CDXML reading doesn't use Bold or undirectional hash bonds
/// for stereochemistry"): ChemDraw's simplified non-directional "Bold"
/// bond display (a plain thick line, used by some chemists in place of
/// a real wedge) must be interpreted as stereo, same as this reader
/// already does for bare "Hash".
#[test]
fn parse_cdxml_bold_bond_is_wedge_up() {
    let cdxml = r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="6" p="10 0"/>
<b B="1" E="2" Order="1" Display="Bold"/>
</fragment></CDXML>"#;
    let (mol, _) = parse_cdxml(cdxml).unwrap();
    let bond = mol.bond(chematic_core::BondIdx(0));
    assert_eq!(
        bond.order,
        BondOrder::Up,
        "Display=\"Bold\" must be read as a wedge-up stereo bond, not silently \
             dropped to a plain single bond"
    );
}

/// "Bold" and "WedgeBegin" both encode the same "coming toward viewer"
/// stereo relative to a bond's B->E atom order -- for identical
/// connectivity/geometry, both display attributes must parse to the
/// same BondOrder.
#[test]
fn parse_cdxml_bold_bond_matches_wedge_begin() {
    let bold = r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="6" p="10 0"/>
<b B="1" E="2" Order="1" Display="Bold"/>
</fragment></CDXML>"#;
    let wedge = r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="6" p="10 0"/>
<b B="1" E="2" Order="1" Display="WedgeBegin"/>
</fragment></CDXML>"#;
    let (mol_bold, _) = parse_cdxml(bold).unwrap();
    let (mol_wedge, _) = parse_cdxml(wedge).unwrap();
    assert_eq!(
        mol_bold.bond(chematic_core::BondIdx(0)).order,
        mol_wedge.bond(chematic_core::BondIdx(0)).order,
    );
}

/// Negative control: bare "Hash"/"Dash" (already-handled non-directional
/// stereo) and a plain undecorated bond must be unaffected by adding
/// the "Bold" arm -- confirms no accidental widening of the match.
#[test]
fn parse_cdxml_bold_addition_does_not_change_other_display_values() {
    for (display, expected) in [
        ("Hash", BondOrder::Down),
        ("Dash", BondOrder::Down),
        ("WedgeEnd", BondOrder::Down),
        ("WedgedHashEnd", BondOrder::Down),
        ("WedgedHashBegin", BondOrder::Up),
    ] {
        let cdxml = format!(
            r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="6" p="10 0"/>
<b B="1" E="2" Order="1" Display="{display}"/>
</fragment></CDXML>"#
        );
        let (mol, _) = parse_cdxml(&cdxml).unwrap();
        assert_eq!(
            mol.bond(chematic_core::BondIdx(0)).order,
            expected,
            "Display=\"{display}\" regressed"
        );
    }
    let plain = r#"<CDXML><fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="6" p="10 0"/>
<b B="1" E="2" Order="1"/>
</fragment></CDXML>"#;
    let (mol, _) = parse_cdxml(plain).unwrap();
    assert_eq!(
        mol.bond(chematic_core::BondIdx(0)).order,
        BondOrder::Single,
        "a bond with no Display attribute at all must stay plain Single"
    );
}

// ── RDKit issue #9359: actual tetrahedral chirality perception ────────
//
// The tests above only ever checked `bond.order` -- they confirm the
// *bond* survives as a wedge, not that a *chirality* was ever perceived
// from it. Before this fix, `flush()` never called any wedge-to-parity
// mechanism at all, for ANY Display value including plain WedgeBegin --
// so `Atom.chirality` was always `Chirality::None` for every CDXML
// molecule, no matter how it was drawn. These tests assert the actual
// RDKit #9359 expectation (`center_atom(mol).GetChiralTag() !=
// CHI_UNSPECIFIED`, i.e. `Atom.chirality != Chirality::None` here), and
// then go further to check *meaning*, not just non-None.

/// RDKit #9359's exact repro fixture (`BOLD_BASE_CDXML`, atom/bond ids
/// preserved), parameterised on the `Display` of the C1-F2 bond (`101`)
/// so the no-wedge/Bold/Hash/WedgeBegin cases share one source of
/// truth. Atom 1 (implicit Carbon, the default when `Element` is
/// omitted) is the only degree-4 atom -- F(2), the ring-continuation
/// C(3)->C(4), Cl(7), and CH3(6) -- matching RDKit's own
/// `center_atom` (`atomic_number==6 and degree==4`) selector.
fn rdkit_9359_cdxml(display: Option<&str>) -> String {
    let bond101 = match display {
        Some(d) => format!(r#"<b id="101" B="1" E="2" Display="{d}"/>"#),
        None => r#"<b id="101" B="1" E="2"/>"#.to_string(),
    };
    format!(
        r#"<CDXML>
<page>
<fragment>
<n id="1" p="270.73 355.40"/>
<n id="2" p="263.44 367.76" Element="9"/>
<n id="3" p="283.22 348.33"/>
<n id="4" p="295.58 355.62"/>
<n id="6" p="258.37 348.11"/>
<n id="7" p="277.79 367.89" Element="17" NumHydrogens="0"><t><s>Cl</s></t></n>
{bond101}
<b id="102" B="1" E="3"/>
<b id="103" B="1" E="7"/>
<b id="104" B="1" E="6"/>
<b id="105" B="3" E="4"/>
</fragment>
</page>
</CDXML>"#
    )
}

const RDKIT_9359_CENTER: chematic_core::AtomIdx = chematic_core::AtomIdx(0);

/// Parses with [`CdxmlParseOptions::infer_nondirectional_stereo`] on --
/// the opt-in path, used by tests specifically exercising Bold/Hash
/// chirality perception (which is off by default, see the
/// default-behavior tests below).
fn parse_infer(cdxml: &str) -> (Molecule, Vec<(f64, f64)>) {
    parse_cdxml_with_options(
        cdxml,
        &CdxmlParseOptions {
            infer_nondirectional_stereo: true,
        },
    )
    .unwrap()
}

#[test]
fn rdkit_9359_no_display_has_no_chirality() {
    // Flat drawing, no wedge at all -- must NOT invent stereo, with or
    // without the opt-in (nothing to infer from either way).
    let (mol, _) = parse_cdxml(&rdkit_9359_cdxml(None)).unwrap();
    assert_eq!(mol.atom(RDKIT_9359_CENTER).chirality, Chirality::None);
    let (mol, _) = parse_infer(&rdkit_9359_cdxml(None));
    assert_eq!(mol.atom(RDKIT_9359_CENTER).chirality, Chirality::None);
}

// ── Default behavior: Bold/Hash do NOT perceive chirality unopted ─────
//
// `Bold` in particular is sometimes drawn by chemists as a plain thick
// line for visual emphasis, not stereo intent (RDKit #9359's own
// discussion leans toward this being opt-in, not silently default) --
// so `parse_cdxml`/`parse_cdxml_all` (no explicit options) must NOT
// manufacture a chirality from it. `bond.order` still faithfully
// records `BondOrder::Up`/`Down` either way (structural fidelity is
// unconditional) -- only chirality PERCEPTION is gated.

#[test]
fn rdkit_9359_bold_default_does_not_perceive_chirality() {
    let (mol, _) = parse_cdxml(&rdkit_9359_cdxml(Some("Bold"))).unwrap();
    assert_eq!(
        mol.atom(RDKIT_9359_CENTER).chirality,
        Chirality::None,
        "Bold must NOT perceive chirality by default (opt-in required)"
    );
    assert_eq!(
        mol.bond(chematic_core::BondIdx(0)).order,
        BondOrder::Up,
        "the bond's own order must still faithfully record Up regardless \
             of whether chirality perception is opted in"
    );
}

#[test]
fn rdkit_9359_hash_default_does_not_perceive_chirality() {
    let (mol, _) = parse_cdxml(&rdkit_9359_cdxml(Some("Hash"))).unwrap();
    assert_eq!(
        mol.atom(RDKIT_9359_CENTER).chirality,
        Chirality::None,
        "Hash must NOT perceive chirality by default (opt-in required)"
    );
    assert_eq!(mol.bond(chematic_core::BondIdx(0)).order, BondOrder::Down);
}

/// Directional wedges are NOT gated by the opt-in -- they're
/// unambiguous stereo intent (a real tapered wedge), so the default
/// (no options) already perceives chirality from them, matching the
/// RDKit #9359 expectation (`center_atom(mol).GetChiralTag() !=
/// CHI_UNSPECIFIED`) out of the box.
#[test]
fn rdkit_9359_wedge_begin_default_perceives_chirality() {
    let (mol, _) = parse_cdxml(&rdkit_9359_cdxml(Some("WedgeBegin"))).unwrap();
    assert_ne!(
        mol.atom(RDKIT_9359_CENTER).chirality,
        Chirality::None,
        "a real directional wedge must perceive chirality without any opt-in"
    );
    assert!(mol.stereo_neighbor_order(RDKIT_9359_CENTER).is_some());
}

// ── Opt-in behavior (`infer_nondirectional_stereo: true`) ─────────────

#[test]
fn rdkit_9359_bold_perceives_chirality_when_opted_in() {
    let (mol, _) = parse_infer(&rdkit_9359_cdxml(Some("Bold")));
    assert_ne!(
        mol.atom(RDKIT_9359_CENTER).chirality,
        Chirality::None,
        "RDKit #9359: Display=\"Bold\" must perceive a real chirality when \
             opted in, matching center_atom(mol).GetChiralTag() != CHI_UNSPECIFIED"
    );
    assert!(mol.stereo_neighbor_order(RDKIT_9359_CENTER).is_some());
}

#[test]
fn rdkit_9359_hash_perceives_chirality_when_opted_in() {
    let (mol, _) = parse_infer(&rdkit_9359_cdxml(Some("Hash")));
    assert_ne!(
        mol.atom(RDKIT_9359_CENTER).chirality,
        Chirality::None,
        "RDKit #9359: Display=\"Hash\" (undirectional) must also perceive \
             chirality when opted in"
    );
}

/// The meaning check the reviewer required, not just non-None: Bold
/// ("coming toward viewer") and Hash ("going away from viewer") are
/// chemically OPPOSITE for the exact same 2D layout -- if perception
/// were just guessing, there'd be no reason for these to differ let
/// alone invert consistently.
#[test]
fn rdkit_9359_bold_and_hash_are_opposite_configurations() {
    let (mol_bold, _) = parse_infer(&rdkit_9359_cdxml(Some("Bold")));
    let (mol_hash, _) = parse_infer(&rdkit_9359_cdxml(Some("Hash")));
    let bold_chirality = mol_bold.atom(RDKIT_9359_CENTER).chirality;
    let hash_chirality = mol_hash.atom(RDKIT_9359_CENTER).chirality;
    assert_ne!(bold_chirality, Chirality::None);
    assert_ne!(hash_chirality, Chirality::None);
    assert_ne!(
        bold_chirality, hash_chirality,
        "Bold (toward viewer) and Hash (away from viewer) on the identical \
             2D layout must give opposite configurations, not the same one"
    );

    // Independent oracle cross-check: accurate CIP assignment must
    // agree that these two are opposite, not just chematic's own raw
    // Chirality enum by coincidence.
    let cip_bold = chematic_chem::assign_cip(&mol_bold).get(RDKIT_9359_CENTER);
    let cip_hash = chematic_chem::assign_cip(&mol_hash).get(RDKIT_9359_CENTER);
    assert!(
        cip_bold.is_some() && cip_hash.is_some(),
        "{cip_bold:?} {cip_hash:?}"
    );
    assert_ne!(
        cip_bold, cip_hash,
        "independent CIP oracle must also see Bold/Hash as opposite R/S"
    );
}

/// "Bold" and "WedgeBegin" both encode the same "coming toward viewer"
/// convention -- for the identical layout they must produce the SAME
/// chirality, not just both-non-None independently.
#[test]
fn rdkit_9359_bold_matches_wedge_begin_configuration() {
    let (mol_bold, _) = parse_infer(&rdkit_9359_cdxml(Some("Bold")));
    let (mol_wedge, _) = parse_cdxml(&rdkit_9359_cdxml(Some("WedgeBegin"))).unwrap();
    assert_eq!(
        mol_bold.atom(RDKIT_9359_CENTER).chirality,
        mol_wedge.atom(RDKIT_9359_CENTER).chirality,
    );
}

/// B/E-reversed variant, non-directional (required by review, and the
/// review's own required *expected value* corrected from an earlier,
/// wrong assumption): `Bold`/`Hash`/`Dash` have NO Begin/End
/// convention in the CDXML spec at all -- unlike `WedgeBegin`/
/// `WedgeEnd`, whose Begin atom is a real, drawn narrow end. Which
/// atom a `<b>` element happens to list first as `B` for a `Bold`
/// bond is an arbitrary artifact of how the file was written, not a
/// stereo-relevant signal -- so swapping B/E for the SAME `Bold`
/// drawing must PRESERVE the perceived configuration, not invert it.
#[test]
fn rdkit_9359_bold_be_reversed_preserves_configuration() {
    let (mol_forward, _) = parse_infer(&rdkit_9359_cdxml(Some("Bold")));
    let reversed = rdkit_9359_cdxml(Some("Bold")).replace(
        r#"<b id="101" B="1" E="2" Display="Bold"/>"#,
        r#"<b id="101" B="2" E="1" Display="Bold"/>"#,
    );
    let (mol_reversed, _) = parse_infer(&reversed);
    let forward = mol_forward.atom(RDKIT_9359_CENTER).chirality;
    let rev = mol_reversed.atom(RDKIT_9359_CENTER).chirality;
    assert_ne!(forward, Chirality::None);
    assert_eq!(
        forward, rev,
        "Bold has no Begin/End convention -- reversing B/E for the identical \
             drawing must NOT change the perceived configuration"
    );
}

/// Same B/E-reversal check for `Hash` (the other non-directional
/// display), independently -- not inferred from the Bold case alone.
#[test]
fn rdkit_9359_hash_be_reversed_preserves_configuration() {
    let (mol_forward, _) = parse_infer(&rdkit_9359_cdxml(Some("Hash")));
    let reversed = rdkit_9359_cdxml(Some("Hash")).replace(
        r#"<b id="101" B="1" E="2" Display="Hash"/>"#,
        r#"<b id="101" B="2" E="1" Display="Hash"/>"#,
    );
    let (mol_reversed, _) = parse_infer(&reversed);
    let forward = mol_forward.atom(RDKIT_9359_CENTER).chirality;
    let rev = mol_reversed.atom(RDKIT_9359_CENTER).chirality;
    assert_ne!(forward, Chirality::None);
    assert_eq!(
        forward, rev,
        "Hash has no Begin/End convention -- reversing B/E for the identical \
             drawing must NOT change the perceived configuration"
    );
}

/// Negative control, directional: `WedgeBegin` DOES have a real
/// Begin/End convention (the Begin atom is the drawn narrow end), so
/// switching to `WedgeEnd` for the SAME B/E atom order must genuinely
/// invert the configuration (the reference atom flips) -- proves the
/// preserve-under-swap behavior above is specific to the
/// non-directional set, not a blanket "wedges don't encode direction"
/// regression. Uses the default parser (no opt-in needed, directional
/// wedges are unaffected by the flag either way).
#[test]
fn rdkit_9359_wedge_begin_vs_wedge_end_inverts_configuration() {
    let (mol_begin, _) = parse_cdxml(&rdkit_9359_cdxml(Some("WedgeBegin"))).unwrap();
    let (mol_end, _) = parse_cdxml(&rdkit_9359_cdxml(Some("WedgeEnd"))).unwrap();
    let begin = mol_begin.atom(RDKIT_9359_CENTER).chirality;
    let end = mol_end.atom(RDKIT_9359_CENTER).chirality;
    assert_ne!(begin, Chirality::None);
    assert_ne!(end, Chirality::None);
    assert_ne!(
        begin, end,
        "WedgeBegin vs WedgeEnd (same B/E order) must invert -- the Begin \
             atom is a real, drawn reference, unlike Bold/Hash"
    );
}

/// Reflection (negate all Y coordinates) must invert the perceived
/// configuration; rotation + translation must leave it unchanged.
#[test]
fn rdkit_9359_bold_reflection_inverts_rotation_translation_preserves() {
    let base = rdkit_9359_cdxml(Some("Bold"));
    let (mol_base, coords_base) = parse_infer(&base);
    let base_chirality = mol_base.atom(RDKIT_9359_CENTER).chirality;
    assert_ne!(base_chirality, Chirality::None);

    // Reflection: negate Y for every atom's `p` attribute.
    let reflect_re = regex_lite_reflect_y(&base);
    let (mol_reflected, _) = parse_infer(&reflect_re);
    assert_eq!(
        mol_reflected.atom(RDKIT_9359_CENTER).chirality,
        invert(base_chirality),
        "mirroring the 2D layout (negate Y) must invert the configuration"
    );

    // Rotation (37 degrees) + translation (+1000, -500): must be a
    // no-op for chirality -- recompute coords directly (not via the
    // CDXML string) since only chematic_perception's own math needs
    // exercising here, not the parser's coordinate formatting.
    let theta: f64 = 37.0_f64.to_radians();
    let (sin_t, cos_t) = theta.sin_cos();
    let rotated_translated: Vec<(f64, f64)> = coords_base
        .iter()
        .map(|&(x, y)| {
            let (rx, ry) = (x * cos_t - y * sin_t, x * sin_t + y * cos_t);
            (rx + 1000.0, ry - 500.0)
        })
        .collect();
    let mut mol_rt = mol_base.clone();
    chematic_perception::apply_local_parity_from_wedges(&mut mol_rt, &rotated_translated);
    assert_eq!(
        mol_rt.atom(RDKIT_9359_CENTER).chirality,
        base_chirality,
        "rotation + translation must not change the perceived configuration"
    );
}

fn invert(c: Chirality) -> Chirality {
    match c {
        Chirality::Clockwise => Chirality::CounterClockwise,
        Chirality::CounterClockwise => Chirality::Clockwise,
        Chirality::None => Chirality::None,
        other => other,
    }
}

/// Negates the Y coordinate of every atom's `p="x y"` attribute in a
/// CDXML fragment string -- test-only helper, not a general XML editor.
fn regex_lite_reflect_y(cdxml: &str) -> String {
    let mut out = String::with_capacity(cdxml.len());
    let mut rest = cdxml;
    while let Some(start) = rest.find("p=\"") {
        out.push_str(&rest[..start + 3]);
        rest = &rest[start + 3..];
        let end = rest.find('"').expect("unterminated p attribute");
        let coord = &rest[..end];
        let mut parts = coord.split_whitespace();
        let x: f64 = parts.next().unwrap().parse().unwrap();
        let y: f64 = parts.next().unwrap().parse().unwrap();
        out.push_str(&format!("{x} {}", -y));
        rest = &rest[end..];
    }
    out.push_str(rest);
    out
}

// ── Ambiguity: non-directional wedge between two stereocenters ────────

/// Both endpoints of a single "Bold" bond independently qualify as
/// stereocenter candidates (degree 4 each) -- a non-directional wedge
/// has no narrow/wide-end convention to say which one it actually
/// describes, so BOTH must abstain (fail closed) rather than each
/// independently claiming a configuration derived from contradictory
/// height assignments of the same bond.
#[test]
fn bold_bond_between_two_stereocenters_abstains_both_sides() {
    let cdxml = r#"<CDXML><fragment>
<n id="1" p="0 0"/>
<n id="2" p="10 0"/>
<n id="3" p="0 10" Element="9"/>
<n id="4" p="-10 0" Element="17"/>
<n id="5" p="0 -10" Element="35"/>
<n id="6" p="20 10" Element="9"/>
<n id="7" p="20 -10" Element="17"/>
<n id="8" p="10 -20" Element="35"/>
<b B="1" E="2" Display="Bold"/>
<b B="1" E="3"/>
<b B="1" E="4"/>
<b B="1" E="5"/>
<b B="2" E="6"/>
<b B="2" E="7"/>
<b B="2" E="8"/>
</fragment></CDXML>"#;
    let (mol, _) = parse_infer(cdxml);
    let a1 = chematic_core::AtomIdx(0);
    let a2 = chematic_core::AtomIdx(1);
    assert_eq!(mol.neighbors(a1).count(), 4);
    assert_eq!(mol.neighbors(a2).count(), 4);
    assert_eq!(
        mol.atom(a1).chirality,
        Chirality::None,
        "ambiguous non-directional wedge between two stereocenters must abstain"
    );
    assert_eq!(mol.atom(a2).chirality, Chirality::None);
    assert_eq!(
        mol.bond(chematic_core::BondIdx(0)).order,
        BondOrder::Up,
        "the bond's own order still faithfully records Bold -> Up; only \
             chirality PERCEPTION abstains, structural fidelity is unconditional"
    );
}

/// Negative control: the exact same connectivity/geometry, but with a
/// real directional "WedgeBegin" instead of "Bold" -- the Begin atom is
/// unambiguously the reference by CDXML's own convention, so this is
/// NOT ambiguous and must keep perceiving chirality normally (proves
/// the downgrade in the test above is specific to non-directional
/// displays, not a blanket "shared bond between two stereocenters"
/// suppression that would also break the common, unambiguous case).
#[test]
fn wedge_begin_between_two_stereocenters_is_not_ambiguous() {
    let cdxml = r#"<CDXML><fragment>
<n id="1" p="0 0"/>
<n id="2" p="10 0"/>
<n id="3" p="0 10" Element="9"/>
<n id="4" p="-10 0" Element="17"/>
<n id="5" p="0 -10" Element="35"/>
<n id="6" p="20 10" Element="9"/>
<n id="7" p="20 -10" Element="17"/>
<n id="8" p="10 -20" Element="35"/>
<b B="1" E="2" Display="WedgeBegin"/>
<b B="1" E="3"/>
<b B="1" E="4"/>
<b B="1" E="5"/>
<b B="2" E="6"/>
<b B="2" E="7"/>
<b B="2" E="8"/>
</fragment></CDXML>"#;
    let (mol, _) = parse_cdxml(cdxml).unwrap();
    let a1 = chematic_core::AtomIdx(0);
    assert_ne!(
        mol.atom(a1).chirality,
        Chirality::None,
        "a directional WedgeBegin between two stereocenters is not ambiguous \
             and must still perceive chirality"
    );
    assert_eq!(
        mol.bond(chematic_core::BondIdx(0)).order,
        BondOrder::Up,
        "a directional wedge must never be downgraded by the ambiguity check"
    );
}

#[test]
fn write_cdxml_roundtrip_charge_isotope_hcount() {
    let cdxml = r#"<CDXML><fragment>
<n id="1" Element="7" Charge="1" Isotope="15" NumHydrogens="2" p="0 0"/>
</fragment></CDXML>"#;
    let (mol, coords) = parse_cdxml(cdxml).unwrap();
    let written = write_cdxml(&mol, &coords);
    let (mol2, _) = parse_cdxml(&written).unwrap();
    let atom2 = mol2.atom(chematic_core::AtomIdx(0));
    assert_eq!(atom2.charge, 1);
    assert_eq!(atom2.isotope, Some(15));
    assert_eq!(atom2.hydrogen_count, Some(2));
}

#[test]
fn bounded_parser_rejects_input_line_and_atom_limits() {
    assert!(matches!(
        parse_cdxml_with_limits(
            ETHANOL_CDXML,
            &CdxmlParseLimits {
                max_input_bytes: 8,
                ..Default::default()
            }
        ),
        Err(CdxmlError::ResourceLimit {
            resource: "input bytes",
            ..
        })
    ));
    assert!(matches!(
        parse_cdxml_with_limits(
            ETHANOL_CDXML,
            &CdxmlParseLimits {
                max_lines: 1,
                ..Default::default()
            }
        ),
        Err(CdxmlError::ResourceLimit {
            resource: "lines",
            ..
        })
    ));
    assert!(matches!(
        parse_cdxml_with_limits(
            ETHANOL_CDXML,
            &CdxmlParseLimits {
                max_atoms: 1,
                ..Default::default()
            }
        ),
        Err(CdxmlError::ResourceLimit {
            resource: "atoms",
            ..
        })
    ));
}

#[test]
fn bounded_parser_rejects_bond_and_fragment_limits() {
    assert!(matches!(
        parse_cdxml_all_with_limits(
            ETHANOL_CDXML,
            &CdxmlParseLimits {
                max_bonds: 1,
                ..Default::default()
            }
        ),
        Err(CdxmlError::ResourceLimit {
            resource: "bonds",
            ..
        })
    ));
    assert!(matches!(
        parse_cdxml_all_with_limits(
            TWO_FRAGMENT_CDXML,
            &CdxmlParseLimits {
                max_fragments: 1,
                ..Default::default()
            }
        ),
        Err(CdxmlError::ResourceLimit {
            resource: "fragments",
            ..
        })
    ));
}
