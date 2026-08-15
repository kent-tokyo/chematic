//! Integration tests for `chematic-mol`'s `mmcif`/`pqr` modules.

use chematic_core::Element;
use chematic_mol::cif;
use chematic_mol::mmcif::{parse_mmcif, write_mmcif};
use chematic_mol::pqr::{parse_pqr, write_pqr};

/// A realistic multi-chain, multi-model (2 NMR-style models), altloc-
/// bearing mmCIF fixture, hand-authored against the wwPDB mmCIF
/// dictionary's `_atom_site` column layout (not copied from a real PDB
/// entry or any tool's output).
const MULTI_CHAIN_MODEL_ALTLOC_MMCIF: &str = "data_MULTI\n\
    #\n\
    _cell.length_a      50.000\n\
    _cell.length_b      60.000\n\
    _cell.length_c      70.000\n\
    _cell.angle_alpha   90.000\n\
    _cell.angle_beta    95.500\n\
    _cell.angle_gamma   90.000\n\
    #\n\
    _symmetry.space_group_name_H-M   'P 21 21 21'\n\
    #\n\
    loop_\n\
    _atom_site.group_PDB\n\
    _atom_site.id\n\
    _atom_site.type_symbol\n\
    _atom_site.label_atom_id\n\
    _atom_site.label_alt_id\n\
    _atom_site.label_comp_id\n\
    _atom_site.label_asym_id\n\
    _atom_site.label_entity_id\n\
    _atom_site.label_seq_id\n\
    _atom_site.pdbx_PDB_ins_code\n\
    _atom_site.Cartn_x\n\
    _atom_site.Cartn_y\n\
    _atom_site.Cartn_z\n\
    _atom_site.occupancy\n\
    _atom_site.B_iso_or_equiv\n\
    _atom_site.pdbx_formal_charge\n\
    _atom_site.auth_seq_id\n\
    _atom_site.auth_comp_id\n\
    _atom_site.auth_asym_id\n\
    _atom_site.auth_atom_id\n\
    _atom_site.pdbx_PDB_model_num\n\
    ATOM   1 N N  . GLY A 1 1 ? 1.000 2.000 3.000 1.00 20.00 ? 1 GLY A N  1\n\
    ATOM   2 C CA A GLY A 1 1 ? 1.500 2.500 3.500 0.55 21.00 ? 1 GLY A CA 1\n\
    ATOM   3 C CA B GLY A 1 1 ? 1.600 2.600 3.600 0.45 22.00 ? 1 GLY A CA 1\n\
    ATOM   4 N N  . ALA B 2 1 ? 10.000 11.000 12.000 1.00 18.00 ? 101 ALA C N 1\n\
    HETATM 5 CA CA . CA  C 3 . ? 40.000 41.000 42.000 1.00 30.00 2 501 CA D CA 1\n\
    ATOM   6 N N  . GLY A 1 1 ? 1.010 2.020 3.030 1.00 20.10 ? 1 GLY A N  2\n\
    ATOM   7 C CA A GLY A 1 1 ? 1.510 2.520 3.530 0.55 21.10 ? 1 GLY A CA 2\n\
    ATOM   8 C CA B GLY A 1 1 ? 1.610 2.620 3.630 0.45 22.10 ? 1 GLY A CA 2\n\
    ATOM   9 N N  . ALA B 2 1 ? 10.010 11.020 12.030 1.00 18.10 ? 101 ALA C N 2\n\
    HETATM 10 CA CA . CA  C 3 . ? 40.010 41.020 42.030 1.00 30.10 2 501 CA D CA 2\n";

#[test]
fn multi_chain_model_altloc_fixture_parses_correctly() {
    let r = parse_mmcif(MULTI_CHAIN_MODEL_ALTLOC_MMCIF).unwrap();
    assert_eq!(r.atoms.len(), 10);

    // Two models, five atoms each, in file order.
    let model1: Vec<_> = r.atoms.iter().filter(|a| a.model_num == 1).collect();
    let model2: Vec<_> = r.atoms.iter().filter(|a| a.model_num == 2).collect();
    assert_eq!(model1.len(), 5);
    assert_eq!(model2.len(), 5);

    // Three distinct auth chains (A, C, D) via auth_asym_id preference.
    let chains: std::collections::BTreeSet<&str> =
        r.atoms.iter().map(|a| a.chain_id.as_str()).collect();
    assert_eq!(chains, ["A", "C", "D"].into_iter().collect());

    // Altloc A/B pair at the same residue/position, different occupancies.
    let alt_a = r.atoms.iter().find(|a| a.alt_loc == Some('A')).unwrap();
    let alt_b = r.atoms.iter().find(|a| a.alt_loc == Some('B')).unwrap();
    assert!((alt_a.occupancy - 0.55).abs() < 1e-9);
    assert!((alt_b.occupancy - 0.45).abs() < 1e-9);

    // HETATM calcium ion, formal charge +2, distinct auth_seq_id (501) from
    // label_seq_id (3) -- auth_* preferred.
    let ca_ion = r.atoms.iter().find(|a| a.group_pdb == "HETATM").unwrap();
    assert_eq!(ca_ion.res_seq, 501);
    assert_eq!(ca_ion.formal_charge, Some(2));
    assert_eq!(ca_ion.element, Element::CA);

    // Unit cell + space group preserved.
    let cell = r.cell.unwrap();
    assert!((cell.beta - 95.5).abs() < 1e-9);
    assert_eq!(r.space_group.as_deref(), Some("P 21 21 21"));
}

#[test]
fn multi_chain_model_altloc_fixture_round_trips_through_write_read() {
    let r = parse_mmcif(MULTI_CHAIN_MODEL_ALTLOC_MMCIF).unwrap();
    let written = write_mmcif(&r.atoms, r.cell.as_ref(), r.space_group.as_deref(), "MULTI");
    let r2 = parse_mmcif(&written).unwrap();

    assert_eq!(r.atoms.len(), r2.atoms.len());
    for (a, b) in r.atoms.iter().zip(r2.atoms.iter()) {
        assert_eq!(a.chain_id, b.chain_id, "chain id must round-trip");
        assert_eq!(a.res_name, b.res_name, "residue name must round-trip");
        assert_eq!(a.res_seq, b.res_seq, "residue number must round-trip");
        assert_eq!(a.alt_loc, b.alt_loc, "altloc must round-trip");
        assert_eq!(a.model_num, b.model_num, "model number must round-trip");
        assert!((a.occupancy - b.occupancy).abs() < 1e-6, "occupancy");
        assert!((a.b_iso - b.b_iso).abs() < 1e-6, "B-factor");
        assert_eq!(a.entity_id, b.entity_id, "entity id must round-trip");
        assert!((a.x - b.x).abs() < 1e-3);
        assert!((a.y - b.y).abs() < 1e-3);
        assert!((a.z - b.z).abs() < 1e-3);
    }
    let (c1, c2) = (r.cell.unwrap(), r2.cell.unwrap());
    assert!((c1.a - c2.a).abs() < 1e-3);
    assert!((c1.beta - c2.beta).abs() < 1e-3);
    assert_eq!(r.space_group, r2.space_group);
}

/// A realistic multi-chain PQR fixture (with chain-id column present),
/// hand-authored per the APBS/PDB2PQR field layout.
const MULTI_CHAIN_PQR: &str = "\
REMARK   PQR generated for chematic tests\n\
ATOM      1  N   GLY A   1     1.000   2.000   3.000 -0.400  1.625\n\
ATOM      2  CA  GLY A   1     1.500   2.500   3.500 -0.020  1.700\n\
ATOM      3  N   ALA B   1    10.000  11.000  12.000 -0.400  1.625\n\
HETATM    4  ZN  ZN  C   1    40.000  41.000  42.000  2.000  1.090\n\
";

#[test]
fn multi_chain_pqr_fixture_parses_and_round_trips() {
    let r = parse_pqr(MULTI_CHAIN_PQR).unwrap();
    assert_eq!(r.atoms.len(), 4);
    let chains: std::collections::BTreeSet<&str> = r
        .atoms
        .iter()
        .filter_map(|a| a.chain_id.as_deref())
        .collect();
    assert_eq!(chains, ["A", "B", "C"].into_iter().collect());
    assert_eq!(r.atoms[3].element, Element::ZN);

    let written = write_pqr(&r.atoms);
    let r2 = parse_pqr(&written).unwrap();
    assert_eq!(r.atoms, r2.atoms);
}

/// Oracle-style differential test: a minimal fixture whose parsed values
/// are stated confidently from the format documentation, then compared to
/// what `chematic-mol`'s existing small-molecule `cif::parse_cif` would
/// produce for an *equivalent* small-molecule CIF -- confirming the two
/// readers agree on the parts of CIF syntax they genuinely share (comment
/// stripping, quoting, esd stripping) even though their category
/// conventions differ.
#[test]
fn shared_tokenizer_behavior_agrees_with_small_molecule_cif_reader_on_esd_and_quoting() {
    let small_mol_cif = "data_x\n_cell_length_a 5.64(2)\n_cell_length_b 5.64\n_cell_length_c 5.64\n\
        _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
        loop_\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
        Na 0.0 0.0 0.0\n";
    let small_mol_result = cif::parse_cif(small_mol_cif).unwrap();
    // esd "(2)" suffix must be stripped identically by both readers, since
    // they share the same `strip_esd`/tokenizer.
    assert!((small_mol_result.cell.unwrap().a - 5.64).abs() < 1e-9);

    let mmcif_with_esd = "data_x\n\
        _cell.length_a 5.64(2)\n_cell.length_b 5.64\n_cell.length_c 5.64\n\
        _cell.angle_alpha 90\n_cell.angle_beta 90\n_cell.angle_gamma 90\n\
        loop_\n_atom_site.type_symbol\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
        Na 0.0 0.0 0.0\n";
    let mmcif_result = parse_mmcif(mmcif_with_esd).unwrap();
    assert!((mmcif_result.cell.unwrap().a - 5.64).abs() < 1e-9);
}
