//! RDKit 2026.09.1 bit-vector oracle for a ring triple in an aromatic graph.

use super::*;

#[test]
fn rdkit_accepted_aromatic_triple_matches_all_five_fingerprints() {
    let mol =
        chematic_smiles::parse("Cc1cccc(Cn2:C(c3ccc(O)cc#3):c(C)c3cc(O)cc[c:42]23)c1").unwrap();
    for (name, actual, expected) in [
        (
            "Morgan",
            rdkit_morgan_ecfp4_bitvec(&mol).unwrap(),
            0x9aaa205f415c1b55,
        ),
        ("Pattern", rdkit_pattern_fp(&mol), 0x9f2729a74d533e1c),
        ("RDK", rdkit_rdk_fp(&mol), 0x37c52445ea97d54f),
        ("Torsion", rdkit_torsion_fp(&mol), 0x317dedccfa08cee5),
        ("Atom Pair", rdkit_atom_pair_fp(&mol), 0xd9ca2b8f33db5a37),
    ] {
        assert_eq!(fingerprint_test_digest(&actual), expected, "{name}");
    }
}
