//! MMFF94 atom types RDKit 2026.03.6 gives the residual rows of #734/#754
//! batch 10 (exposed 10k and ChEMBL 5k), pinned from
//! `MMFFGetMoleculeProperties` on `AddHs(MolFromSmiles(smiles))`; heavy
//! atoms only (the hydrogens follow them).

fn heavy_types(smiles: &str) -> Vec<u8> {
    let mol = chematic_smiles::parse(smiles).expect("parses");
    let heavy = mol.atom_count();
    let h = chematic_chem::add_hydrogens(&mol);
    let types = chematic_ff::mmff94_numeric::assign_mmff94_numeric_types(&h).expect("types");
    types[..heavy].to_vec()
}

#[test]
fn fullerene_bis_adduct_cage_kekule() {
    let smiles = "COCCOCCOCCN1CC23C4=C5C6=C7c8c9c%10c%11c%12c%13c%14c(c2c2c%15c%16c%17c%18c%19c(c5c5c%20c%21c%22c%23c(c8C%22C65)c%10c5c%11c6c%13c8c(c%15%14)c%16c%10c%18c%11c(c%20%19)c%21c%13c%23c5c5c%13c%11c%10c8c65)C%17C42)C%12C9C73C1COCCOCCOC";
    let rdkit: [u8; 82] = [
        1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 8, 1, 1, 2, 2, 2, 2, 37, 37, 37, 2, 2, 37, 37, 37, 37, 37,
        37, 37, 37, 37, 2, 2, 2, 2, 2, 2, 2, 2, 37, 37, 1, 1, 37, 2, 2, 2, 37, 37, 37, 37, 37, 37,
        37, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 37, 37, 2, 1, 1, 1, 1, 1, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1,
    ];
    assert_eq!(heavy_types(smiles), rdkit);
}

#[test]
fn neutral_n_next_to_imine_and_triaminated_carbon_is_ngd_plus() {
    let smiles = "C1=CC=C(C=C1)N2N=C3N(C2N3C4=CC=CC=C4)C5=CC=CC=C5";
    let rdkit: [u8; 24] = [
        37, 37, 37, 37, 37, 37, 40, 9, 3, 56, 20, 40, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37,
        37,
    ];
    assert_eq!(heavy_types(smiles), rdkit);
}

#[test]
fn porphyrin_imine_n_in_accepted_macrocycle_is_nm() {
    let smiles =
        "CC1=C2NC(=C1CCC(O)=O)C=C3N=C(C=C4NC(=CC5=NC(=C2)C(=C5C)C=C)C(=C4C)C=C)C(=C3CCC(O)=O)C";
    let rdkit: [u8; 42] = [
        1, 64, 63, 39, 63, 64, 1, 1, 3, 6, 7, 2, 2, 62, 2, 2, 63, 39, 63, 2, 2, 62, 2, 2, 2, 2, 1,
        2, 2, 64, 64, 1, 2, 2, 2, 2, 1, 1, 3, 6, 7, 1,
    ];
    assert_eq!(heavy_types(smiles), rdkit);
}

#[test]
fn sulfonyl_azide_n_is_imine_azo() {
    let smiles = "CC1(C)O[C@H]2[C@@H]3OS(=O)(=O)O[C@@H]3CO[C@@]2(COS(=O)(=O)N=[N+]=[N-])O1";
    let rdkit: [u8; 24] = [
        1, 1, 1, 6, 1, 1, 6, 18, 32, 32, 6, 1, 1, 6, 1, 1, 6, 18, 32, 32, 9, 53, 47, 6,
    ];
    assert_eq!(heavy_types(smiles), rdkit);
}

#[test]
fn sulfonylimine_n_is_nso() {
    let smiles =
        "Cc1ccc(S(=O)(=O)N=S2(=O)O[C@@H]3[C@@H](CO[C@@]4(COS(N)(=O)=O)OC(C)(C)O[C@@H]34)O2)cc1";
    let rdkit: [u8; 32] = [
        1, 37, 37, 37, 37, 18, 32, 32, 48, 18, 32, 6, 1, 1, 1, 6, 1, 1, 6, 18, 43, 32, 32, 6, 1, 1,
        1, 6, 1, 6, 37, 37,
    ];
    assert_eq!(heavy_types(smiles), rdkit);
}

#[test]
fn imidazolium_c2_is_cim_plus() {
    let smiles = "CCc1[nH]c(O)nc1C(=O)c1ccc(-n2cc[n+](C)c2C)cc1";
    let rdkit: [u8; 23] = [
        1, 1, 63, 39, 63, 6, 66, 64, 3, 7, 37, 37, 37, 37, 81, 78, 78, 81, 1, 80, 1, 37, 37,
    ];
    assert_eq!(heavy_types(smiles), rdkit);
}

#[test]
fn cationic_tricyclic_amidinium_n_are_ncn_plus() {
    let smiles = "Cc1cn2c([O-])c3c(nc2n1)[n+](COCCO)cn3C";
    let rdkit: [u8; 20] = [
        1, 64, 63, 39, 2, 35, 2, 3, 9, 63, 66, 55, 1, 6, 1, 1, 6, 57, 55, 1,
    ];
    assert_eq!(heavy_types(smiles), rdkit);
}

#[test]
fn kekule_written_bis_quinolinium_uses_rdkit_kekule_structure() {
    // Kekule-written copy (RDKit `MolToSmiles(kekuleSmiles=True)`): RDKit
    // perceives aromaticity again and types from its canonical Kekule form.
    let smiles = "C1=CC=C2C(=C1)C1=CC=[N+]2CC2=CC(=CC=C2)C2=CC=CC(=C2)C[N+]2=CC=C(NCCCCCCCCCCN1)C1=CC=CC=C12";
    let rdkit: [u8; 46] = [
        37, 37, 37, 37, 37, 37, 37, 37, 37, 58, 1, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37,
        1, 54, 3, 2, 2, 40, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 40, 37, 37, 37, 37, 37, 37,
    ];
    assert_eq!(heavy_types(smiles), rdkit);
}

#[test]
fn kekule_written_fullerene_adduct_uses_rdkit_symmetrized_rings() {
    // Kekule-written copy (RDKit `MolToSmiles(kekuleSmiles=True)`): RDKit
    // perceives aromaticity again and types from its canonical Kekule form.
    let smiles = "COCCOCCOCCN1CC23C4=C5C6=C7C8=C9C%10=C%11C%12=C%13C%14=C%15C%16=C%17C%14=C%14C%12=C%12C%10=C%10C%18=C%19C%12=C%14C%12=C%17C%14C%17=C%12C%19=C%12C%18=C(C6=C9%10)C5C5C%12=C%17C(=C6C9=C2C2C4=C7C4=C8C%11=C%13C7=C4C2C(=C%157)C9=C%16C6%14)C53C1COCCOCCOC";
    let rdkit: [u8; 82] = [
        1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 8, 1, 1, 37, 37, 37, 37, 37, 37, 37, 37, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 37, 37, 37, 37, 2, 2, 37, 37, 1, 1, 37, 37, 2, 2, 2, 2,
        1, 37, 37, 37, 37, 37, 37, 37, 37, 1, 2, 2, 2, 2, 1, 1, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1,
    ];
    assert_eq!(heavy_types(smiles), rdkit);
}
