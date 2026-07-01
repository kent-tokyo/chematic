use super::*;
use chematic_smiles::parse;

fn mol(s: &str) -> Molecule {
    parse(s).unwrap()
}

// --- Existing tests (must remain green) ---------------------------------

#[test]
fn test_alkanes() {
    assert_eq!(name(&mol("C")).unwrap(), "methane");
    assert_eq!(name(&mol("CC")).unwrap(), "ethane");
    assert_eq!(name(&mol("CCC")).unwrap(), "propane");
    assert_eq!(name(&mol("CCCC")).unwrap(), "butane");
    assert_eq!(name(&mol("CCCCC")).unwrap(), "pentane");
    assert_eq!(name(&mol("CCCCCC")).unwrap(), "hexane");
}

#[test]
fn test_alkenes_alkynes() {
    assert_eq!(name(&mol("C=C")).unwrap(), "ethene");
    assert_eq!(name(&mol("CC=C")).unwrap(), "propene");
    assert_eq!(name(&mol("C#C")).unwrap(), "ethyne");
    assert_eq!(name(&mol("CC#C")).unwrap(), "propyne");
}

#[test]
fn test_cycloalkanes() {
    assert_eq!(name(&mol("C1CC1")).unwrap(), "cyclopropane");
    assert_eq!(name(&mol("C1CCC1")).unwrap(), "cyclobutane");
    assert_eq!(name(&mol("C1CCCC1")).unwrap(), "cyclopentane");
    assert_eq!(name(&mol("C1CCCCC1")).unwrap(), "cyclohexane");
}

#[test]
fn test_alcohol() {
    assert_eq!(name(&mol("CO")).unwrap(), "methanol");
    assert_eq!(name(&mol("CCO")).unwrap(), "ethanol");
    assert_eq!(name(&mol("CCCO")).unwrap(), "propan-1-ol");
}

#[test]
fn test_amine() {
    assert_eq!(name(&mol("CN")).unwrap(), "methan-1-amine");
    assert_eq!(name(&mol("CCN")).unwrap(), "ethan-1-amine");
}

#[test]
fn test_haloalkane() {
    assert_eq!(name(&mol("CCCl")).unwrap(), "chloroethane");
    assert_eq!(name(&mol("CCBr")).unwrap(), "bromoethane");
    assert_eq!(name(&mol("CF")).unwrap(), "fluoromethane");
    assert_eq!(name(&mol("CI")).unwrap(), "iodomethane");
}

#[test]
fn test_not_supported() {
    assert!(name(&mol("CC.CC")).is_err()); // disconnected
}

#[test]
fn test_empty() {
    use chematic_core::MoleculeBuilder;
    let mol = MoleculeBuilder::new().build();
    assert_eq!(name(&mol), Err(IupacError::Empty));
}

// --- New: benzene & aromatic heterocycles --------------------------------

#[test]
fn test_benzene() {
    assert_eq!(name(&mol("c1ccccc1")).unwrap(), "benzene");
}

#[test]
fn test_aromatic_heterocycles() {
    assert_eq!(name(&mol("c1ccncc1")).unwrap(), "pyridine");
    assert_eq!(name(&mol("c1ccoc1")).unwrap(), "furan");
    assert_eq!(name(&mol("c1ccsc1")).unwrap(), "thiophene");
    assert_eq!(name(&mol("c1cc[nH]c1")).unwrap(), "pyrrole");
    assert_eq!(name(&mol("c1cnc[nH]1")).unwrap(), "imidazole");
}

// --- New: ketones with position locant -----------------------------------

#[test]
fn test_ketones() {
    assert_eq!(name(&mol("CC(=O)C")).unwrap(), "propan-2-one");
    assert_eq!(name(&mol("CC(=O)CC")).unwrap(), "butan-2-one");
    assert_eq!(name(&mol("CCC(=O)CC")).unwrap(), "pentan-3-one");
    assert_eq!(name(&mol("CCCC(=O)C")).unwrap(), "pentan-2-one");
}

// --- New: carboxylic acids -----------------------------------------------

#[test]
fn test_carboxylic_acids() {
    assert_eq!(name(&mol("CC(=O)O")).unwrap(), "ethanoic acid");
    assert_eq!(name(&mol("CCC(=O)O")).unwrap(), "propanoic acid");
    assert_eq!(name(&mol("C(=O)O")).unwrap(), "methanoic acid");
}

// --- New: esters ---------------------------------------------------------

#[test]
fn test_esters() {
    assert_eq!(name(&mol("CC(=O)OC")).unwrap(), "methyl ethanoate");
    assert_eq!(name(&mol("C(=O)OC")).unwrap(), "methyl methanoate");
    assert_eq!(name(&mol("CC(=O)OCC")).unwrap(), "ethyl ethanoate");
}

// --- New: amides ---------------------------------------------------------

#[test]
fn test_amides() {
    assert_eq!(name(&mol("CC(=O)N")).unwrap(), "ethanamide");
    assert_eq!(name(&mol("C(=O)N")).unwrap(), "methanamide");
    assert_eq!(name(&mol("CCC(=O)N")).unwrap(), "propanamide");
}

// ---- New: branched alkanes (v0.1.101) ------------------------------------

#[test]
fn test_branched_alkanes() {
    assert_eq!(name(&mol("CC(C)C")).unwrap(), "2-methylpropane");
    assert_eq!(name(&mol("CC(C)CC")).unwrap(), "2-methylbutane");
    assert_eq!(name(&mol("CC(C)(C)C")).unwrap(), "2,2-dimethylpropane");
    assert_eq!(name(&mol("CCCC(C)CC")).unwrap(), "3-methylhexane");
}

#[test]
fn test_branched_alkane_lowest_locant() {
    // CCC(C)C = 2-methylbutane (not 3-methylbutane — lower locant wins).
    assert_eq!(name(&mol("CCC(C)C")).unwrap(), "2-methylbutane");
}

// ---- New: substituted benzenes (v0.1.101) --------------------------------

#[test]
fn test_substituted_benzenes() {
    assert_eq!(name(&mol("c1ccccc1O")).unwrap(), "phenol");
    assert_eq!(name(&mol("c1ccccc1N")).unwrap(), "aniline");
    assert_eq!(name(&mol("c1ccccc1Cl")).unwrap(), "chlorobenzene");
    assert_eq!(name(&mol("c1ccccc1Br")).unwrap(), "bromobenzene");
}

#[test]
fn test_substituted_benzene_carbonyl() {
    assert_eq!(name(&mol("c1ccccc1C=O")).unwrap(), "benzaldehyde");
    assert_eq!(name(&mol("c1ccccc1C(=O)O")).unwrap(), "benzoic acid");
}

// ---- New: nitriles (v0.1.101) -------------------------------------------

#[test]
fn test_nitriles() {
    assert_eq!(name(&mol("CC#N")).unwrap(), "ethanenitrile");
    assert_eq!(name(&mol("CCC#N")).unwrap(), "propanenitrile");
}

// ---- New Round 2 tests (v0.1.102) ---------------------------------------

#[test]
fn test_thiols() {
    assert_eq!(name(&mol("CS")).unwrap(), "methanethiol");
    assert_eq!(name(&mol("CCS")).unwrap(), "ethanethiol");
    assert_eq!(name(&mol("CCCS")).unwrap(), "propanethiol");
}

#[test]
fn test_alcohol_locants() {
    assert_eq!(name(&mol("CCCCO")).unwrap(), "butan-1-ol");
    assert_eq!(name(&mol("CC(O)C")).unwrap(), "propan-2-ol");
    assert_eq!(name(&mol("CCC(O)C")).unwrap(), "butan-2-ol");
}

#[test]
fn test_disubstituted_benzene() {
    // Para-chlorophenol: OH and Cl are 3 bonds apart in the ring (positions 1 and 4).
    assert_eq!(name(&mol("Oc1ccc(Cl)cc1")).unwrap(), "4-chlorophenol");
    // Meta-chlorophenol: OH and Cl are 2 bonds apart (positions 1 and 3).
    assert_eq!(name(&mol("c1ccc(O)cc1Cl")).unwrap(), "3-chlorophenol");
}

#[test]
fn test_methylcycloalkane() {
    assert_eq!(name(&mol("CC1CCCCC1")).unwrap(), "methylcyclohexane");
    assert_eq!(name(&mol("CC1CCCC1")).unwrap(), "methylcyclopentane");
    assert_eq!(name(&mol("CC1CCC1")).unwrap(), "methylcyclobutane");
}

// ---- New Round 3 tests (v0.1.103) ----------------------------------------

#[test]
fn test_ethers() {
    assert_eq!(name(&mol("COC")).unwrap(), "methoxymethane");
    assert_eq!(name(&mol("COCC")).unwrap(), "methoxyethane");
    assert_eq!(name(&mol("CCOCC")).unwrap(), "ethoxyethane");
    assert_eq!(name(&mol("COCCC")).unwrap(), "1-methoxypropane");
}

#[test]
fn test_trimethylbenzene() {
    assert_eq!(
        name(&mol("Cc1cccc(C)c1C")).unwrap(),
        "1,2,3-trimethylbenzene"
    );
    assert_eq!(
        name(&mol("Cc1ccc(C)cc1C")).unwrap(),
        "1,2,4-trimethylbenzene"
    );
    assert_eq!(
        name(&mol("Cc1cc(C)cc(C)c1")).unwrap(),
        "1,3,5-trimethylbenzene"
    );
}

#[test]
fn test_secondary_amine() {
    assert_eq!(name(&mol("CCNCC")).unwrap(), "N-ethylethanamine");
    assert_eq!(name(&mol("CNCC")).unwrap(), "N-methylethanamine");
    assert_eq!(name(&mol("CN(C)C")).unwrap(), "N,N-dimethylmethanamine");
}

// ---- New Round 9 tests (v0.1.109) ----------------------------------------

#[test]
fn test_branched_aldehyde() {
    assert_eq!(name(&mol("CC(C)C=O")).unwrap(), "2-methylpropanal");
    assert_eq!(name(&mol("CCC(C)C=O")).unwrap(), "2-methylbutanal");
}

#[test]
fn test_branched_amide() {
    assert_eq!(name(&mol("CC(C)C(=O)N")).unwrap(), "2-methylpropanamide");
    assert_eq!(name(&mol("CCC(C)C(=O)N")).unwrap(), "2-methylbutanamide");
}

// ---- New Round 8 tests (v0.1.108) ----------------------------------------

#[test]
fn test_branched_ester() {
    assert_eq!(
        name(&mol("CC(C)C(=O)OC")).unwrap(),
        "methyl 2-methylpropanoate"
    );
    assert_eq!(
        name(&mol("CC(C)C(=O)OCC")).unwrap(),
        "ethyl 2-methylpropanoate"
    );
}

#[test]
fn test_branched_ketone() {
    assert_eq!(name(&mol("CC(=O)C(C)C")).unwrap(), "3-methylbutan-2-one");
    assert_eq!(
        name(&mol("CC(=O)C(C)(C)C")).unwrap(),
        "3,3-dimethylbutan-2-one"
    );
}

// ---- New Round 7 tests (v0.1.107) ----------------------------------------

#[test]
fn test_secondary_thiol() {
    assert_eq!(name(&mol("CCC(S)C")).unwrap(), "butane-2-thiol");
    assert_eq!(name(&mol("CCCC(S)C")).unwrap(), "pentane-2-thiol");
}

#[test]
fn test_branched_carboxylic_acid() {
    assert_eq!(name(&mol("CC(C)C(=O)O")).unwrap(), "2-methylpropanoic acid");
    assert_eq!(name(&mol("CCC(C)C(=O)O")).unwrap(), "2-methylbutanoic acid");
    assert_eq!(
        name(&mol("CC(C)(C)C(=O)O")).unwrap(),
        "2,2-dimethylpropanoic acid"
    );
}

// ---- New Round 6 tests (v0.1.106) ----------------------------------------

#[test]
fn test_alkene_locants() {
    assert_eq!(name(&mol("CC=CC")).unwrap(), "but-2-ene");
    assert_eq!(name(&mol("C=CCC")).unwrap(), "but-1-ene");
    assert_eq!(name(&mol("CC=CCC")).unwrap(), "pent-2-ene");
    assert_eq!(name(&mol("C=CCCC")).unwrap(), "pent-1-ene");
}

#[test]
fn test_alkyne_locants() {
    assert_eq!(name(&mol("CC#CC")).unwrap(), "but-2-yne");
    assert_eq!(name(&mol("C#CCC")).unwrap(), "but-1-yne");
}

#[test]
fn test_amine_locants() {
    assert_eq!(name(&mol("CCCN")).unwrap(), "propan-1-amine");
    assert_eq!(name(&mol("CCC(N)C")).unwrap(), "butan-2-amine");
    assert_eq!(name(&mol("CC(N)CCC")).unwrap(), "pentan-2-amine");
}

// ---- New Round 5 tests (v0.1.105) ----------------------------------------

#[test]
fn test_haloalkane_locants() {
    // n=3: terminal → "1-chloropropane"
    assert_eq!(name(&mol("CCCCl")).unwrap(), "1-chloropropane");
    // n=4: terminal → "1-chlorobutane"
    assert_eq!(name(&mol("CCCCCl")).unwrap(), "1-chlorobutane");
    // n=4: internal → "2-chlorobutane"
    assert_eq!(name(&mol("CCC(Cl)C")).unwrap(), "2-chlorobutane");
    // n=5: internal → "2-chloropentane"
    assert_eq!(name(&mol("CCCC(Cl)C")).unwrap(), "2-chloropentane");
    // di-halo: ClCCCl = 2C → "1,2-dichloroethane"; ClCCCCl = 3C → "1,3-dichloropropane"
    assert_eq!(name(&mol("ClCCCl")).unwrap(), "1,2-dichloroethane");
    assert_eq!(name(&mol("ClCCCCl")).unwrap(), "1,3-dichloropropane");
}

#[test]
fn test_cycloalkanol() {
    assert_eq!(name(&mol("OC1CCC1")).unwrap(), "cyclobutanol");
    assert_eq!(name(&mol("OC1CCCC1")).unwrap(), "cyclopentanol");
    assert_eq!(name(&mol("OC1CCCCC1")).unwrap(), "cyclohexanol");
}

// ---- New Round 4 tests (v0.1.104) ----------------------------------------

#[test]
fn test_disubstituted_benzene_non_principal() {
    // Two halogens (para)
    assert_eq!(
        name(&mol("Clc1ccc(Br)cc1")).unwrap(),
        "1-bromo-4-chlorobenzene"
    );
    assert_eq!(
        name(&mol("Clc1ccc(F)cc1")).unwrap(),
        "1-chloro-4-fluorobenzene"
    );
    // Two methyls: ortho (Cc1ccccc1C) and para (Cc1ccc(C)cc1)
    assert_eq!(name(&mol("Cc1ccccc1C")).unwrap(), "1,2-dimethylbenzene");
    assert_eq!(name(&mol("Cc1ccc(C)cc1")).unwrap(), "1,4-dimethylbenzene");
    // Methyl + halogen (para)
    assert_eq!(
        name(&mol("Cc1ccc(Cl)cc1")).unwrap(),
        "1-chloro-4-methylbenzene"
    );
}

#[test]
fn test_propyl_substituent() {
    // 11C: longest chain = octane (8C), propyl substituent at C4
    assert_eq!(name(&mol("CCCC(CCC)CCCC")).unwrap(), "4-propyloctane");
}

#[test]
fn test_dimethylcycloalkane() {
    assert_eq!(
        name(&mol("CC1CCC(C)CC1")).unwrap(),
        "1,4-dimethylcyclohexane"
    );
    assert_eq!(name(&mol("CC1CCCC1C")).unwrap(), "1,2-dimethylcyclopentane");
    assert_eq!(
        name(&mol("CC1CCC(C)C1")).unwrap(),
        "1,3-dimethylcyclopentane"
    );
}

// ---- Sprint 4 tests: N-heterocycles, sulfides, naphthalene ------------------

#[test]
fn test_saturated_n_rings() {
    assert_eq!(name(&mol("C1CCNCC1")).unwrap(), "piperidine");
    assert_eq!(name(&mol("C1CCNC1")).unwrap(), "pyrrolidine");
    assert_eq!(name(&mol("C1CCN1")).unwrap(), "azetidine");
}

#[test]
fn test_morpholine() {
    // Morpholine: 6-membered ring with O and N
    assert_eq!(name(&mol("C1COCCN1")).unwrap(), "morpholine");
    assert_eq!(name(&mol("C1CNCCO1")).unwrap(), "morpholine");
}

#[test]
fn test_piperazine() {
    // Piperazine: 4C + 2N in 6-membered ring
    assert_eq!(name(&mol("C1CNCCN1")).unwrap(), "piperazine");
}

#[test]
fn test_naphthalene() {
    assert_eq!(name(&mol("c1ccc2ccccc2c1")).unwrap(), "naphthalene");
}

#[test]
fn test_sulfide_naming() {
    assert_eq!(name(&mol("CSC")).unwrap(), "methyl methyl sulfide");
    assert_eq!(name(&mol("CSCC")).unwrap(), "ethyl methyl sulfide");
    assert_eq!(name(&mol("CCSCC")).unwrap(), "ethyl ethyl sulfide");
}

// ---- C2: spiro / bicyclo naming (v0.4.x) --------------------------------

#[test]
fn test_spiro_naming() {
    // spiro[4.4]nonane — two cyclopentane rings sharing one atom (9 C total)
    // SMILES: C1CCC2(C1)CCCC2 = 9 atoms
    let m = mol("C1CCC2(C1)CCCC2");
    assert_eq!(name(&m).unwrap(), "spiro[4.4]nonane");
    // spiro[4.5]decane — cyclopentane + cyclohexane (10 C total)
    let m3 = mol("C1CCC2(C1)CCCCC2");
    assert_eq!(name(&m3).unwrap(), "spiro[4.5]decane");
    // spiro[2.3]hexane — cyclopropane + cyclobutane (6 C total)
    // Ring1 = 4-membered (bridge=3) via C1..C1, Ring2 = 3-membered (bridge=2) via C2..C2
    let m4 = mol("C1CC2(C1)CC2");
    assert_eq!(name(&m4).unwrap(), "spiro[2.3]hexane");
}

#[test]
fn test_bicyclo_naming() {
    // bicyclo[2.2.2]octane (8 C: 2 bridgeheads + bridges 2,2,2)
    // SMILES: C1CC2CCC1CC2 = 8 atoms, bridges 2,2,2
    let m = mol("C1CC2CCC1CC2");
    assert_eq!(name(&m).unwrap(), "bicyclo[2.2.2]octane");
    // bicyclo[2.2.1]heptane = norbornane (7 C: bridges 2,2,1)
    // SMILES: C1CC2CCC1C2 = 7 atoms
    let m2 = mol("C1CC2CCC1C2");
    assert_eq!(name(&m2).unwrap(), "bicyclo[2.2.1]heptane");
}
