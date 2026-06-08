//! ChEMBL-style SMILES roundtrip integration test.
//!
//! For each SMILES:
//!   1. Parse → mol1
//!   2. canonical_smiles(mol1) → smiles1
//!   3. Parse(smiles1) → mol2
//!   4. canonical_smiles(mol2) → smiles2
//!   5. Assert atom/bond counts match and smiles1 == smiles2 (canonical stability).

use chematic_smiles::{canonical_smiles, parse};

/// 50 representative SMILES covering a broad range of chemical space.
const CHEMBL_SMILES: &[&str] = &[
    // --- Simple aromatics ---
    "c1ccccc1",         // benzene
    "c1ccc2ccccc2c1",   // naphthalene
    "c1ccc2[nH]ccc2c1", // indole
    // --- Heteroaromatics ---
    "c1ccncc1",   // pyridine
    "c1cn[nH]c1", // pyrazole
    "c1cc[nH]c1", // pyrrole
    "c1ccoc1",    // furan
    "c1ccsc1",    // thiophene
    "c1cnc[nH]1", // imidazole
    "c1ccno1",    // isoxazole
    // --- Stereocenters ---
    "[C@@H](F)(Cl)Br",                        // chiral carbon (CCW)
    "[C@H](F)(Cl)Br",                         // chiral carbon (CW)
    "N[C@@H](C)C(=O)O",                       // L-alanine-like
    "N[C@H](C)C(=O)O",                        // D-alanine-like
    "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", // glucose
    // --- E/Z double bonds ---
    "F/C=C/F",    // trans-difluoroethylene
    "F/C=C\\F",   // cis-difluoroethylene
    "C(/F)=C/Cl", // E-fluorochloroethylene
    // --- Multi-fragment (salts/ions) ---
    "[Na+].[Cl-]",            // sodium chloride
    "[NH4+].[Cl-]",           // ammonium chloride
    "[K+].[O-]S(=O)(=O)[O-]", // potassium sulfate fragment
    // --- Common drugs ---
    "CC(=O)Oc1ccccc1C(=O)O",               // aspirin
    "Cn1cnc2c1c(=O)n(c(=O)n2C)C",          // caffeine
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",          // ibuprofen
    "CC(=O)Nc1ccc(O)cc1",                  // paracetamol (acetaminophen)
    "OC(=O)c1ccccc1O",                     // salicylic acid
    "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C", // testosterone-like steroid
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",        // theophylline
    "c1ccc(cc1)C(=O)O",                    // benzoic acid
    // --- Bracket atoms with charges ---
    "[NH4+]", // ammonium ion
    "[O-]",   // oxide ion
    "[OH3+]", // hydronium
    "[Fe+2]", // ferrous iron
    "[Ca+2]", // calcium ion
    "[Mg+2]", // magnesium ion
    // --- Fused ring systems ---
    "c1ccc2cc3ccccc3cc2c1",        // anthracene
    "c1ccc2ccc3cccc4ccc(c1)c2c34", // pyrene
    "C1CC2CCCC2CC1",               // decalin
    "O=C1NC(=O)c2ccccc21",         // isatoic anhydride-like
    // --- Aliphatic chains ---
    "CCCCCCCC",  // octane
    "CC(C)CCC",  // isohexane
    "CC(=O)O",   // acetic acid
    "CCCC(=O)O", // butanoic acid
    "OCC(O)CO",  // glycerol-like
    // --- Nitrogen/oxygen heterocycles ---
    "C1CCNCC1",     // piperidine
    "C1CCOCC1",     // tetrahydropyran (THP)
    "C1CCOC1",      // tetrahydrofuran (THF)
    "C1CCNC1",      // pyrrolidine
    "O=C1CCCN1",    // 2-pyrrolidinone (gamma-lactam)
    "C1CN2CCCC2C1", // octahydroindolizine fragment
    "n1ccnc1",      // imidazole-ish tautomer
];

#[test]
fn roundtrip_consistency() {
    let mut failures: Vec<String> = Vec::new();

    for &s in CHEMBL_SMILES {
        // Step 1: parse original
        let mol1 = match parse(s) {
            Ok(m) => m,
            Err(e) => {
                failures.push(format!("PARSE FAILED '{}': {}", s, e));
                continue;
            }
        };

        // Step 2: canonical SMILES of mol1
        let smiles1 = canonical_smiles(&mol1);

        // Step 3: re-parse canonical
        let mol2 = match parse(&smiles1) {
            Ok(m) => m,
            Err(e) => {
                failures.push(format!(
                    "RE-PARSE FAILED '{}' (canonical='{}'):  {}",
                    s, smiles1, e
                ));
                continue;
            }
        };

        // Step 4: canonical SMILES of mol2 (must be identical to smiles1)
        let smiles2 = canonical_smiles(&mol2);

        // Assert atom count preserved
        if mol1.atom_count() != mol2.atom_count() {
            failures.push(format!(
                "ATOM COUNT MISMATCH '{}': {} vs {} (canonical='{}')",
                s,
                mol1.atom_count(),
                mol2.atom_count(),
                smiles1
            ));
        }

        // Assert bond count preserved
        if mol1.bond_count() != mol2.bond_count() {
            failures.push(format!(
                "BOND COUNT MISMATCH '{}': {} vs {} (canonical='{}')",
                s,
                mol1.bond_count(),
                mol2.bond_count(),
                smiles1
            ));
        }

        // Assert canonical stability
        if smiles1 != smiles2 {
            failures.push(format!(
                "CANONICAL UNSTABLE '{}': '{}' -> '{}'",
                s, smiles1, smiles2
            ));
        }
    }

    if !failures.is_empty() {
        panic!(
            "{} roundtrip failure(s):\n{}",
            failures.len(),
            failures.join("\n")
        );
    }
}
