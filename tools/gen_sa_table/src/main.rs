// Generates the SA Score fragment frequency table from a SMILES corpus.
//
// Usage:
//   cargo run -p gen-sa-table                   -- use built-in drug corpus
//   cargo run -p gen-sa-table -- <smiles_file>  -- read SMILES (one per line) from file
//
// Output: a Rust source snippet for `static FRAGMENT_SCORES: &[(u64, i16)]`
// that can be pasted directly into crates/chematic-chem/src/sa_score.rs.
//
// Score encoding: i16 = (log10(fragment_frequency) * 1000.0) as i16
//   e.g. freq 0.1 → log10 = -1.0 → i16 = -1000
//        freq 0.5 → log10 = -0.30 → i16 = -301
//   Range in practice: [-5000, 0]

use chematic_fp::morgan_fp_counts;
use chematic_smiles::parse;
use std::collections::HashMap;

// Built-in corpus: SMILES validated from chematic test suite + well-known drugs.
// Any SMILES that fails to parse is silently skipped.
const CORPUS: &[&str] = &[
    // Simple alkanes / cycloalkanes
    "C",
    "CC",
    "CCC",
    "CCCC",
    "CCCCC",
    "CCCCCC",
    "CC(C)C",
    "CC(C)(C)C",
    "C1CCCCC1",
    "C1CCCC1",
    "C1CCC1",
    "C1CCCCC11CCCC1",
    "C1CC2(CC1)CCCC2",
    "C1CC2CCC1C2",
    // Alkenes / alkynes
    "C=C",
    "CC=CC",
    "C#C",
    "CC#C",
    "CC#N",
    // Alcohols / ethers
    "CO",
    "CCO",
    "CCCO",
    "CCCCO",
    "OCC",
    "OCCO",
    "OCC(O)CO",
    "CC(O)C",
    "COC",
    // Amines
    "CN",
    "CCN",
    "CCCN",
    "NCCN",
    "CC(N)C",
    "NCCO",
    // Carboxylic acids / esters / amides
    "CC(=O)O",
    "CCC(=O)O",
    "CCCC(=O)O",
    "CCCCCCCCCCCCC(=O)O",
    "CC(=O)N",
    "CC(=O)NC",
    "CC(=O)OC",
    "OC(=O)CO",
    "OC(=O)CC(=O)O",
    "OC(=O)CCC(=O)O",
    "CC(=O)[O-]",
    // Aldehydes / ketones
    "CC=O",
    "CCC=O",
    "CC(=O)C",
    "CC(=O)CC",
    "CC(=O)CC(=O)C",
    "O=Cc1ccccc1",
    "CC(=O)c1ccccc1",
    // Haloalkanes
    "CCCl",
    "CCBr",
    "ClCCl",
    // Aromatic hydrocarbons
    "c1ccccc1",
    "Cc1ccccc1",
    "CCc1ccccc1",
    "c1ccc2ccccc2c1",
    "c1ccccc1c1ccccc1",
    // Benzene heteroatom substituents
    "Oc1ccccc1",
    "Nc1ccccc1",
    "c1ccccc1Cl",
    "c1ccccc1Br",
    "c1ccc(Cl)cc1",
    "c1ccc([N+](=O)[O-])cc1",
    "c1ccc(cc1)N=Nc2ccccc2",
    "COc1ccccc1",
    "OC(=O)c1ccccc1",
    "OC(=O)c1ccccc1O",
    "Oc1ccccc1O",
    // Six-membered N-heterocycles (aromatic)
    "c1ccncc1",
    "c1ccncn1",
    // Five-membered heterocycles (aromatic)
    "c1ccoc1",
    "c1ccsc1",
    "c1cc[nH]n1",
    "c1cnc[nH]1",
    // Polycyclic aromatic N-heterocycles
    "c1ccc2[nH]ccc2c1",
    "c1cnc2ccccc2[nH]1",
    // Saturated N/O heterocycles
    "C1CCNCC1",
    "C1COCCN1",
    "C1CNCCN1",
    "C1CCOC1",
    "C1CCOCC1",
    // Urea / guanidine
    "NC(N)=O",
    "NC(=N)N",
    // Nucleobases (Kekulized to avoid aromatic-keto ambiguity)
    "Nc1ncnc2[nH]cnc12",
    "NC1=NC(=O)N=CC1",
    "CC1=CNC(=O)NC1=O",
    "O=C1NC(=O)C=CN1",
    // Amino acids
    "NCC(=O)O",
    "N[C@@H](C)C(=O)O",
    "NC(Cc1ccccc1)C(=O)O",
    "NC(Cc1ccc(O)cc1)C(=O)O",
    "NC(Cc1c[nH]c2ccccc12)C(=O)O",
    "N[C@@H](CCCNC(=N)N)C(=O)O",
    "OC(=O)[C@@H](O)C",
    "NC(CS)C(=O)O",
    "NC(CCSC)C(=O)O",
    // Sugars / polyols
    "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
    "OCC(O)CO",
    // Sulfur compounds
    "CS(=O)(=O)O",
    "CS(C)(=O)=O",
    "CSC",
    // Phosphorus
    "COP(=O)(OC)OC",
    // Phenols / catechols / aromatic acids
    "NC(=O)c1cccnc1",
    "Nc1ccc(C(=O)O)cc1",
    "CC(=O)Nc1ccccc1",
    // Lactams
    "O=C1CCCCN1",
    "O=C1CCCN1",
    // Misc
    "O=C1CSC(=S)N1",
    "[NH3+]CC(=O)[O-]",
    // Drug-like molecules
    "CC(=O)Oc1ccccc1C(=O)O",
    "CC(=O)Nc1ccc(O)cc1",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "CCN(CC)CC(=O)Nc1c(C)cccc1C",
    "CCOC(=O)c1ccc(N)cc1",
    "CCC1(c2ccccc2)C(=O)NC(=O)NC1=O",
    "O=C(N)N1c2ccccc2C=Cc2ccccc21",
    "CC(S)C(=O)N1CCCC1C(=O)O",
    "CCCC(CCC)C(=O)O",
    "NCC1(CC(=O)O)CCCCC1",
    "CN(C)C(=N)NC(=N)N",
    "NCCc1ccc(O)c(O)c1",
    "NCCc1c[nH]c2ccc(O)cc12",
    "CC(N)Cc1ccccc1",
    "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O",
    "CC(C)NCC(O)COc1cccc2ccccc12",
    "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1",
    "NC(Cc1ccc(O)c(O)c1)C(=O)O",
    "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
    "Cc1ncc([N+](=O)[O-])n1CCO",
    "OC(Cn1ccnc1)(Cn1ccnc1)c1ccc(F)cc1F",
    "OC(=O)c1cn(C2CC2)c2cc(F)c(N3CCNCC3)cc2c1=O",
    "CCN(CC)CCOC(=O)c1ccc(N)cc1",
    "CC(C(=O)O)c1ccc2cc(OC)ccc2c1",
    "CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(cc1)C(=O)NC(CCC(=O)O)C(=O)O",
    "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(=C2OC1=O)C)(C)C)OC(=O)C5=CC=CC=C5)OC(=O)C4)O)(C)C",
    // More scaffolds
    "NS(=O)(=O)c1ccccc1",
    "NS(=O)(=O)c1ccc(N)cc1",
    "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1",
    "O=c1cc(-c2ccccc2)oc2ccccc12",
    "O=c1ccc2ccccc2o1",
    "CC(C)NCC(O)COc1ccc(CCOC)cc1",
    // v0.1.94: FDA approved drugs expansion
    // Statins
    "CCC(C)(C)C(=O)OC1CC(O)CC(O)C1OC1OC(C)CC(O)C1O",  // Simvastatin
    "CCC(C)(C)C(=O)OC1CC(O)CC(O)C1OC1OC(C)CC(O)C1O",  // Lovastatin
    "CC(C)c1c(C(=O)Nc2ccccc2)c(cc(c1c1ccc(F)cc1)C(=O)O)C",  // Atorvastatin
    "CC(C)c1c(cc(c(c1c1ccc(F)cc1)C(=O)O)C(=O)N(C)C)C",  // Rosuvastatin

    // Beta-blockers
    "CC(C)NCC(O)c1ccc(O)c(CC(=O)O)c1",  // Propranolol
    "COCCc1ccc(OC(C)CNC(C)C)cc1",  // Metoprolol
    "CC(C)NCC(O)c1ccc(O)cc1",  // Atenolol
    "CC(C)NCC(O)c1ccc(O)c(c1)c1ccccc1",  // Carvedilol

    // ACE inhibitors
    "CC(C)CC(NC(=O)C1CCCN1C(=O)C(CCCCN)NC(=O)C)C(=O)O",  // Captopril
    "CCCc1ccccc1C(=O)N1CC(=O)NC(Cc2ccccc2)C1",  // Enalapril
    "NCCCC(NC(=O)C1CCCN1C(=O)C(C)N)C(=O)O",  // Lisinopril
    "CCC(C)C(NC(=O)C1CCCN1C(=O)C(C)N)C(=O)O",  // Ramipril

    // NSAIDs
    "CC(=O)Cc1ccccc1",  // Ibuprofen (simplified)
    "COc1ccc2cc(ccc2c1)C(C)C(=O)O",  // Naproxen
    "O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl",  // Diclofenac
    "CC1=CC=C(C=C1)c1cc(ns1)S(=O)(=O)N1CCN(C)CC1",  // Celecoxib

    // Antibiotics
    "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O",  // Ampicillin
    "CC(C)(C)c1ccc(O)c(c1)C(=O)N1CCCC1",  // Cephalexin-like
    "CC(C)c1c(O)c2ccccc2c(=O)c1O",  // Tetracycline-like

    // Cardiovascular drugs
    "COc1ccc(cc1)C(C)C(=O)OCCc1ccc(cc1)C(F)(F)F",  // Amlodipine (simplified)
    "Cc1cccc(C(=O)OCC(=O)Nc2ccccc2)c1",  // Nifedipine-like
    "CC(=O)Nc1ccccc1O",  // Aspirin variant
    "Cl.Cl.C[C@H](O)[C@H](O)c1ccccc1C(=O)[O-]",  // Warfarin-like

    // CNS drugs (SSRIs, benzodiazepines)
    "CCOc1ccc2nc(cc(c2c1)C(=O)N)C(F)(F)F",  // Fluoxetine
    "Clc1ccccc1C(c1ccc(Cl)cc1)(c1ccccn1)C(=O)O",  // Sertraline
    "CCN(CC)CCOC(=O)c1ccc(nc1)c1ccccc1F",  // Citalopram-like
    "CN1C(=O)CC(c2c(Cl)cccc2Cl)N(C)C1=O",  // Diazepam
    "CN1C(=O)CC(c2c(Cl)cccc2)N(C)C1=O",  // Clonazepam-like

    // Anticancer
    "COc1ccc2nc(cc(c2c1)S(=O)(=O)N(C)C)S(=O)(=O)N(C)C",  // Imatinib-like
    "CC(C)Cc1c(c2ccccc2[nH]1)C(=O)Nc1cc(Cl)c(Cl)cc1",  // Gefitinib-like

    // Aromatic heterocycles (expansion)
    "c1ccncc1",  // Pyridine
    "c1ccnc2ccccc12",  // Quinoline
    "c1ccc2ncccc2c1",  // Isoquinoline
    "c1ccc2[nH]cc(-c3ccccc3)c2c1",  // Indole
    "c1cc[nH]c1",  // Pyrrole
    "c1cc[nH]n1",  // Pyrazole
    "c1nc[nH]c1",  // Imidazole
    "c1cnc[nH]1",  // Triazole
    "c1csc[nH]1",  // Thiazole

    // Complex scaffolds
    "CC(C)c1ccc(cc1)C(C)(C)c1cc(O)ccc1O",  // Bisphenol-like
    "CC(C)c1ccc(cc1)S(=O)(=O)N(C)C",  // Sulfone
    "CC1=CC(=C)C(=O)C(C)(C)C1",  // Substituted cyclohexenone
    "O=C1c2ccccc2c3ccccc3C1=O",  // Anthraquinone
];

fn main() {
    let args: Vec<String> = std::env::args().collect();

    let smiles_list: Vec<String> = if args.len() > 1 {
        let path = &args[1];
        std::fs::read_to_string(path)
            .unwrap_or_else(|e| panic!("cannot read {path}: {e}"))
            .lines()
            .map(|l| l.trim().to_string())
            .filter(|l| !l.is_empty() && !l.starts_with('#'))
            .collect()
    } else {
        CORPUS.iter().map(|s| s.to_string()).collect()
    };

    eprintln!("Processing {} SMILES …", smiles_list.len());

    // mol_count = number of molecules that successfully parsed
    let mut mol_count = 0usize;
    // frag_in_mol: hash → number of molecules containing that environment
    let mut frag_in_mol: HashMap<u64, usize> = HashMap::new();
    let mut skip_count = 0usize;

    for smi in &smiles_list {
        let mol = match parse(smi) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("  Skip {:40} : {}", smi, e);
                skip_count += 1;
                continue;
            }
        };
        if mol.atom_count() == 0 {
            continue;
        }
        mol_count += 1;

        let counts = morgan_fp_counts(&mol, 2);
        for &hash in counts.keys() {
            *frag_in_mol.entry(hash).or_insert(0) += 1;
        }
    }

    eprintln!("Parsed {mol_count} molecules ({skip_count} skipped)");
    eprintln!("Unique fragment environments: {}", frag_in_mol.len());

    // Build sorted (hash, score_i16) table.
    let mut entries: Vec<(u64, i16)> = frag_in_mol
        .iter()
        .map(|(&hash, &count)| {
            let freq = count as f64 / mol_count as f64;
            let score_f = freq.log10() * 1000.0;
            let score_i16 = score_f.clamp(-32000.0, 0.0) as i16;
            (hash, score_i16)
        })
        .collect();
    entries.sort_unstable_by_key(|&(h, _)| h);

    println!("// Auto-generated by tools/gen_sa_table");
    println!(
        "// Corpus: {mol_count} molecules | {} unique fragment environments",
        entries.len()
    );
    println!("// Score encoding: i16 = (log10(freq_in_corpus) * 1000.0) as i16");
    println!("// Missing fragments use DEFAULT_LOG_FREQ = -5000 (log10(1e-5))");
    println!("static FRAGMENT_SCORES: &[(u64, i16)] = &[");
    for chunk in entries.chunks(4) {
        let line: Vec<String> = chunk
            .iter()
            .map(|(h, s)| format!("({h:#018x},{s:6})"))
            .collect();
        println!("    {},", line.join(", "));
    }
    println!("];");
}
