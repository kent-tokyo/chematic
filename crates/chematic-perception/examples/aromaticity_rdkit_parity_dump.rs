//! Diag/aromaticity-rdkit-parity: dump chematic's aromaticity/kekulization
//! output for a frozen, deliberately-constructed SMILES fixture corpus, one
//! JSON row per fixture, for cross-checking against RDKit in
//! `scripts/aromaticity_rdkit_parity_diagnosis.py`.
//!
//! Diagnostic only. Calls only existing public APIs
//! (`chematic_smiles::parse`, `chematic_core::kekulize`/`apply_kekule`,
//! `chematic_perception::{apply_aromaticity, apply_aromaticity_ex,
//! apply_aromaticity_rdkit_parity_experimental, find_sssr,
//! augmented_ring_set, count_aromatic_rings}`) exactly as an external caller
//! would. Does not modify any production module.
//!
//! Each fixture is fed to chematic in ONE of two notations, chosen per
//! molecule to test the more interesting direction:
//! - Aromatic (lowercase) notation, when a valid one exists: exercises BOTH
//!   `chematic_core::kekulize` (aromatic -> Kekulé bond orders, "kekulization"
//!   mechanism) AND `apply_aromaticity` re-perception from that Kekulé form
//!   ("atom/bond aromatic flags" mechanism), starting from the *same*
//!   in-memory molecule -- no second SMILES string needed.
//! - Kekulé (uppercase) notation, for molecules with no valid aromatic
//!   notation at all (most exocyclic-multiple-bond and boron cases): only
//!   the re-perception direction applies.
//!
//! Run:
//! ```text
//! cargo run -p chematic-perception --example aromaticity_rdkit_parity_dump \
//!     > validation/results/aromaticity_rdkit_parity_fixture_dump.jsonl
//! .venv/bin/python scripts/aromaticity_rdkit_parity_diagnosis.py
//! ```

use chematic_perception::{
    AromaticityAlgorithm, apply_aromaticity, apply_aromaticity_ex,
    apply_aromaticity_rdkit_parity_experimental, assign_aromaticity_ex, augmented_ring_set,
    count_aromatic_rings, find_sssr,
};
use serde_json::json;

/// (id, category, description, smiles). Frozen -- see
/// `docs/rfcs/aromaticity_rdkit_parity_rfc.md` §4 for what each category is
/// deliberately targeting and why the specific molecule was chosen.
const FIXTURES: &[(&str, &str, &str, &str)] = &[
    // --- baseline sanity ---
    (
        "benzene",
        "baseline",
        "simplest aromatic positive control",
        "c1ccccc1",
    ),
    (
        "cyclohexane",
        "negative_control",
        "fully saturated, trivial non-aromatic",
        "C1CCCCC1",
    ),
    (
        "cyclohexadiene_1_3",
        "negative_control",
        "conjugated diene but ring not fully conjugated (2 sp3 atoms)",
        "C1=CCCC=C1",
    ),
    (
        "cyclooctatetraene",
        "negative_control_antiaromatic_electron_count",
        "8-membered fully conjugated, 4n (not 4n+2) pi count",
        "C1=CC=CC=CC=C1",
    ),
    // --- heteroaromatics: classic Huckel edge cases, mixed pyrrole-/pyridine-type N ---
    (
        "pyridine",
        "heteroaromatic",
        "single pyridine-type N (1pi, no H)",
        "c1ccncc1",
    ),
    (
        "pyrimidine",
        "heteroaromatic",
        "two pyridine-type N",
        "c1cncnc1",
    ),
    (
        "furan",
        "heteroaromatic",
        "O lone-pair donor (2pi)",
        "c1ccoc1",
    ),
    (
        "thiophene",
        "heteroaromatic",
        "S lone-pair donor (2pi)",
        "c1ccsc1",
    ),
    (
        "pyrrole",
        "heteroaromatic",
        "pyrrole-type N, needs 1 implicit H",
        "c1cc[nH]c1",
    ),
    (
        "n_methylpyrrole",
        "heteroaromatic",
        "pyrrole-type N, substituted (no H) -- regression guard for the\
         pyrrole/pyridine Kekule-erasure bug fixed in aromaticity.rs",
        "Cn1cccc1",
    ),
    (
        "imidazole",
        "heteroaromatic",
        "one pyrrole-type N ([nH]) + one pyridine-type N in the SAME ring",
        "c1c[nH]cn1",
    ),
    (
        "pyrazole",
        "heteroaromatic",
        "adjacent pyrrole-type N + pyridine-type N",
        "c1cc[nH]n1",
    ),
    (
        "oxazole",
        "heteroaromatic",
        "pyridine-type N + O lone-pair donor",
        "c1ocnc1",
    ),
    (
        "thiazole",
        "heteroaromatic",
        "pyridine-type N + S lone-pair donor",
        "c1cscn1",
    ),
    (
        "isoxazole",
        "heteroaromatic",
        "pyridine-type N + O, N-O adjacency",
        "c1cnoc1",
    ),
    (
        "triazole_1_2_3",
        "heteroaromatic",
        "2 pyridine-type N + 1 pyrrole-type N",
        "c1cn[nH]n1",
    ),
    (
        "tetrazole",
        "heteroaromatic",
        "3 pyridine-type N + 1 pyrrole-type N",
        "c1nnn[nH]1",
    ),
    // --- fused / polycyclic aromatics ---
    (
        "naphthalene",
        "fused_polycyclic",
        "simplest fused bicyclic, 2 valid Kekule structures",
        "c1ccc2ccccc2c1",
    ),
    (
        "anthracene",
        "fused_polycyclic",
        "linear tricyclic",
        "c1ccc2cc3ccccc3cc2c1",
    ),
    (
        "quinoline",
        "fused_polycyclic",
        "benzo-fused pyridine",
        "c1ccc2ncccc2c1",
    ),
    (
        "isoquinoline",
        "fused_polycyclic",
        "benzo-fused pyridine, N position 2",
        "c1ccc2cnccc2c1",
    ),
    (
        "indole",
        "fused_polycyclic",
        "benzo-fused pyrrole",
        "c1ccc2[nH]ccc2c1",
    ),
    (
        "indolizine",
        "fused_polycyclic_bridged_sssr_artifact",
        "bridgehead-N bicyclic (5+6 fused sharing N); CLAUDE.md flags this shape \
         as a known SSSR fundamental-cycle artifact (a 9-ring instead of 5+6)",
        "c1ccn2ccccc12",
    ),
    (
        "purine",
        "fused_polycyclic",
        "bicyclic, 4 N total, mixed pyrrole/pyridine-type",
        "c1ncc2[nH]cnc2n1",
    ),
    (
        "azulene",
        "fused_polycyclic",
        "non-alternant 5+7 fused PAH; canonical form has an explicit \
         non-aromatic fusion bond between two aromatic rings",
        "c1ccc2cccc-2cc1",
    ),
    // --- charged aromatics ---
    (
        "tropylium_cation",
        "charged_aromatic",
        "7-ring carbocation, 6pi/7 centers",
        "c1ccc[cH+]cc1",
    ),
    (
        "cyclopentadienyl_anion",
        "charged_aromatic",
        "5-ring carbanion, 6pi/5 centers",
        "c1cc[cH-]c1",
    ),
    (
        "imidazolium",
        "charged_aromatic",
        "protonated imidazole, both N pyrrole-type",
        "c1c[nH+]c[nH]1",
    ),
    (
        "pyridinium",
        "charged_aromatic",
        "protonated pyridine",
        "c1cc[nH+]cc1",
    ),
    (
        "pyrylium",
        "charged_aromatic",
        "O-centered aromatic cation",
        "c1cc[o+]cc1",
    ),
    // --- Se / Te / P / B ---
    (
        "selenophene",
        "chalcogen_se_te",
        "Se lone-pair donor; supported only under AromaticityAlgorithm::RdkitLike, not default Huckel",
        "c1cc[se]c1",
    ),
    (
        "tellurophene",
        "chalcogen_se_te",
        "Te lone-pair donor; same RdkitLike-only gate as selenophene",
        "c1cc[te]c1",
    ),
    (
        "phosphole",
        "phosphorus_boron",
        "P lone-pair donor, RDKit's default model treats it as aromatic; \
         chematic's own doc comment says P-containing rings are NOT supported",
        "c1cc[pH]c1",
    ),
    (
        "borole",
        "phosphorus_boron",
        "B empty p-orbital, 4pi antiaromatic-count ring; no valid aromatic \
         notation exists (RDKit itself keeps this Kekule) -- fed directly as Kekule",
        "C1=CC=CB1",
    ),
    (
        "borazine",
        "phosphorus_boron",
        "inorganic benzene analog (alternating B/N, no ring double bonds at \
         all); RDKit's default model does NOT mark it aromatic",
        "B1NBNBN1",
    ),
    // --- exocyclic multiple bonds on ring atoms ---
    (
        "pyridone_2",
        "exocyclic_multiple_bond",
        "aromatic per RDKit despite exocyclic ring C=O -- carbonyl C \
         contributes 0pi, NH lone pair supplies 2pi",
        "O=c1cccc[nH]1",
    ),
    (
        "tropone",
        "exocyclic_multiple_bond",
        "7-ring, exocyclic C=O, NO heteroatom lone-pair donor -- still 6pi \
         over the other 6 ring carbons; fed as Kekule (unambiguous structure)",
        "O=C1C=CC=CC=C1",
    ),
    (
        "benzoquinone_1_4",
        "exocyclic_multiple_bond",
        "textbook non-aromatic: two exocyclic C=O carbons each contribute 0pi, \
         only 4pi remain from the two ring C=C bonds",
        "O=C1C=CC(=O)C=C1",
    ),
    (
        "cyclopentadienone",
        "exocyclic_multiple_bond",
        "5-ring, exocyclic C=O, remaining 4 ring pi electrons -> antiaromatic count",
        "O=C1C=CC=C1",
    ),
    (
        "thiophene_1_oxide",
        "exocyclic_multiple_bond",
        "exocyclic S=O ties up the ring S's lone pair -- aromaticity.rs has an \
         explicit comment about this sulfoxide/sulfone rule",
        "O=S1C=CC=C1",
    ),
];

fn bond_order_str(o: chematic_core::BondOrder) -> String {
    format!("{o:?}")
}

fn main() {
    let mut n_ok = 0usize;
    let mut n_fail = 0usize;

    for &(id, category, description, smiles) in FIXTURES {
        let raw = match chematic_smiles::parse(smiles) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("FATAL: fixture {id:?} failed to parse {smiles:?}: {e}");
                n_fail += 1;
                continue;
            }
        };

        // --- kekulization direction: aromatic bonds (if any) -> Single/Double ---
        let kekule_result = chematic_core::kekulize(&raw);
        let (kekulize_ok, kekulize_error, kek_mol) = match kekule_result {
            Ok(k) => (true, None, chematic_core::apply_kekule(&raw, &k)),
            Err(e) => (false, Some(e.to_string()), raw.clone()),
        };

        // --- re-perception direction: Kekule bond orders -> aromatic flags ---
        // `*_model_aromatic_atom_count` records the RAW Huckel-model verdict
        // (assign_aromaticity_ex, pre-`build_molecule_from_model`) directly,
        // independent of whatever the final rebuilt Molecule's atom.aromatic
        // flags end up saying -- load-bearing for detecting cases where the
        // final flags merely preserve an already-true input flag rather than
        // being confirmed by this model (see RFC finding on
        // `build_molecule_from_model`'s promote-only atom/bond loops).
        let huckel_model_aromatic_atom_count =
            assign_aromaticity_ex(&kek_mol, AromaticityAlgorithm::Huckel).aromatic_atom_count();
        let rdkitlike_model_aromatic_atom_count =
            assign_aromaticity_ex(&kek_mol, AromaticityAlgorithm::RdkitLike).aromatic_atom_count();
        let default_mol = apply_aromaticity(&kek_mol);
        let rdkitlike_mol = apply_aromaticity_ex(&kek_mol, AromaticityAlgorithm::RdkitLike);
        let experimental = apply_aromaticity_rdkit_parity_experimental(&raw);
        let (experimental_ok, experimental_error, experimental_mol) = match experimental {
            Ok(m) => (true, None, Some(m)),
            Err(e) => (false, Some(e.to_string()), None),
        };

        // --- ring perception ---
        let ring_set = find_sssr(&raw);
        let raw_sssr_rings: Vec<Vec<u32>> = ring_set
            .rings()
            .iter()
            .map(|ring| ring.iter().map(|a| a.0).collect())
            .collect();
        let augmented = augmented_ring_set(&raw, ring_set.rings());
        let augmented_ring_count = augmented.len();
        let aromatic_ring_count = count_aromatic_rings(&raw);

        let atoms: Vec<_> = raw
            .atoms()
            .map(|(idx, atom)| {
                json!({
                    "idx": idx.0,
                    "element": atom.element.symbol(),
                    "charge": atom.charge,
                    "default_aromatic": default_mol.atom(idx).aromatic,
                    "rdkitlike_aromatic": rdkitlike_mol.atom(idx).aromatic,
                    "experimental_aromatic": experimental_mol.as_ref().map(|m| m.atom(idx).aromatic),
                })
            })
            .collect();

        let bonds: Vec<_> = raw
            .bonds()
            .map(|(idx, bond)| {
                json!({
                    "idx": idx.0,
                    "a1": bond.atom1.0,
                    "a2": bond.atom2.0,
                    "default_aromatic": default_mol.bond(idx).order == chematic_core::BondOrder::Aromatic,
                    "default_order": bond_order_str(default_mol.bond(idx).order),
                    "rdkitlike_aromatic": rdkitlike_mol.bond(idx).order == chematic_core::BondOrder::Aromatic,
                    "experimental_aromatic": experimental_mol.as_ref()
                        .map(|m| m.bond(idx).order == chematic_core::BondOrder::Aromatic),
                    "kekulized_order": bond_order_str(kek_mol.bond(idx).order),
                })
            })
            .collect();

        let row = json!({
            "id": id,
            "category": category,
            "description": description,
            "smiles": smiles,
            "kekulize_ok": kekulize_ok,
            "kekulize_error": kekulize_error,
            "experimental_ok": experimental_ok,
            "experimental_error": experimental_error,
            "raw_sssr_ring_count": raw_sssr_rings.len(),
            "raw_sssr_rings": raw_sssr_rings,
            "augmented_ring_count": augmented_ring_count,
            "count_aromatic_rings": aromatic_ring_count,
            "huckel_model_aromatic_atom_count": huckel_model_aromatic_atom_count,
            "rdkitlike_model_aromatic_atom_count": rdkitlike_model_aromatic_atom_count,
            "atoms": atoms,
            "bonds": bonds,
        });
        println!("{row}");
        n_ok += 1;
    }

    eprintln!("dumped {n_ok} fixtures, {n_fail} parse failures");
    if n_fail > 0 {
        std::process::exit(1);
    }
}
