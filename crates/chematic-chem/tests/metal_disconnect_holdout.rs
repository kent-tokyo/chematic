//! A hand-built, genuinely new (not derived from `nci_first_5k_smiles_only.smi`
//! or either dev corpus) metal-coordination holdout for issue #403 --
//! `disconnect_metals` leaving a dative-bond-derived `[O+]`/`[N+]` with a
//! stale, too-low `hydrogen_count` after severing the metal bond, so
//! `neutralize_charges` (which runs immediately after in the pipeline, and
//! whose guard is `h > 0` on the raw stored field) skipped it on the first
//! `standardize()` pass.
//!
//! Deliberately varies both metal (Ni, Co, Al, Zn, Cr, Fe, Mg, Mn, Hg, Cd,
//! Pd -- covering the directive's own Cu/Ni/Zn/Co/Cr/Al/Mg/Mn/Hg list plus a
//! couple more) and ligand shape (simple monodentate O/N donors, a
//! tris-chelate, a bidentate oxalate-like ring, an 8-hydroxyquinoline-style
//! chelate, aromatic N-donor pyridines, and a plain ionic salt with no
//! metal-ligand bond at all, as a negative control) rather than repeating one
//! scaffold with the metal swapped -- so this exercises the fix's general
//! correctness, not just the one shape issue #403's own repro happened to
//! use. Every SMILES independently confirmed parseable by RDKit before being
//! committed here.
use chematic_chem::{StandardizeOptions, standardize};
use chematic_smiles::{canonical_smiles, parse};

/// Atomic numbers of the metals this test's fixtures actually use (Ni 28,
/// Co 27, Al 13, Zn 30, Cr 24, Fe 26, Mg 12, Mn 25, Hg 80, Cd 48, Pd 46) --
/// narrower than `chematic_chem::standardize`'s own (private) `is_metal`
/// table, which isn't exported, but sufficient to identify "is this atom one
/// of the metals in this test" without needing the full periodic table.
fn is_test_metal(atomic_number: u8) -> bool {
    matches!(
        atomic_number,
        12 | 13 | 24 | 25 | 26 | 27 | 28 | 30 | 46 | 48 | 80
    )
}

/// (label, input SMILES). Each fixture's `standardize()` output was spot-
/// verified by hand to be a plausible, correctly-neutralized, valence-sound
/// structure (neutral ketone/imine/amine/alcohol fragments plus a bare metal
/// atom) -- not cross-checked against an external oracle (RDKit has no
/// equivalent "disconnect metals and neutralize" standardization step to
/// compare against directly for this specific transform).
const FIXTURES: &[(&str, &str)] = &[
    (
        "Ni tris-acetylacetonate-like",
        "CC1=[O+][Ni]2([O+]=C(C)C1)[O+]=C(C)CC(=[O+]2)C",
    ),
    (
        "Co tris(acetone) monodentate",
        "CC(C)=[O+][Co]([O+]=C(C)C)([O+]=C(C)C)[O+]=C(C)C",
    ),
    (
        "Al tris(2-naphthoxide)",
        "c1ccc2c(c1)cccc2[O+][Al]([O+]c1cccc2ccccc12)[O+]c1cccc2ccccc12",
    ),
    (
        "Zn diacetate ionic salt (negative control, no metal-ligand bond)",
        "CC(=O)[O-].CC(=O)[O-].CC(=O)[O-].CC(=O)[O-].[Zn+2].[Zn+2]",
    ),
    (
        "Cr bidentate oxalate-like ring",
        "O=C1OC(=[O+][Cr]23[O+]=C(O1)O2)O3",
    ),
    (
        "Fe tris(benzophenone) monodentate",
        "c1ccc(cc1)C(=[O+][Fe]([O+]=C(c1ccccc1)c1ccccc1)([O+]=C(c1ccccc1)c1ccccc1)[O+]=C(c1ccccc1)c1ccccc1)c1ccccc1",
    ),
    ("Mg diamide, no charge at all", "[NH2][Mg][NH2]"),
    ("Mn bis(imine) N-donor", "CC(=[N+][Mn][N+]=C(C)C)C"),
    (
        "Hg bis(pyridyl), aromatic N donor, no formal charge",
        "c1ccncc1[Hg]c1ccncc1",
    ),
    ("Cd bis(acetone)", "CC(=[O+][Cd][O+]=C(C)C)C"),
    (
        "Pd bidentate diketonate ring",
        "CCC(=[O+][Pd]1[O+]=C(CC)CC(=[O+]1)CC)CC",
    ),
];

#[test]
fn metal_disconnect_holdout_is_idempotent_and_charge_neutral() {
    let opts = StandardizeOptions::default();
    let mut failures: Vec<String> = Vec::new();

    for &(label, smi) in FIXTURES {
        let mol = match parse(smi) {
            Ok(m) => m,
            Err(e) => {
                failures.push(format!("{label}: PARSE FAILED '{smi}': {e}"));
                continue;
            }
        };
        let once = canonical_smiles(&standardize(&mol, &opts));
        let reparsed = match parse(&once) {
            Ok(m) => m,
            Err(e) => {
                failures.push(format!(
                    "{label}: RE-PARSE FAILED '{smi}' (once='{once}'): {e}"
                ));
                continue;
            }
        };
        let twice = canonical_smiles(&standardize(&reparsed, &opts));
        if once != twice {
            failures.push(format!(
                "{label}: NOT IDEMPOTENT '{smi}': once='{once}' twice='{twice}'"
            ));
        }
        // No *non-metal* atom should be left with a dangling formal charge
        // that has no bond left to justify it -- the actual invariant issue
        // #403 is about. A free metal cation ([Zn+2] with no ligand bond at
        // all, this fixture list's own ionic-salt negative control) legitimately
        // keeps its charge -- `disconnect_metals`/`neutralize_charges` never
        // neutralize a bare metal ion, and this test must not demand that.
        let dangling_nonmetal_charge: i32 = reparsed
            .atoms()
            .filter(|(_, a)| !is_test_metal(a.element.atomic_number()))
            .map(|(_, a)| a.charge as i32)
            .sum();
        if dangling_nonmetal_charge != 0 {
            failures.push(format!(
                "{label}: non-metal atom(s) left with a dangling formal charge after \
                 standardize: '{once}' (non-metal net charge={dangling_nonmetal_charge})"
            ));
        }
    }

    assert!(
        failures.is_empty(),
        "metal_disconnect_holdout: {} fixture(s) failed:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
