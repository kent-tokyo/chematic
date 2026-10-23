//! Isotope-separated total composition (Open Babel #1464).
use chematic_core::{Molecule, implicit_hcount};
use std::collections::BTreeMap;
use std::fmt::Write;

/// Total molecular formula with explicit isotope labels and net charge.
///
/// Uses the same hydrogen/model convention as [`crate::calc_mol_formula`],
/// adding bracket/implicit Hs and counting explicit hydrogen atoms once.
/// Within each element, unlabelled atoms precede mass-number-ordered labels:
/// `C[13C]H3O2-`, `[13C]H5[15N]`, `H[2H]O`. C and H families come first,
/// then other symbols alphabetically. D/T are written `[2H]`/`[3H]`.
/// Isotope zero means unspecified. Dummy atoms remain `*` (or `[13*]`).
/// Net charge follows existing formula syntax: `+`, `-`, `+2`, `-3`.
/// No isotope-abundance, chemical-validity, radical or structural identity is
/// implied. Atom maps, stereo and isotope locations are not encoded. The
/// existing isotope-blind formula and core explicit-atom formula are unchanged.
pub fn calc_mol_formula_with_isotopes(mol: &Molecule) -> String {
    chematic_smiles::with_rdkit_model_molecule(mol, formula)
}

fn formula(mol: &Molecule) -> String {
    let mut counts = BTreeMap::<(&'static str, Option<u16>), u64>::new();
    let mut charge = 0i64;
    for (idx, atom) in mol.atoms() {
        let symbol = if atom.wildcard {
            "*"
        } else {
            atom.element.symbol()
        };
        let isotope = atom.isotope.filter(|&mass| mass > 0);
        *counts.entry((symbol, isotope)).or_default() += 1;
        charge += i64::from(atom.charge);
        let hydrogens = if atom.wildcard {
            atom.hydrogen_count.unwrap_or(0)
        } else {
            implicit_hcount(mol, idx)
        };
        if hydrogens > 0 {
            *counts.entry(("H", None)).or_default() += u64::from(hydrogens);
        }
    }
    let mut text = String::new();
    for family in ["C", "H"] {
        for (&(symbol, isotope), &count) in &counts {
            if symbol == family {
                term(&mut text, symbol, isotope, count);
            }
        }
    }
    for (&(symbol, isotope), &count) in &counts {
        if symbol != "C" && symbol != "H" {
            term(&mut text, symbol, isotope, count);
        }
    }
    if text.is_empty() {
        text.push_str("H0");
    }
    match charge {
        0 => {}
        1 => text.push('+'),
        -1 => text.push('-'),
        n if n > 0 => write!(text, "+{n}").expect("String formatting cannot fail"),
        n => write!(text, "-{}", -n).expect("String formatting cannot fail"),
    }
    text
}

fn term(text: &mut String, symbol: &str, isotope: Option<u16>, count: u64) {
    if let Some(mass) = isotope {
        write!(text, "[{mass}{symbol}]").expect("String formatting cannot fail");
    } else {
        text.push_str(symbol);
    }
    if count > 1 {
        write!(text, "{count}").expect("String formatting cannot fail");
    }
}
