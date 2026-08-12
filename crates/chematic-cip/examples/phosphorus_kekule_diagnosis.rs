//! Milestone 4C-0: mechanical diagnosis of the 9 "wrong" phosphorus rows in the
//! post-M4B-2 residual (`validation/cip_m4b2_post_port_residual.jsonl`). Diagnosis
//! only, no implementation change.
//!
//! All 9 rows are cyclophosphazene-family rings (P3N3 or P4N4, alternating P=N
//! bonds around the ring, not flagged aromatic by either chematic or RDKit). This
//! script tests Kekule-respelling invariance: for each row, it flips every ring
//! bond Single<->Double (a chemically neutral resonance-structure respelling of the
//! *same* molecule -- verified externally via matching InChI in a companion Python
//! check, see `docs/rfcs/cip_accurate_rfc.md` M4C-0 section) and re-runs
//! `assign_cip_accurate_experimental` on both forms.
//!
//! Result: chematic's answer is identical across both spellings for all 9 rows,
//! both with and without MANCUDE treatment engaged (confirming the MANCUDE gate
//! does not fire on these non-aromatic rings at all). Cross-checked against RDKit
//! separately: modern `rdCIPLabeler` FLIPS its answer under the same respelling on
//! all 9 rows, while legacy `_CIPCode` is stable and agrees with chematic on all 9.
//!
//! Usage: cargo run -p chematic-cip --release --example phosphorus_kekule_diagnosis

use chematic_cip::{AccurateCipAssignment, CipBudget};
use chematic_core::{BondOrder, CipCode, Molecule};

const ROWS: &[(&str, u32, &str)] = &[
    (
        "C1CCN(P2(N3CCCC3)=N[P@@](N3CCCC3)(N3CC3)=N[P@](N3CCCC3)(N3CC3)=N2)C1",
        11,
        "R",
    ),
    (
        "C1CCN([P@@]2(N3CC3)=NP(N3CC3)(N3CC3)=N[P@@](N3CCCCC3)(N3CC3)=N2)CC1",
        4,
        "S",
    ),
    (
        "C1CN(P2(N3CCOCC3)=N[P@](N3CCOCC3)(N3CC3)=N[P@@](N3CCOCC3)(N3CC3)=N2)CCO1",
        22,
        "S",
    ),
    (
        "C1CN([P@@]2(N3CC3)=NP(N3CC3)(N3CC3)=N[P@@](N3CCOCC3)(N3CC3)=N2)CCO1",
        16,
        "S",
    ),
    (
        "CN(C)P1(N(C)C)=NP(N(C)C)(N(C)C)=N[P@@](N(C)C)(N2CC2)=N[P@](N(C)C)(N2CC2)=N1",
        16,
        "R",
    ),
    ("CNP1(NC)=N[P@](NC)(N2CC2)=N[P@](NC)(N2CC2)=N1", 6, "S"),
    ("CN[P@@]1(N)=NP(N2CC2)(N2CC2)=N[P@](N)(NC)=N1", 2, "R"),
    (
        "Cl[P@@]1(N2CCCC2)=NP(N2CC2)(N2CC2)=N[P@@](Cl)(N2CCCC2)=N1",
        16,
        "S",
    ),
    ("N[P@]1(Cl)=NP(N2CC2)(N2CC2)=N[P@](N)(Cl)=N1", 12, "S"),
];

/// Find the P/N ring (>=3 P atoms) and flip every one of its bonds Single<->Double.
fn respell_pn_ring(mol: &Molecule) -> Molecule {
    let sssr = chematic_perception::find_sssr(mol);
    let ring = sssr
        .rings()
        .iter()
        .find(|ring| {
            ring.iter()
                .filter(|&&a| mol.atom(a).element.symbol() == "P")
                .count()
                >= 3
        })
        .expect("P/N ring present");

    let n = ring.len();
    let (bidx0, bond0) = mol
        .bond_between(ring[0], ring[1])
        .expect("ring bond exists");
    let flipped0 = flip(bond0.order);
    let mut result = mol.with_bond_order(bidx0, flipped0);
    for i in 1..n {
        let a = ring[i];
        let b = ring[(i + 1) % n];
        let (bidx, bond) = result.bond_between(a, b).expect("ring bond exists");
        let flipped = flip(bond.order);
        result = result.with_bond_order(bidx, flipped);
    }
    result
}

fn flip(order: BondOrder) -> BondOrder {
    match order {
        BondOrder::Single => BondOrder::Double,
        BondOrder::Double => BondOrder::Single,
        other => other,
    }
}

fn code_at(
    assignment: &Result<AccurateCipAssignment, chematic_cip::CipCompareError>,
    idx: u32,
) -> Option<CipCode> {
    assignment
        .as_ref()
        .ok()
        .and_then(|a| a.assignments.iter().find(|(ai, _)| ai.0 == idx))
        .map(|(_, code)| *code)
}

fn main() {
    let budget = CipBudget::default_budget();
    println!(
        "{:<4} {:<8} {:<10} {:<10} {:<10} {:<10} invariant?",
        "atom", "chematic", "orig_with", "orig_woMC", "respelled", "respelledWo"
    );
    let mut all_invariant = true;
    for &(smi, atom_idx, chematic_expected) in ROWS {
        let mol = chematic_smiles::parse(smi).expect("valid SMILES");
        let respelled = respell_pn_ring(&mol);

        let orig_with = chematic_cip::assign_cip_accurate_experimental(&mol, budget);
        let orig_without =
            chematic_cip::assign_cip_accurate_experimental_without_mancude(&mol, budget);
        let resp_with = chematic_cip::assign_cip_accurate_experimental(&respelled, budget);
        let resp_without =
            chematic_cip::assign_cip_accurate_experimental_without_mancude(&respelled, budget);

        let ow = code_at(&orig_with, atom_idx);
        let owo = code_at(&orig_without, atom_idx);
        let rw = code_at(&resp_with, atom_idx);
        let rwo = code_at(&resp_without, atom_idx);

        let invariant = ow == rw && owo == rwo && ow == owo;
        all_invariant &= invariant;

        println!(
            "{:<4} {:<8} {:<10?} {:<10?} {:<10?} {:<10?} {}",
            atom_idx,
            chematic_expected,
            ow,
            owo,
            rw,
            rwo,
            if invariant { "yes" } else { "NO" }
        );
    }
    println!(
        "\nAll 9 rows Kekule-invariant under P/N ring respelling: {}",
        all_invariant
    );
    println!(
        "(MANCUDE gate never fires on these rings: orig_with == orig_woMC in every row above)"
    );
}
