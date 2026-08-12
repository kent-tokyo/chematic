//! Milestone 4C-1: diagnosis (not a fix) of the 2 remaining phosphorus `skip:tied`
//! rows in the post-M4B-2 residual, distinct from Milestone 4C-0's 9 "wrong" rows.
//! Both rows are in the same molecule, `CNP1(NC)=N[P@](NC)(N2CC2)=NP(NC)(NC)=N
//! [P@@](NC)(N2CC2)=N1` (an 8-membered P4N4 cyclophosphazene), tied at atoms 6 and 19.
//!
//! Result: chematic ties consistently (both atoms, both Kekule spellings) -- confirmed
//! via `assign_cip_accurate_experimental`. Neither RDKit oracle has a stable answer for
//! this molecule either: both modern `rdCIPLabeler` and legacy `_CIPCode` flip under
//! the same Kekule respell that leaves chematic's tie unchanged (see
//! `docs/rfcs/cip_accurate_rfc.md` Milestone 4C-1 for the full writeup, including the
//! digraph trace of *why* chematic ties -- a chain-length-1 degenerate case for Rule
//! 4b's Like/Unlike operator, not a 3+-way tie or a `nearest_embedded` ambiguity).
//!
//! Usage: cargo run -p chematic-cip --release --example phosphorus_tied_diagnosis

use chematic_cip::CipBudget;
use chematic_core::{BondOrder, Molecule};

const SMILES: &str = "CNP1(NC)=N[P@](NC)(N2CC2)=NP(NC)(NC)=N[P@@](NC)(N2CC2)=N1";
const ATOMS: [u32; 2] = [6, 19];

/// Same helper as `phosphorus_kekule_diagnosis.rs`: flip every P/N ring bond
/// Single<->Double, a chemically neutral Kekule respelling of the identical molecule
/// (structural identity cross-checked externally via matching InChI, see the RFC).
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
    let mut result = mol.with_bond_order(bidx0, flip(bond0.order));
    for i in 1..n {
        let a = ring[i];
        let b = ring[(i + 1) % n];
        let (bidx, bond) = result.bond_between(a, b).expect("ring bond exists");
        result = result.with_bond_order(bidx, flip(bond.order));
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

fn main() {
    let mol = chematic_smiles::parse(SMILES).expect("valid SMILES");
    let respelled = respell_pn_ring(&mol);
    let budget = CipBudget::default_budget();

    println!("=== chematic answer, both spellings, both atoms ===");
    let mut all_tied = true;
    for (label, m) in [("orig", &mol), ("respelled", &respelled)] {
        let result = chematic_cip::assign_cip_accurate_experimental(m, budget);
        for &atom_idx in &ATOMS {
            let assigned = result.as_ref().ok().and_then(|a| {
                a.assignments
                    .iter()
                    .find(|(ai, _)| ai.0 == atom_idx)
                    .map(|(_, c)| format!("{c:?}"))
            });
            let skipped = result.as_ref().ok().and_then(|a| {
                a.skipped
                    .iter()
                    .find(|(ai, _)| ai.0 == atom_idx)
                    .map(|(_, r)| format!("{r:?}"))
            });
            if assigned.is_some() {
                all_tied = false;
            }
            println!("  {label} atom{atom_idx}: assigned={assigned:?} skipped={skipped:?}");
        }
    }
    println!("\nchematic ties on both atoms, both spellings: {all_tied}");
    println!(
        "(digraph trace -- see docs/rfcs/cip_accurate_rfc.md Milestone 4C-1 -- shows this is a\n \
         chain-length-1 degenerate case for Rule 4b's Like/Unlike operator: both tied\n \
         branches' nearest embedded stereocenter is the SAME underlying atom, reached via\n \
         two different ring paths, so there is no second chain element to compare)"
    );
}
