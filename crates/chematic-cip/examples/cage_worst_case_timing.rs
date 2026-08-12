//! One-off timing measurement for the Milestone 4A-2 PR: the 5 distinct
//! three-armed-cage molecules are exactly this crate's own documented
//! `needs_pass2_or_3` worst-case bucket (see `docs/rfcs/cip_accurate_rfc.md`'s
//! MANCUDE-Decision-A0 entry, `36,198` comparisons on an adamantane-cage amide).
//! Prints each molecule's wall-clock time for `assign_cip_accurate_experimental`,
//! run 20x each and reporting the minimum (least noisy) time, so a before/after
//! comparison against `main` isn't dominated by one-shot JIT/allocator noise.

use std::time::Instant;

use chematic_cip::{CipBudget, assign_cip_accurate_experimental};

const CAGE_MOLECULES: &[&str] = &[
    "CCCCCCCC/C=C/CCCCCCCC(=O)OCc1cc(=O)c(OC(=O)[C@]23C[C@H]4C[C@H](C[C@H](C4)C2)C3)co1",
    "COc1ccccc1N1CCN(CCCCNC(=O)C23C[C@H]4C[C@@H](C2)C[C@@H](C3)C4)CC1",
    "O=C(NCCCCN1CCCC(/C=C\\c2ccccc2)C1)C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3",
    "O=C(NCCCCN1CCN(c2cccc3ccccc23)CC1)C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3",
    "O=c1cc(COC2CCOCC2)occ1OC(=O)[C@]12C[C@H]3C[C@H](C[C@H](C3)C1)C2",
];

fn main() {
    let mut worst_us: u128 = 0;
    for smi in CAGE_MOLECULES {
        let mol = chematic_smiles::parse(smi).expect("valid SMILES");
        let mut best_ns = u128::MAX;
        for _ in 0..20 {
            let start = Instant::now();
            let _ = assign_cip_accurate_experimental(&mol, CipBudget::default_budget())
                .expect("assignment succeeds");
            best_ns = best_ns.min(start.elapsed().as_nanos());
        }
        let best_us = best_ns / 1000;
        worst_us = worst_us.max(best_us);
        println!("{best_us:>8} us  (best of 20)  {smi}");
    }
    println!("worst_of_cage_molecules_us={worst_us}");
}
