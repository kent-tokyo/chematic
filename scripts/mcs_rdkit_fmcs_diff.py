"""MCS vs. RDKit FMCS differential measurement (99-point directive Phase 2).

Pairs consecutive molecules within each established corpus (deterministic, no new
randomness), runs chematic's find_mcs and RDKit's rdFMCS.FindMCS with equivalent
default config, and compares the resulting MCS atom count. A generous per-side
timeout separates "exhaustive" results (both sides proved optimal) from "timed out"
ones; only exhaustive-vs-exhaustive mismatches are treated as informative.

Both engines' defaults are NOT directly comparable out of the box -- see the two
confounds normalized below, both config/parsing mismatches rather than algorithm
differences:

1. RDKit's own default bondCompare (CompareOrder) is lenient (single matches
   aromatic); chematic's similarly-named default (BondCompare::OrderOrAromatic) is
   actually a strict exact-match, semantically equivalent to RDKit's
   CompareOrderExact. Passing bondCompare=CompareOrderExact to RDKit here matches
   chematic's true default semantics.
2. chematic's own parser AND canonical-SMILES writer preserve literal Kekule bond
   notation verbatim (round-trip fidelity), while RDKit's parser always normalizes
   ring bonds to aromatic bond type + lowercase atoms regardless of input notation.
   Normalizing both molecules through RDKit's own Chem.MolToSmiles before parsing
   into EITHER engine removes this confound.

Known, unfixed, real divergence source (as of 2026-08-29, tracked as a Parked Item,
not fixed autonomously -- changes chematic's default MCS atom-comparator semantics):
chematic's AtomCompare::Elements requires both atomic number AND aromaticity flag to
match; RDKit's identically-named CompareElements matches on atomic number alone,
regardless of aromaticity (confirmed via live oracle, including cases where both
input molecules fully agree on aromaticity -- RDKit still omits the atom-level
aromaticity primitive entirely, encoding aromaticity only via bond-type queries).
This causes chematic's find_mcs to return None or an undersized MCS relative to
RDKit whenever two matched atoms' aromaticity disagrees. See
memory/mcs_atom_compare_elements_aromaticity_parked.md for the parked writeup.
"""
import sys
import time

import chematic
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

TIMEOUT_S = 2.0


def consecutive_pairs(smiles_list, max_pairs):
    pairs = []
    i = 0
    while i + 1 < len(smiles_list) and len(pairs) < max_pairs:
        pairs.append((smiles_list[i], smiles_list[i + 1]))
        i += 2
    return pairs


def run_pair(smi1, smi2):
    rd1 = Chem.MolFromSmiles(smi1)
    rd2 = Chem.MolFromSmiles(smi2)
    if rd1 is None or rd2 is None:
        return None
    norm1 = Chem.MolToSmiles(rd1)
    norm2 = Chem.MolToSmiles(rd2)
    rd1 = Chem.MolFromSmiles(norm1)
    rd2 = Chem.MolFromSmiles(norm2)
    try:
        cm1 = chematic.from_smiles(norm1)
        cm2 = chematic.from_smiles(norm2)
    except Exception:
        return None

    rd_result = rdFMCS.FindMCS(
        [rd1, rd2], timeout=int(TIMEOUT_S), bondCompare=rdFMCS.BondCompare.CompareOrderExact
    )
    rd_atoms = rd_result.numAtoms
    rd_timed_out = rd_result.canceled

    mcs, chem_timed_out = chematic.find_mcs_checked(
        [cm1, cm2], timeout_ms=int(TIMEOUT_S * 1000)
    )
    chem_atoms = mcs.heavy_atoms if mcs is not None else 0

    return {
        "smi1": smi1,
        "smi2": smi2,
        "rd_atoms": rd_atoms,
        "chem_atoms": chem_atoms,
        "rd_timed_out": rd_timed_out,
        "chem_timed_out": chem_timed_out,
    }


def main():
    files = [
        "scripts/descriptor_census_corpus.smi",
        "scripts/chembl_accuracy_corpus_4999.smi",
        "scripts/nci_first_5k_smiles_only.smi",
    ]
    max_pairs_per_corpus = int(sys.argv[1]) if len(sys.argv) > 1 else 200

    for f in files:
        with open(f) as fh:
            lines = [l.split()[0] for l in fh if l.strip()]
        pairs = consecutive_pairs(lines, max_pairs_per_corpus)

        exhaustive_match = 0
        exhaustive_total = 0
        either_timed_out = 0
        mismatches = []
        errors = 0

        t0 = time.time()
        for smi1, smi2 in pairs:
            try:
                r = run_pair(smi1, smi2)
            except Exception:
                errors += 1
                continue
            if r is None:
                errors += 1
                continue
            if r["rd_timed_out"] or r["chem_timed_out"]:
                either_timed_out += 1
                continue
            exhaustive_total += 1
            if r["rd_atoms"] == r["chem_atoms"]:
                exhaustive_match += 1
            else:
                mismatches.append(r)
        elapsed = time.time() - t0

        print(f"\n=== {f} ===")
        print(f"pairs: {len(pairs)}  errors: {errors}  either_timed_out: {either_timed_out}")
        print(f"exhaustive: {exhaustive_match}/{exhaustive_total} match  ({elapsed:.1f}s)")
        for m in mismatches[:10]:
            direction = "chematic>RDKit" if m["chem_atoms"] > m["rd_atoms"] else "chematic<RDKit"
            print(f"  MISMATCH ({direction}): rd={m['rd_atoms']} chem={m['chem_atoms']}  {m['smi1']!r} | {m['smi2']!r}")
        if len(mismatches) > 10:
            print(f"  ... and {len(mismatches) - 10} more mismatches")


if __name__ == "__main__":
    main()
