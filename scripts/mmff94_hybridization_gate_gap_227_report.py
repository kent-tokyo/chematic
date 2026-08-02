#!/usr/bin/env python3
"""Issue #227 Phase 0.3: joins chematic's `total_degree(atom) > 3`
hybridization-gate approximation against RDKit's real `getHybridization() !=
SP2` gate (per-atom, on the 265-molecule Wave 1 corpus), classifying every
ring C/N heavy atom into exactly one of four exclusive buckets:

  - same_decision: chematic's gate verdict matches RDKit's real verdict
  - rdkit_rejects_chematic_accepts: RDKit's real hybridization rejects the
    atom's ring (non-SP2), but chematic's total_degree<=3 approximation
    does not catch it (the approximation under-triggers)
  - rdkit_accepts_chematic_rejects: chematic's total_degree>3 approximation
    fires, but the atom is genuinely SP2 per RDKit (the approximation
    over-triggers)
  - oracle_unavailable: RDKit or chematic could not process this molecule
    (parse/Kekulize failure)

unclassified must be 0 -- every ring C/N atom present in either dump must
land in exactly one bucket above.

Run: .venv/bin/python scripts/mmff94_hybridization_gate_gap_227_report.py
"""

import json

RDKIT_PATH = "validation/results/mmff94_hybridization_gate_gap_227_rdkit.jsonl"
CHEMATIC_PATH = "validation/results/mmff94_hybridization_gate_gap_227_chematic.jsonl"


def load(path):
    rows = {}
    with open(path) as f:
        for line in f:
            row = json.loads(line)
            rows[row["name"]] = row
    return rows


def main():
    rdkit_rows = load(RDKIT_PATH)
    chematic_rows = load(CHEMATIC_PATH)

    names = sorted(set(rdkit_rows) | set(chematic_rows))

    buckets = {
        "same_decision": 0,
        "rdkit_rejects_chematic_accepts": 0,
        "rdkit_accepts_chematic_rejects": 0,
        "oracle_unavailable": 0,
    }
    total_atoms = 0
    examples = {k: [] for k in buckets}

    for name in names:
        r = rdkit_rows.get(name)
        c = chematic_rows.get(name)

        if r is None or r.get("status") != "ok" or c is None or c.get("status") != "ok":
            # Molecule-level oracle unavailability: count every ring C/N
            # atom chematic *did* manage to enumerate (if any) as
            # oracle_unavailable so the ledger stays exhaustive even when
            # one side failed outright.
            n = len(c["atoms"]) if (c and c.get("status") == "ok") else 0
            if n == 0 and r and r.get("status") == "ok":
                n = len(r["atoms"])
            buckets["oracle_unavailable"] += n
            total_atoms += n
            if n and len(examples["oracle_unavailable"]) < 5:
                examples["oracle_unavailable"].append(name)
            continue

        r_by_idx = {a["index"]: a for a in r["atoms"]}
        c_by_idx = {a["index"]: a for a in c["atoms"]}
        all_idx = set(r_by_idx) | set(c_by_idx)

        for idx in all_idx:
            total_atoms += 1
            ra = r_by_idx.get(idx)
            ca = c_by_idx.get(idx)
            if ra is None or ca is None:
                # Present in one engine's ring-atom set but not the
                # other's -- a real, distinct finding (ring-perception
                # divergence, out of Phase 0.3's scope), not silently
                # dropped: counted as oracle_unavailable since neither
                # side has a directly comparable gate verdict here.
                buckets["oracle_unavailable"] += 1
                if len(examples["oracle_unavailable"]) < 5:
                    examples["oracle_unavailable"].append(f"{name}#{idx} (ring-set mismatch)")
                continue

            rdkit_reject = ra["gate_fires_reject"]
            chematic_reject = ca["gate_fires_reject"]

            if rdkit_reject == chematic_reject:
                buckets["same_decision"] += 1
            elif rdkit_reject and not chematic_reject:
                buckets["rdkit_rejects_chematic_accepts"] += 1
                if len(examples["rdkit_rejects_chematic_accepts"]) < 8:
                    examples["rdkit_rejects_chematic_accepts"].append(
                        f"{name}#{idx} ({ca['element']}, td={ca['total_degree']}, "
                        f"rdkit_hyb={ra['hybridization']})"
                    )
            else:
                buckets["rdkit_accepts_chematic_rejects"] += 1
                if len(examples["rdkit_accepts_chematic_rejects"]) < 8:
                    examples["rdkit_accepts_chematic_rejects"].append(
                        f"{name}#{idx} ({ca['element']}, td={ca['total_degree']}, "
                        f"rdkit_hyb={ra['hybridization']})"
                    )

    unclassified = total_atoms - sum(buckets.values())

    print(json.dumps({"total_ring_c_n_atoms": total_atoms, "buckets": buckets, "unclassified": unclassified}, indent=2))
    print()
    for k, v in examples.items():
        if v:
            print(f"-- {k} examples --")
            for e in v:
                print(f"  {e}")


if __name__ == "__main__":
    main()
