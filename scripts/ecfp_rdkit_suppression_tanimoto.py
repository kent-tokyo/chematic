#!/usr/bin/env python3
"""Phase B Tanimoto-correlation-vs-RDKit before/after (non-gating reference
measurement). Same sampling frame as PR #120's
`ecfp_rdkit_environment_parity.py::tanimoto_correlation`: deterministic
sample of up to 300 molecules (seed=42) from the population, all pairwise
(i<j) Tanimoto similarities within the sample, Pearson correlation between
chematic's pairwise-similarity series and RDKit's own.

Computed twice on the identical sample/seed:
  - "before": chematic's existing `ecfp4_rdkit_invariants()` path (no
    suppression) vs RDKit's real `default.folded_on_bits`.
  - "after": this PR's new `ecfp4_rdkit_environment_experimental()`
    (suppression) path vs the same RDKit reference.

Reports both correlations and whether suppression improved, worsened, or
left the correlation with RDKit's real fingerprint unchanged. This is a
reference measurement, not a merge gate.

Usage:
    python scripts/ecfp_rdkit_suppression_tanimoto.py \
        --chematic <morgan_suppression_tanimoto_dump.jsonl> \
        --rdkit-oracle <gen_ecfp_rdkit_environment_oracle.py --rows-out output> \
        --summary-out <out.json>
"""

from __future__ import annotations

import argparse
import json
import random
import statistics
import sys


def load_jsonl(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def _tanimoto(a: frozenset, b: frozenset) -> float:
    if not a and not b:
        return 1.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


def tanimoto_correlation(chem_bit_field, chem_rows, rd_rows, sample_size=300, seed=42):
    """`chem_bit_field` selects which chematic series to correlate:
    "baseline_folded_bits" or "suppression_folded_bits"."""
    pairs = []
    for chem, rd in zip(chem_rows, rd_rows):
        if not chem.get("parse_ok") or not rd.get("parse_ok"):
            continue
        chem_bits = frozenset(chem[chem_bit_field])
        rd_bits = frozenset(rd["default"]["folded_on_bits"])
        pairs.append((chem_bits, rd_bits))

    rng = random.Random(seed)
    sample = pairs if len(pairs) <= sample_size else rng.sample(pairs, sample_size)

    chem_sims, rd_sims = [], []
    n = len(sample)
    for i in range(n):
        for j in range(i + 1, n):
            chem_sims.append(_tanimoto(sample[i][0], sample[j][0]))
            rd_sims.append(_tanimoto(sample[i][1], sample[j][1]))

    pearson_r = round(statistics.correlation(chem_sims, rd_sims), 4) if len(chem_sims) >= 2 else None
    return {
        "pearson_r": pearson_r,
        "population_molecules": len(pairs),
        "sample_molecules": n,
        "pair_count": len(chem_sims),
        "seed": seed,
        "gating": False,
    }


def run(chematic_rows, rdkit_rows):
    if len(chematic_rows) != len(rdkit_rows):
        print(
            f"PIPELINE ERROR: row count mismatch chematic={len(chematic_rows)} "
            f"rdkit={len(rdkit_rows)}",
            file=sys.stderr,
        )
        sys.exit(1)
    for idx, (chem, rd) in enumerate(zip(chematic_rows, rdkit_rows)):
        if chem.get("row_id") != idx or rd.get("row_id") != idx:
            print(f"PIPELINE ERROR at position {idx}: row_id out of sync", file=sys.stderr)
            sys.exit(1)
        if chem.get("smiles") != rd.get("smiles"):
            print(
                f"PIPELINE ERROR at row {idx}: chematic smiles={chem.get('smiles')!r} "
                f"!= rdkit smiles={rd.get('smiles')!r}",
                file=sys.stderr,
            )
            sys.exit(1)

    before = tanimoto_correlation("baseline_folded_bits", chematic_rows, rdkit_rows)
    after = tanimoto_correlation("suppression_folded_bits", chematic_rows, rdkit_rows)

    delta = None
    verdict = None
    if before["pearson_r"] is not None and after["pearson_r"] is not None:
        delta = round(after["pearson_r"] - before["pearson_r"], 4)
        if delta > 0.0001:
            verdict = "improved"
        elif delta < -0.0001:
            verdict = "worsened"
        else:
            verdict = "unchanged"

    return {
        "before_baseline_ecfp4_rdkit_invariants": before,
        "after_suppression_ecfp4_rdkit_environment_experimental": after,
        "delta": delta,
        "verdict": verdict,
        "reference": "RDKit default.folded_on_bits (rdFingerprintGenerator.GetMorganGenerator, includeRedundantEnvironments=False)",
        "gating": False,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chematic", required=True)
    p.add_argument("--rdkit-oracle", required=True)
    p.add_argument("--summary-out", default=None)
    args = p.parse_args()

    chematic_rows = load_jsonl(args.chematic)
    rdkit_rows = load_jsonl(args.rdkit_oracle)
    summary = run(chematic_rows, rdkit_rows)

    print(json.dumps(summary, indent=2))
    if args.summary_out:
        with open(args.summary_out, "w") as f:
            json.dump(summary, f, indent=2, sort_keys=True)
        print(f"summary written to {args.summary_out}")


if __name__ == "__main__":
    main()
