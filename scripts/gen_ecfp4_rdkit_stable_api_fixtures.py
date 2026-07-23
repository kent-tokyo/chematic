#!/usr/bin/env python3
"""Generates `validation/ecfp4_rdkit_stable_api_fixtures.json`: the single shared
fixture+expectation corpus that chematic's Rust, Python, and WASM surfaces for the
promoted `rdkit_morgan_ecfp4_experimental`/`rdkit_morgan_fingerprint` API all test
against, byte-for-byte identically. See `docs/ecfp4_bitexact_api_rfc.md` for the
diagnosis this promotes and `crates/chematic-fp/src/rdkit_morgan_ecfp4.rs` /
`crates/chematic-fp/src/rdkit_morgan_config.rs` for the production code.

All success-case expected values come directly from a live RDKit oracle -- never from
chematic's own output -- so cross-language parity tests that compare against this file
are simultaneously checking "does chematic match RDKit" (not just "do the three
bindings agree with each other, possibly on a shared wrong answer").

Fixture SMILES are drawn from the diagnosis workstream's
`scripts/ecfp4_bitexact_matrix_fixtures.csv` (ids 1-33; same molecules, not the same
file -- that CSV's tags reflect its own diagnosis-time classification, which changed
after `fix/kekulize-charge-aware-k1` fixed 6 of them), plus one fixture (34) that still
genuinely fails RDKit-parity aromaticity preprocessing on this branch -- required so
the shared corpus has real error-path coverage, not just success-path coverage.

Usage:
    .venv/bin/python scripts/gen_ecfp4_rdkit_stable_api_fixtures.py \\
        --out validation/ecfp4_rdkit_stable_api_fixtures.json
"""

from __future__ import annotations

import argparse
import json
import subprocess

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")

RDKIT_PINNED_COMMIT = "8afba32ec539dcb2369bc84549d802aca3f7eb39"  # Release_2026_03_4

PRODUCTION_RADIUS = 2
PRODUCTION_FP_SIZE = 2048
RADIUS_SWEEP = [0, 1, 2, 3]
FP_SIZE_SWEEP = [128, 256, 512, 1024, 2048]

_PARSE_PARAMS = Chem.SmilesParserParams()
_PARSE_PARAMS.removeHs = False  # keep isotope-labeled/explicit atoms as real graph atoms

# id, smiles -- same 33 molecules as scripts/ecfp4_bitexact_matrix_fixtures.csv (ids
# 1-33), all of which succeed on this branch (verified via
# `cargo run -p chematic-fp --release --features diagnostics --example
# rdkit_ecfp4_bitexact_matrix_dump`, 33/33 "success" post-K1).
SUCCESS_FIXTURES = [
    ("1", "c1ccccc1"),
    ("2", "CC"),
    ("3", "CCC"),
    ("4", "c1ccncc1"),
    ("5", "C1=CC=NC=C1"),
    ("6", "c1ccc2ccccc2c1"),
    ("7", "C1=CC2=CC=CC=C2C=C1"),
    ("8", "c1ccoc1"),
    ("9", "C1=CC=CO1"),
    ("10", "c1ccsc1"),
    ("11", "C1=CC=CS1"),
    ("12", "CC.c1ccccc1"),
    ("13", "[Na+].[Cl-].c1ccccc1"),
    ("14", "[NH4+].[Cl-]"),
    ("15", "[13CH4]"),
    ("16", "[13CH3]C(=O)O"),
    ("17", "[2H]OC"),
    ("18", "c1cc([2H])ccc1"),
    ("19", "CC(=O)[O-]"),
    ("20", "C[NH3+]"),
    ("21", "[NH4+]"),
    ("22", "c1ccc(cc1)[N+](=O)[O-]"),
    ("23", "[O-]S(=O)(=O)[O-]"),
    ("24", "c1ccc[cH+]cc1"),
    ("25", "c1c[nH+]c[nH]1"),
    ("26", "c1cc[nH+]cc1"),
    ("27", "c1cc[o+]cc1"),
    ("28", "c1cc[te]c1"),
    ("29", "c1cc[pH]c1"),
    ("30", "C[C@H](N)C(=O)O"),
    ("31", "C[C@@H](N)C(=O)O"),
    ("32", "C/C=C/C"),
    ("33", "C/C=C\\C"),
]

# Bridgehead-N purine-like ring -- the one pre-existing chematic_core::kekulize()
# failure K1 did NOT fix (see rdkit_morgan_ecfp4.rs's
# `kekule_bridgehead_n_purine_reports_kekulization_failed_not_a_fallback_result` test
# and chematic_perception::rdkit_parity's own
# `production_api_reports_kekulize_failure_not_panic`). Real error-path fixture -- not
# a fabricated/synthetic failure.
ERROR_FIXTURE = ("34", "Cc1cn2c(=O)c3ncn(COCCO)c3nc2n1C")


def rd_parse(smi):
    return Chem.MolFromSmiles(smi, _PARSE_PARAMS)


def make_gen(radius, fp_size=2048):
    return rdFingerprintGenerator.GetMorganGenerator(
        radius=radius,
        fpSize=fp_size,
        includeChirality=False,
        useBondTypes=True,
        includeRedundantEnvironments=False,
    )


def rd_sparse_counts(mol, radius):
    gen = make_gen(radius=radius)
    fp = gen.GetSparseCountFingerprint(mol)
    return {rid & 0xFFFFFFFF: c for rid, c in fp.GetNonzeroElements().items()}


def rd_raw_bit_info(mol, radius):
    """raw_id -> [(atom_idx, radius), ...] across radii 0..=radius (RDKit's real
    BitInfoMap on the unfolded fingerprint)."""
    gen = make_gen(radius=radius)
    ao = rdFingerprintGenerator.AdditionalOutput()
    ao.AllocateBitInfoMap()
    gen.GetSparseFingerprint(mol, additionalOutput=ao)
    out = {}
    for raw_id, envs in ao.GetBitInfoMap().items():
        rid = raw_id & 0xFFFFFFFF
        out[rid] = sorted((int(a), int(r)) for a, r in envs)
    return out


def rd_folded_bit_info(mol, radius, fp_size):
    gen = make_gen(radius=radius, fp_size=fp_size)
    ao = rdFingerprintGenerator.AdditionalOutput()
    ao.AllocateBitInfoMap()
    gen.GetFingerprint(mol, additionalOutput=ao)
    out = {}
    for bit, envs in ao.GetBitInfoMap().items():
        out[int(bit)] = sorted((int(a), int(r)) for a, r in envs)
    return out


def rd_folded_on_bits(mol, radius, fp_size):
    return sorted(make_gen(radius=radius, fp_size=fp_size).GetFingerprint(mol).GetOnBits())


def build_success_fixture(fid, smi):
    mol = rd_parse(smi)
    if mol is None:
        raise SystemExit(f"RDKit failed to parse fixture {fid}: {smi}")

    sparse_counts = rd_sparse_counts(mol, PRODUCTION_RADIUS)
    raw_bit_info = rd_raw_bit_info(mol, PRODUCTION_RADIUS)
    folded_bit_info = rd_folded_bit_info(mol, PRODUCTION_RADIUS, PRODUCTION_FP_SIZE)
    folded_bits = rd_folded_on_bits(mol, PRODUCTION_RADIUS, PRODUCTION_FP_SIZE)

    radius_axis = []
    for r in RADIUS_SWEEP:
        radius_axis.append({
            "radius": r,
            "folded_bits": rd_folded_on_bits(mol, r, PRODUCTION_FP_SIZE),
        })

    fp_size_axis = []
    for n in FP_SIZE_SWEEP:
        fp_size_axis.append({
            "fp_size": n,
            "folded_bits": rd_folded_on_bits(mol, PRODUCTION_RADIUS, n),
        })

    return {
        "id": fid,
        "smiles": smi,
        "expect": "ok",
        "sparse_counts": {str(k): v for k, v in sorted(sparse_counts.items())},
        "raw_bit_info": {str(k): v for k, v in sorted(raw_bit_info.items())},
        "folded_bit_info": {str(k): v for k, v in sorted(folded_bit_info.items())},
        "folded_bits": folded_bits,
        "radius_axis": radius_axis,
        "fp_size_axis": fp_size_axis,
    }


def build_error_fixture(fid, smi):
    return {
        "id": fid,
        "smiles": smi,
        "expect": "error",
        "error_kind": "Aromaticity",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    fixtures = [build_success_fixture(fid, smi) for fid, smi in SUCCESS_FIXTURES]
    fixtures.append(build_error_fixture(*ERROR_FIXTURE))

    try:
        chematic_commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip()
    except Exception:
        chematic_commit = "unknown"

    doc = {
        "schema_version": "1",
        "description": (
            "Shared cross-language fixture+expectation corpus for "
            "chematic_fp::rdkit_morgan_ecfp4_experimental / rdkit_morgan_fingerprint. "
            "Success-case values are generated directly from a live RDKit oracle, not "
            "from chematic -- Rust/Python/WASM tests all read this one file and must "
            "match it exactly."
        ),
        "rdkit_version": rdBase.rdkitVersion,
        "rdkit_pinned_commit": RDKIT_PINNED_COMMIT,
        "generator_script": "scripts/gen_ecfp4_rdkit_stable_api_fixtures.py",
        "generated_at_chematic_commit": chematic_commit,
        "production_config": {"radius": PRODUCTION_RADIUS, "fp_size": PRODUCTION_FP_SIZE},
        "radius_sweep": RADIUS_SWEEP,
        "fp_size_sweep": FP_SIZE_SWEEP,
        "fixtures": fixtures,
    }

    with open(args.out, "w") as f:
        json.dump(doc, f, indent=2, sort_keys=False)
    print(f"wrote {args.out} ({len(fixtures)} fixtures, "
          f"rdkit=={rdBase.rdkitVersion}, chematic@{chematic_commit})")


if __name__ == "__main__":
    main()
