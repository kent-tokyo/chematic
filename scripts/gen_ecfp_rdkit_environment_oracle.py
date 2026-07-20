#!/usr/bin/env python3
"""
RDKit Morgan/ECFP environment-parity oracle generator for the M2 diagnostic
PR (see scripts/ecfp_rdkit_environment_parity.py for the comparator that
consumes this).

For every input SMILES, dumps RDKit's own Morgan-fingerprint intermediate
state at the granularity the comparator needs: per-atom connectivity
invariant (radius 0), unfolded sparse (atom,radius)->raw-id attribution at
every radius (via the *unfolded* GetSparseFingerprint's AdditionalOutput --
NOT the folded GetFingerprint's, which has 2048-modulo collisions), the
folded 2048-bit fingerprint + its own bitInfo, and sparse counts.

Two Morgan-generator variants are captured per molecule:
  - "default": includeRedundantEnvironments=False -- RDKit's real gate/
    default behavior (radius=2, fpSize=2048, includeChirality=False,
    useBondTypes=True, onlyNonzeroInvariants=False,
    includeRingMembership=True [RDKit's own generator default, recorded
    explicitly rather than left implicit]).
  - "full": same options but includeRedundantEnvironments=True -- chematic's
    own semantics (chematic never suppresses redundant environments today),
    so this is the correct comparison base for the radius-1/radius-2
    identifier-mismatch stages. The difference between "full" and "default"
    bitInfo pairs is RDKit's OWN actual suppression decision -- ground truth
    for the redundant-environment-mismatch stage, not inferred from ball
    stability.

Primary-source confirmation (verified interactively against RDKit 2026.03.3,
not re-derived here): a Morgan generator built with radius=2, fpSize=2048,
includeChirality=False, useBondTypes=True, onlyNonzeroInvariants=False,
includeRedundantEnvironments=False, folded via GetFingerprint(), produces a
bit-identical on-bit set to legacy AllChem.GetMorganFingerprintAsBitVect(mol,
2, 2048) -- confirming the modern generator API with these exact options is
semantically the same API this project's existing reference tooling
(scripts/ecfp_rdkit_invariants_fingerprint_ref.py) already targets.

Usage:
    .venv/bin/python scripts/gen_ecfp_rdkit_environment_oracle.py \\
        --corpus ~/Downloads/SMILES.csv \\
        --fixtures scripts/ecfp_rdkit_edge_fixtures.csv scripts/ecfp_rdkit_morgan_edge_fixtures.csv \\
        --rows-out /path/to/oracle_rows.jsonl \\
        --manifest-out validation/ecfp_rdkit_environment_parity_manifest.json

    .venv/bin/python scripts/gen_ecfp_rdkit_environment_oracle.py --verify-determinism \\
        --corpus scripts/ecfp_rdkit_morgan_edge_fixtures.csv
"""

import argparse
import datetime
import hashlib
import json
import os
import subprocess
import sys

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

SCHEMA_VERSION = "1"

MORGAN_OPTIONS = {
    "radius": 2,
    "fpSize": 2048,
    "countSimulation": False,
    "includeChirality": False,
    "useBondTypes": True,
    "onlyNonzeroInvariants": False,
    "includeRingMembership": True,  # RDKit generator's own real default -- recorded, not implicit
}

_PARSE_PARAMS = Chem.SmilesParserParams()
_PARSE_PARAMS.removeHs = False


def parse_smiles(smi):
    """Matches scripts/ecfp_rdkit_invariant_parity.py's parsing: chematic
    keeps explicit H atoms as real graph atoms, so RDKit's default (which
    silently strips 'trivial' explicit H) must be disabled to keep atom_idx
    correspondence comparable."""
    return Chem.MolFromSmiles(smi, _PARSE_PARAMS)


def make_generator(include_redundant):
    return rdFingerprintGenerator.GetMorganGenerator(
        radius=MORGAN_OPTIONS["radius"],
        fpSize=MORGAN_OPTIONS["fpSize"],
        countSimulation=MORGAN_OPTIONS["countSimulation"],
        includeChirality=MORGAN_OPTIONS["includeChirality"],
        useBondTypes=MORGAN_OPTIONS["useBondTypes"],
        onlyNonzeroInvariants=MORGAN_OPTIONS["onlyNonzeroInvariants"],
        includeRingMembership=MORGAN_OPTIONS["includeRingMembership"],
        includeRedundantEnvironments=include_redundant,
    )


def _bit_info_dict(ao):
    raw = {str(k): [list(pair) for pair in v] for k, v in ao.GetBitInfoMap().items()}
    return dict(sorted(raw.items(), key=lambda kv: int(kv[0])))


def sparse_dump(gen, mol):
    """True unfolded raw ids + their (atom,radius) attribution -- no
    2048-modulo fold collision, the correct source for identifier/membership
    comparison.

    RDKit quirk (confirmed empirically, not documented): SparseBitVect's
    GetOnBits() returns each raw 32-bit hash as a SIGNED int, while the same
    generator's AdditionalOutput bitInfoMap keys the identical hash as
    UNSIGNED (e.g. -2048234256 vs 2246733040 for the same environment --
    2246733040 - 2**32 == -2048234256). Normalize on-bits to unsigned here so
    "sparse_on_bits" and "sparse_bit_info" always key the same representation
    and the comparator never has to special-case this.
    """
    ao = rdFingerprintGenerator.AdditionalOutput()
    ao.AllocateBitInfoMap()
    sparse_fp = gen.GetSparseFingerprint(mol, additionalOutput=ao)
    on_bits = sorted(bit & 0xFFFFFFFF for bit in sparse_fp.GetOnBits())
    return {
        "sparse_on_bits": on_bits,
        "sparse_bit_info": _bit_info_dict(ao),
    }


def default_dump(gen, mol):
    """sparse_dump plus the real folded 2048-bit gate-shape outputs."""
    out = sparse_dump(gen, mol)
    ao = rdFingerprintGenerator.AdditionalOutput()
    ao.AllocateBitInfoMap()
    folded_fp = gen.GetFingerprint(mol, additionalOutput=ao)
    count_fp = gen.GetSparseCountFingerprint(mol)
    out["folded_on_bits"] = sorted(folded_fp.GetOnBits())
    out["folded_bit_info"] = _bit_info_dict(ao)
    out["sparse_counts"] = {str(k): v for k, v in sorted(count_fp.GetNonzeroElements().items())}
    return out


def atom_ball(mol, radius, atom_idx):
    """Sorted atom indices within `radius` bond-hops of `atom_idx`, derived
    from RDKit's own Chem.FindAtomEnvironmentOfRadiusN (bond indices -> atom
    endpoints) -- the hash-independent environment-membership ground truth,
    RDKit's counterpart to chematic's own graph-BFS atom_ball
    (ecfp_diagnostics.rs).

    Two non-default arguments are required to make this a genuine CUMULATIVE
    ball comparable to chematic's BFS (confirmed empirically, not documented
    behavior at a glance -- FindAtomEnvironmentOfRadiusN's defaults do
    something different from what its docstring's one-liner suggests):
      - enforceSize=False: the *default* (True) returns an EMPTY bond list
        once the environment stops growing at exactly this radius (RDKit's
        internal "nothing new, treat as redundant" signal) -- e.g. ethane's
        atom 0 at radius=2 returns [] under the default, even though its
        real (saturated) 2-hop ball is still {0, 1}. enforceSize=False makes
        it return every bond within `radius` hops, monotonically growing,
        matching chematic's own BFS semantics.
      - useHs=True: the *default* (False) silently excludes bonds to
        explicit H atoms from the result -- e.g. plain methane's carbon
        (`[H]C([H])([H])[H]`, all 4 H explicit) returns [] at radius=1 under
        the default, even though the carbon has 4 real bonds. chematic keeps
        explicit H atoms as real graph atoms with no special-casing, so the
        RDKit side must include them too for the comparison to be apples-to-
        apples.
    """
    atoms = {atom_idx}
    for bond_idx in Chem.FindAtomEnvironmentOfRadiusN(
        mol, radius, atom_idx, useHs=True, enforceSize=False
    ):
        bond = mol.GetBondWithIdx(bond_idx)
        atoms.add(bond.GetBeginAtomIdx())
        atoms.add(bond.GetEndAtomIdx())
    return sorted(atoms)


def atom_balls(mol, radius):
    return {str(a.GetIdx()): atom_ball(mol, radius, a.GetIdx()) for a in mol.GetAtoms()}


def oracle_row(row_id, source, smi):
    rd = parse_smiles(smi)
    if rd is None:
        return {"row_id": row_id, "smiles": smi, "source": source, "parse_ok": False}

    gen_default = make_generator(include_redundant=False)
    gen_full = make_generator(include_redundant=True)

    return {
        "row_id": row_id,
        "smiles": smi,
        "source": source,
        "parse_ok": True,
        "atom_count": rd.GetNumAtoms(),
        "atomic_numbers": [a.GetAtomicNum() for a in rd.GetAtoms()],
        "connectivity_invariants": list(rdMolDescriptors.GetConnectivityInvariants(rd, True)),
        "atom_balls": {"1": atom_balls(rd, 1), "2": atom_balls(rd, 2)},
        "default": default_dump(gen_default, rd),
        "full": sparse_dump(gen_full, rd),
    }


def load_sources(corpus_paths, fixture_paths):
    """[(source_label, smiles), ...] -- duplicates across files are kept as
    separate rows (the comparator's own accounting is index-based, not
    dedup'd), source_label distinguishes corpus vs which fixture file."""
    sources = []
    for p in corpus_paths:
        with open(os.path.expanduser(p)) as f:
            for line in f:
                smi = line.strip()
                if smi:
                    sources.append(("corpus", smi))
    for p in fixture_paths:
        label = f"fixture:{os.path.basename(p)}"
        with open(p) as f:
            for line in f:
                smi = line.strip()
                if smi:
                    sources.append((label, smi))
    return sources


def sha256_file(path):
    h = hashlib.sha256()
    with open(os.path.expanduser(path), "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def generate_rows(sources):
    return [oracle_row(i, source, smi) for i, (source, smi) in enumerate(sources)]


def write_jsonl(rows, path):
    with open(path, "w") as f:
        for row in rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")


def _repo_root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def diagnostic_source_commit_sha():
    return (
        subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=_repo_root()).decode().strip()
    )


def diagnostic_source_tree_dirty():
    """True if any TRACKED file differs from HEAD -- ignores untracked files
    on purpose, since the validation/* artifacts this script itself is about
    to produce are legitimately untracked/uncommitted at generation time
    (they get committed in a follow-up commit). What must NOT be dirty is
    the source (Rust/Python/fixtures) that produced this run's numbers."""
    diff = subprocess.check_output(["git", "diff", "--name-only", "HEAD"], cwd=_repo_root())
    return len(diff.strip()) > 0


def build_manifest(args, stamp, row_count, input_count):
    return {
        "schema_version": SCHEMA_VERSION,
        "diagnostic_source_commit_sha": diagnostic_source_commit_sha(),
        "diagnostic_source_tree_dirty": diagnostic_source_tree_dirty(),
        "oracle_generator_script_sha256": sha256_file(os.path.abspath(__file__)),
        "rdkit_version": rdBase.rdkitVersion,
        "python_version": sys.version,
        "corpus_paths": args.corpus,
        "corpus_sha256": {p: sha256_file(p) for p in args.corpus},
        "fixture_paths": args.fixtures,
        "fixture_sha256": {p: sha256_file(p) for p in args.fixtures},
        "morgan_options": MORGAN_OPTIONS,
        "morgan_options_note": (
            "GetMorganGenerator(radius=2, fpSize=2048, includeChirality=False, "
            "useBondTypes=True, onlyNonzeroInvariants=False, includeRingMembership=True, "
            "includeRedundantEnvironments=False) + GetFingerprint() confirmed bit-identical "
            "to legacy AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048) -- verified "
            "interactively against RDKit 2026.03.3, on-bit sets equal."
        ),
        "generated_at_utc": stamp,
        "command_line": sys.argv,
        "oracle_row_count": row_count,
        "input_count": input_count,
    }


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--corpus", nargs="*", default=[])
    p.add_argument("--fixtures", nargs="*", default=[])
    p.add_argument("--rows-out", default=None)
    p.add_argument("--manifest-out", default=None)
    p.add_argument(
        "--stamp", default=None, help="generated_at_utc override (default: current UTC time)"
    )
    p.add_argument("--verify-determinism", action="store_true")
    args = p.parse_args()

    sources = load_sources(args.corpus, args.fixtures)
    if not sources:
        print("no input SMILES (pass --corpus/--fixtures)", file=sys.stderr)
        sys.exit(1)

    if args.verify_determinism:
        rows_a = generate_rows(sources)
        rows_b = generate_rows(sources)
        text_a = "\n".join(json.dumps(r, sort_keys=True) for r in rows_a)
        text_b = "\n".join(json.dumps(r, sort_keys=True) for r in rows_b)
        if text_a != text_b:
            print(
                "DETERMINISM CHECK FAILED: two generations of the same input differ",
                file=sys.stderr,
            )
            sys.exit(1)
        print(f"determinism OK: {len(rows_a)} rows byte-identical across two runs")
        return

    rows = generate_rows(sources)
    parse_fail = sum(1 for r in rows if not r["parse_ok"])
    print(f"input={len(sources)} rows={len(rows)} rdkit_parse_fail={parse_fail}")

    if args.rows_out:
        write_jsonl(rows, args.rows_out)
        print(f"rows written to {args.rows_out}")

    if args.manifest_out:
        stamp = args.stamp or datetime.datetime.now(datetime.timezone.utc).isoformat()
        manifest = build_manifest(args, stamp, len(rows), len(sources))
        with open(args.manifest_out, "w") as f:
            json.dump(manifest, f, sort_keys=True, indent=2)
            f.write("\n")
        print(f"manifest written to {args.manifest_out}")


if __name__ == "__main__":
    main()
