#!/usr/bin/env python3
"""Compare H-folded native queries to unmerged RDKit queries on explicit-H targets."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import platform
import subprocess

import chematic
from rdkit import Chem, rdBase

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert rdBase.rdkitVersion == "2026.09.1", rdBase.rdkitVersion
    targets = ["C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O",
               "[H][C@@](C)(N)C(=O)O", "[H][C@](C)(N)C(=O)O"]
    patterns = []
    for neighbors in itertools.permutations(["[H]", "C", "N", "C(=O)O"]):
        for tag in ["@", "@@"]:
            patterns.append(f"[C{tag}]({neighbors[0]})({neighbors[1]})({neighbors[2]}){neighbors[3]}")
    cases = [(p, t) for p in patterns for t in targets]
    for p, target in [
        ("[C@]12([H])CC(C)CC1NC2", "[C@H]12CC(C)CC1NC2"),
        ("[C@]12([H])CC(C)CC2NC1", "[C@H]12CC(C)CC2NC1"),
        ("[C@]1([H])2CC(C)CC1NC2", "[C@@H]12CC(C)CC1NC2"),
    ]:
        cases.append((p, target))
        cases.append((p, target.replace("@@", "TEMP").replace("@", "@@").replace("TEMP", "@")))
    mismatches, upstream_differences = [], []
    for pattern, target in cases:
        rd_target = Chem.AddHs(Chem.MolFromSmiles(target))
        expected = rd_target.HasSubstructMatch(Chem.MolFromSmarts(pattern), useChirality=True)
        query = chematic.SmartsQuery(pattern)
        folded = query.merge_hydrogens()
        assert folded.atom_count == query.atom_count - 1
        actual = folded.matches(chematic.from_smiles(target), use_chirality=True)
        upstream = Chem.MolFromSmiles(target).HasSubstructMatch(
            Chem.MolFromSmarts(pattern, mergeHs=True), useChirality=True)
        record = dict(pattern=pattern, target=target, unmerged_explicit_h_reference=expected,
                      chematic_folded=actual, rdkit_folded=upstream)
        if actual != expected:
            mismatches.append(record)
        if upstream != expected:
            upstream_differences.append(record)
    assert not mismatches, mismatches
    assert upstream_differences
    paths = ["crates/chematic-smarts/src/merge_hydrogens.rs",
             "crates/chematic-smarts/src/match_vf2.rs",
             "crates/chematic-smarts/src/match_vf2/chirality.rs",
             "crates/chematic-smarts/src/parser.rs",
             "crates/chematic-smarts/src/parser/rings.rs",
             "crates/chematic-py/src/misc.rs", "scripts/validation/validate_query_hydrogens.py"]
    module = Path(chematic.__file__).resolve()
    report = dict(issue="https://github.com/rdkit/rdkit/issues/9675",
                  source_baseline=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                  artifact_lane="uncommitted source extension; not a release", python=platform.python_version(),
                  rdkit_version=rdBase.rdkitVersion, tested_cases=len(cases),
                  fixture_provenance=dict(id="query-hydrogens-generated-v1", retrieved="2026-10-14",
                                          license="MIT OR Apache-2.0", redistribution="allowed",
                                          rationale="Project-authored minimal chemical queries and generated neighbour permutations; no upstream implementation or bulk data copied",
                                          minimization="One amino-acid stereocentre and three ring-order examples",
                                          oracle="RDKit 2026.09.1, unmerged queries, Chem.AddHs targets, useChirality=True"),
                  cases_sha256=hashlib.sha256(json.dumps(cases, separators=(",", ":")).encode()).hexdigest(),
                  pattern_generation="24 root-neighbour permutations times both tags, four target spellings; three ring patterns times both configurations",
                  chematic_mismatches=mismatches, rdkit_folded_mismatch_count=len(upstream_differences),
                  rdkit_folded_mismatch_examples=upstream_differences[:12],
                  module_sha256=hashlib.sha256(module.read_bytes()).hexdigest(),
                  source_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                  scope="plain unmapped terminal H on single/implicit bonds, complete tetrahedral affected frames; recursive queries rejected")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(dict(tested_cases=len(cases), chematic_mismatches=len(mismatches),
                          rdkit_folded_mismatches=len(upstream_differences))))


if __name__ == "__main__":
    main()
