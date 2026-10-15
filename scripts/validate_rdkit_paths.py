#!/usr/bin/env python3
"""Ordered path vectors against pinned RDKit 2026.09.1, original fixtures."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import platform
import random

import chematic
from rdkit import Chem, RDLogger, rdBase

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT/'validation/rdkit_issues/path-enumeration/fixtures.json'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def orders(n):
    out = [list(range(n))]
    for row in (list(reversed(range(n))), random.Random(9710+n).sample(range(n), n)):
        if row not in out:
            out.append(row)
    return out


def variant(case, added_hs, order):
    native = chematic.from_smiles(case['smiles']) if case['smiles'] else chematic.RWMol().GetMol()
    if added_hs:
        native = native.rdkit_add_hydrogens()
        reference = Chem.AddHs(Chem.MolFromSmiles(case['smiles']))
    else:
        params = Chem.SmilesParserParams()
        params.removeHs = False
        reference = Chem.MolFromSmiles(case['smiles'], params)
    assert reference is not None and native.atom_count == reference.GetNumAtoms()
    native = native.renumber_atoms(order)
    reference = Chem.RenumberAtoms(reference, order) if order else Chem.Mol(reference)
    if not added_hs:
        reference = Chem.RemoveHs(reference)
    return native, reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'validation/results/rdkit-path-enumeration-2026-10-15.json')
    args = parser.parse_args()
    assert rdBase.rdkitVersion == '2026.09.1', rdBase.rdkitVersion
    RDLogger.DisableLog('rdApp.*')
    manifest = json.loads(FIXTURE.read_text())
    assert manifest['sha256'] == digest(json.dumps(manifest['cases'], sort_keys=True, separators=(',', ':')).encode())
    for case in manifest['cases']:
        assert case['sha256'] == digest(case['smiles'].encode())
    for case in manifest['golden']:
        assert chematic.from_smiles(case['smiles']).rdkit_paths_of_length(case['length'], **case['params']) == case['expected']['paths']
    checksum = hashlib.sha256()
    calls = variants = paths = indices = empty = shortest = hs = permuted = closures = 0
    records = []
    for case in manifest['cases']:
        checked = 0
        for added_hs in (False, True):
            base = chematic.from_smiles(case['smiles']) if case['smiles'] else chematic.RWMol().GetMol()
            if added_hs:
                base = base.rdkit_add_hydrogens()
            for order in orders(base.atom_count):
                native, reference = variant(case, added_hs, order)
                n = reference.GetNumAtoms()
                assert native.rdkit_num_atoms() == n
                assert native.rdkit_adjacency_matrix() == Chem.GetAdjacencyMatrix(reference, force=True).tolist()
                before = native.smiles
                native.SetIntProp('path_sentinel', 9710)
                lengths = list(dict.fromkeys([*range(min(n+3, 8)), n+1, n+2]))
                roots = [None, *range(n), n, n+7]
                for length, use_bonds, use_hs, only_shortest, root in itertools.product(lengths, (False, True), (False, True), (False, True), roots):
                    want = [list(row) for row in Chem.FindAllPathsOfLengthN(reference, length, useBonds=use_bonds, useHs=use_hs, rootedAtAtom=-1 if root is None else root, onlyShortestPaths=only_shortest)]
                    got = native.rdkit_paths_of_length(length, use_bonds=use_bonds, use_hs=use_hs, rooted_at_atom=root, only_shortest_paths=only_shortest)
                    assert got == want, (case['id'], added_hs, order, length, use_bonds, use_hs, only_shortest, root, got, want)
                    checksum.update(json.dumps([case['id'], added_hs, order, length, use_bonds, use_hs, only_shortest, root, want], separators=(',', ':')).encode())
                    calls += 1
                    checked += 1
                    paths += len(want)
                    indices += sum(map(len, want))
                    empty += not bool(want)
                    shortest += only_shortest
                    hs += use_hs
                    permuted += order != list(range(len(order)))
                    if not use_bonds:
                        closures += sum(len(set(row)) != len(row) for row in want)
                assert native.smiles == before and native.GetIntProp('path_sentinel') == 9710
                variants += 1
        records.append(dict(id=case['id'], status='pass', exact_comparisons=checked))
    # Huge lengths still iterate in upstream after its frontier is empty, and
    # UINT32_MAX can overflow counters. Test these only in the bounded native API.
    for use_bonds in (False, True):
        for length in (2**32-2, 2**32-1):
            assert chematic.from_smiles('C1CC1').rdkit_paths_of_length(length, use_bonds=use_bonds, max_visits=0) == []
    report = dict(schema_version=1, date='2026-10-15', status='isolated source candidate; uncommitted and unpushed',
                  runtime=dict(python=platform.python_version(), rdkit=rdBase.rdkitVersion),
                  oracle=dict(name='RDKit', version=rdBase.rdkitVersion, release='Release_2026_09_1', commit='fece8caa860bdf6c9c82bdb9253f22a0cf87150c', source='Code/GraphMol/Subgraphs/Subgraphs.cpp', source_sha256='aa8295ce04042f809d93a29f8f6ca0beca2c656b4f3865270237131119a55597', settings='Exact ordered vectors, all roots plus unrooted/outside roots; bond/atom, H and shortest flags; H removal/addition; three distinct seeded atom permutations where available; lengths 0..min(n+2,7), n+1, n+2'),
                  fixture_manifest=str(FIXTURE.relative_to(ROOT)), fixture_manifest_sha256=digest(FIXTURE.read_bytes()),
                  results=dict(fixtures=len(records), golden_cases=len(manifest['golden']), graph_variants=variants, exact_comparisons=calls, paths_compared=paths, indices_compared=indices, empty_results=empty, shortest_mode_cases=shortest, use_hs_true_cases=hs, permuted_cases=permuted, returned_atom_ring_closures=closures, unexpected_mismatches=0, normalized_output_sha256=checksum.hexdigest()), records=records,
                  scope=['Sanitized RDKit model indices; input <=1024 atoms/16384 bonds', 'Rust/Python only, no branched or unique chemical subgraph API', 'Bond-set deduplication with first path retained, including atom output', 'Shortest mode includes graph Hs regardless of useHs; single-atom output includes isolated Hs', 'Source molecule/properties immutable, independent returned vectors', 'Bounded frontiers/work/results fail without partial output', 'UINT32_MAX tested only natively, not against upstream loop overflow; no timing or library-wide parity claim'])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['results'], indent=2))


if __name__ == '__main__':
    main()
