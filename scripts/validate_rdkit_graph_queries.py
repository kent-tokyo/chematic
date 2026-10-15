#!/usr/bin/env python3
"""Exact graph matrices and tie-sensitive paths against pinned RDKit."""
import argparse
import hashlib
import json
from pathlib import Path
import random

import chematic
from rdkit import Chem, rdBase, RDLogger

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT/'validation/rdkit_issues/graph-queries/fixtures.json'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def orders(n):
    result = [list(range(n))]
    reverse = list(reversed(range(n)))
    if reverse not in result:
        result.append(reverse)
    shuffled = list(range(n))
    random.Random(9707+n).shuffle(shuffled)
    if shuffled not in result:
        result.append(shuffled)
    return result


def variant(case, added_hs, order):
    native = chematic.from_smiles(case['smiles']) if case['smiles'] else chematic.RWMol().GetMol()
    if added_hs:
        native = native.rdkit_add_hydrogens()
        reference = Chem.AddHs(Chem.MolFromSmiles(case['smiles']))
    else:
        params = Chem.SmilesParserParams()
        params.removeHs = False
        reference = Chem.MolFromSmiles(case['smiles'], params)
    assert reference is not None
    assert native.atom_count == reference.GetNumAtoms()
    if order is not None:
        native = native.renumber_atoms(order)
        reference = Chem.RenumberAtoms(reference, order) if order else Chem.Mol(reference)
    if not added_hs:
        reference = Chem.RemoveHs(reference)
    return native, reference


def compare(native, reference, checksum):
    n = reference.GetNumAtoms()
    assert native.rdkit_num_atoms() == n
    original = native.smiles
    native.SetIntProp('graph_query_sentinel',9707)
    entries = 0
    for use_bo in (False,True,False):
        want = Chem.GetAdjacencyMatrix(reference,useBO=use_bo,emptyVal=0,force=True).tolist()
        got = native.rdkit_adjacency_matrix(use_bo=use_bo)
        assert got == want, (use_bo,got,want)
        assert len(got) == n and all(len(row)==n for row in got)
        checksum.update(json.dumps(got,separators=(',',':')).encode())
        entries += n*n
    paths = disconnected = tied = 0
    for start in range(n):
        for end in range(n):
            if start == end:
                continue
            want = list(Chem.GetShortestPath(reference,start,end))
            got = native.rdkit_shortest_path(start,end)
            assert got == want, (start,end,got,want)
            paths += 1
            disconnected += not bool(got)
            if got:
                opposite = list(Chem.GetShortestPath(reference,end,start))
                tied += got != list(reversed(opposite))
                assert len(got)-1 == int(Chem.GetDistanceMatrix(reference)[start,end])
            checksum.update(json.dumps([start,end,got],separators=(',',':')).encode())
    assert native.smiles == original
    assert native.GetIntProp('graph_query_sentinel') == 9707
    return entries, paths, disconnected, tied


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'validation/results/rdkit-graph-queries-2026-10-15.json')
    args = parser.parse_args()
    assert rdBase.rdkitVersion == '2026.09.1', rdBase.rdkitVersion
    RDLogger.DisableLog('rdApp.*')
    manifest = json.loads(FIXTURES.read_text())
    assert manifest['sha256'] == digest(json.dumps(manifest['cases'],sort_keys=True,separators=(',',':')).encode())
    checksum = hashlib.sha256()
    variants = entries = paths = disconnected = tied = 0
    records = []
    for case in manifest['cases']:
        assert case['sha256'] == digest(case['smiles'].encode())
        case_variants = 0
        for added_hs in (False,True):
            native, reference = variant(case,added_hs,None)
            for order in orders(native.atom_count):
                native, reference = variant(case,added_hs,order)
                result = compare(native,reference,checksum)
                entries += result[0]; paths += result[1]; disconnected += result[2]; tied += result[3]
                variants += 1; case_variants += 1
        records.append({'id':case['id'],'status':'passed','variants':case_variants})
    # Separate diagnostics: nonzero fills are outside the compatibility profile.
    probe = Chem.MolFromSmiles('CC.C')
    byte_fill = float(Chem.GetAdjacencyMatrix(probe,useBO=True,emptyVal=1,force=True)[0,0])
    assert byte_fill != 1.0 and byte_fill > 0
    cache_probe = Chem.MolFromSmiles('CC.C')
    baseline = Chem.GetAdjacencyMatrix(cache_probe,useBO=True,emptyVal=0,force=True)
    cached = Chem.GetAdjacencyMatrix(cache_probe,useBO=True,emptyVal=1,force=False)
    assert baseline.tolist() == cached.tolist()
    report = {
        'schema_version':1, 'status':'local-source-candidate; uncommitted and unpushed',
        'oracle':{'name':'RDKit','version':rdBase.rdkitVersion,
            'settings':'MolFromSmiles with removeHs=False then RenumberAtoms/RemoveHs, or AddHs then RenumberAtoms; GetAdjacencyMatrix(force=True,emptyVal=0), GetShortestPath'},
        'fixture_manifest':str(FIXTURES.relative_to(ROOT)), 'fixture_sha256':manifest['sha256'],
        'results':{'fixtures':len(records),'graph_variants':variants,'matrix_calls':3*variants,
            'exact_matrix_entries':entries,'exact_shortest_paths':paths,
            'disconnected_paths':disconnected,'direction_dependent_tied_paths':tied,
            'unexpected_mismatches':0,'normalized_output_sha256':checksum.hexdigest()},
        'records':records,
        'observed_rdkit_nonzero_fill':{'weighted_emptyVal_1_diagonal':byte_fill,
            'explanation':'Release_2026_09_1 Matrices.cpp fills double storage using memset(emptyVal); nonzero fill is refused by chematic',
            'cache_probe':'Changing emptyVal without force reuses existing adjacency in this oracle; chematic does not cache graph query outputs'},
        'scope':['Rust/Python only','Sanitized RDKit model atom order','Input cap 1024 atoms / 16384 bonds',
            'Nonzero empty_value refused','Query/zero bonds and unsupported model inputs refused','No WASM or rdkit_compat facade functions'],
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['results'],indent=2))


if __name__ == '__main__':
    main()
