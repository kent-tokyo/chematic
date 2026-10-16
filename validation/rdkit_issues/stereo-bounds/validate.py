#!/usr/bin/env python3
"""Check every fixture atom permutation against hand labels and RDKit."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    import chematic
    import chematic.chematic as native
    import rdkit
    from rdkit import Chem

    cases = json.loads((HERE/'cases.json').read_text())['cases']
    args.output_dir.mkdir(parents=True, exist_ok=True)
    proof = bytearray()
    decisions = permutations = 0
    summaries = []
    for case in cases:
        mol = chematic.from_smiles(case['smiles'])
        reference = Chem.MolFromSmiles(case['smiles'])
        if case['add_hs']:
            mol, reference = mol.add_hydrogens(), Chem.AddHs(reference)
        assert reference.GetNumAtoms() == mol.atom_count
        count = 0
        for order in itertools.permutations(range(mol.atom_count)):
            actual = mol.renumber_atoms(order)
            oracle = Chem.RenumberAtoms(reference, order)
            results = []
            for row in case['paths']:
                path = [order.index(i) for i in row['path']]
                bond = oracle.GetBondBetweenAtoms(path[1],path[2])
                stereo = bond.GetStereo()
                selected = set(bond.GetStereoAtoms())
                if stereo in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY):
                    expected = None
                else:
                    assert stereo in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOE)
                    # Independent unordered-set side parity using RDKit's own
                    # reference atoms, rather than native begin/end normalization.
                    same_side = stereo == Chem.BondStereo.STEREOZ
                    if (path[0] in selected) != (path[3] in selected):
                        same_side = not same_side
                    expected = 'cis' if same_side else 'trans'
                forward = actual.rdkit_alkene_path_stereo(path)
                reverse = actual.rdkit_alkene_path_stereo(path[::-1])
                assert expected == row['relation'] == forward == reverse, (case,order,path,expected,forward,reverse)
                results.append({'path':path,'rdkit_reference_atoms':sorted(selected),
                                'expected':expected,'native_forward':forward,'native_reverse':reverse})
                decisions += 2
            proof.extend((json.dumps({'smiles':case['smiles'],'add_hs':case['add_hs'],
                                     'order':order,'results':results},separators=(',',':'))+'\n').encode())
            permutations += 1
            count += 1
        summaries.append({'smiles':case['smiles'],'add_hs':case['add_hs'],'permutations':count})
    compressed = gzip.compress(proof,mtime=0)
    (args.output_dir/'permutations.jsonl.gz').write_bytes(compressed)
    report = {'schema':'chematic.alkene-path-validation.v1',
              'issue':'https://github.com/rdkit/rdkit/issues/9719',
              'scope':'Declared four-atom relative side queries; not an upstream screenshot reproduction, matrix change, embedding-quality, speed or release claim.',
              'cases':summaries,'permutations':permutations,'path_direction_decisions':decisions,'mismatches':0,
              'rdkit_version':rdkit.__version__,'python':sys.version,
              'candidate_module_sha256':hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
              'fixture_sha256':hashlib.sha256((HERE/'cases.json').read_bytes()).hexdigest(),
              'proof_sha256':hashlib.sha256(compressed).hexdigest(),
              'decoded_proof_sha256':hashlib.sha256(proof).hexdigest()}
    (args.output_dir/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({key:report[key] for key in ['permutations','path_direction_decisions','mismatches']}))


if __name__ == '__main__':
    main()
