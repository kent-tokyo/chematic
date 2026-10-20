#!/usr/bin/env python3
"""Pinned, same-graph atom permutations and RDKit #9629 reproduction."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import random
import sys

import chematic
from chematic import rdkit_compat as Compat
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors as D

ROOT = Path(__file__).resolve().parents[2]
CASES = [
    ('bicyclic-issue9629', 'C1CCN2CCCC2C1'),
    ('bicyclic-small', 'C1CCN2CCC2C1'),
    ('bridged-amine', 'C1CN2CCC1CC2'),
    ('bridged-amine-small', 'C1CC2CCC1N2'),
    ('tetrahedral-cw', 'F[C@](Cl)(Br)I'),
    ('tetrahedral-ccw', 'F[C@@](Cl)(Br)I'),
    ('mapped-implicit-h', 'N[C@@H:7](C)C(=O)O'),
    ('explicit-h', '[H][C@](F)(Cl)Br'),
    ('deuterium', '[2H][C@](F)(Cl)Br'),
    ('alkene-e', 'F/C=C/Cl'),
    ('alkene-z', 'F/C=C\\Cl'),
    ('axial-allene', 'CC=[C@AL1]=C(C)F'),
    ('square-planar', 'Cl[Pt@SP1](Cl)(N)N'),
    ('trigonal-bipyramidal', 'F[P@TB1](Cl)(Br)(I)N'),
    ('octahedral', 'F[Co@OH1](Cl)(Br)(I)(N)O'),
]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def orders(n):
    if n <= 5:
        return list(itertools.permutations(range(n)))
    result = [tuple(range(n)), tuple(reversed(range(n)))]
    result += [tuple(range(i,n)) + tuple(range(i)) for i in range(1,n)]
    rng = random.Random(9629)
    for _ in range(16):
        order = list(range(n))
        rng.shuffle(order)
        result.append(tuple(order))
    return list(dict.fromkeys(result))


def evaluate():
    if rdBase.rdkitVersion != '2026.09.1':
        raise RuntimeError(f'expected RDKit 2026.09.1, got {rdBase.rdkitVersion}')
    params = Chem.SmilesParserParams()
    params.removeHs = False
    fixtures, executed = [], []
    for identifier, smiles in CASES:
        rd_source = Chem.MolFromSmiles(smiles, params)
        assert rd_source is not None
        reference_count = D.CalcNumAtomStereoCenters(rd_source)
        expected_smiles = Chem.MolToSmiles(rd_source, isomericSmiles=True)
        source = chematic.from_smiles(smiles)
        original_rows = source.atom_table
        assert source.atom_count == rd_source.GetNumAtoms()
        assert source.num_stereocenters == reference_count
        permutations = orders(source.atom_count)
        upstream_drift = 0
        identity_record = None
        for order in permutations:
            raw = Chem.RenumberAtoms(rd_source, order)
            raw_count = D.CalcNumAtomStereoCenters(raw)
            repaired = Chem.RenumberAtoms(rd_source, order)
            Chem.GetSSSR(repaired)
            repaired_count = D.CalcNumAtomStereoCenters(repaired)
            out = source.renumber_atoms(order)
            assert out.num_stereocenters == repaired_count == reference_count, (identifier,order)
            assert out.atom_table == [original_rows[i] for i in order], (identifier,order)
            crossread = Chem.MolFromSmiles(out.smiles, params)
            assert crossread is not None
            assert Chem.MolToSmiles(crossread,isomericSmiles=True) == expected_smiles, (identifier,order,out.smiles)
            upstream_drift += raw_count != reference_count
            result = {'fixture':identifier, 'new_to_old':list(order), 'reference_count':reference_count, 'rdkit_unrepaired_count':raw_count, 'rdkit_sssr_repaired_count':repaired_count, 'chematic_count':out.num_stereocenters, 'stereo_crossread_equal':True}
            executed.append(result)
            if order == tuple(range(source.atom_count)):
                identity_record = result
        assert identity_record is not None
        if identifier == 'bicyclic-issue9629':
            assert identity_record['rdkit_unrepaired_count'] == 2
            assert identity_record['chematic_count'] == identity_record['rdkit_sssr_repaired_count'] == 1
        fixtures.append({
            'id':identifier, 'smiles':smiles, 'sha256':digest((smiles+'\n').encode()),
            'provenance':{'source':'minimal abstract molecular graphs and independently written permutation recipe', 'issue':'https://github.com/rdkit/rdkit/issues/9629', 'retrieved':'2026-10-14', 'transformations':'one minimal issue graph; synthetic ring, atom-map, isotope and stereo controls; no upstream dataset or source code copied'},
            'license':'MIT OR Apache-2.0',
            'redistribution':{'status':'allowed','rationale':'minimal abstract structures and original project test recipe'},
            'minimization':'small graphs isolate ring-dependent stereo counts or one stereo representation at a time',
            'oracle':{'engine':'RDKit','version':'2026.09.1','settings':{'legacy_stereo':True,'removeHs':False,'count_reference':'original MolFromSmiles and RenumberAtoms followed by GetSSSR','configuration_reference':'RDKit isomeric canonical SMILES after cross-read'}},
            'expected':{'stereocenter_count':reference_count,'atom_count':source.atom_count,'configuration':'same isomeric canonical structure','count_tolerance':0},
            'permutations':len(permutations),'rdkit_unrepaired_count_drifts':upstream_drift,
            'native_count_mismatches':0,'native_configuration_mismatches':0,
            'identity':identity_record,
        })
    module=Path(chematic.chematic.__file__).resolve()
    source_paths=['crates/chematic-core/src/molecule/renumber.rs','crates/chematic-py/src/mol_methods/transformations.rs','crates/chematic-py/python/chematic/rdkit_compat.py']
    return {
        'schema_version':1,'status':'local-source-candidate',
        'base':'cac480fb23cd8adeec6e5e9bed1f18677f87b392 source snapshot plus preserved RDKit #9675 and Open Babel #2845 candidates; patch isolates this feature',
        'scope':'15 diagnostic structures; no broad stereoperception or docking claim; external coordinate arrays are caller-owned',
        'runtime':{'python':sys.version.split()[0],'rdkit':rdBase.rdkitVersion,'chematic_extension_sha256':digest(module.read_bytes())},
        'permutation_recipe':{'small_n':'all permutations when n<=5','larger_n':'identity, reverse, every cyclic shift, 16 seeded random shuffles, deduplicated','seed':9629},
        'source_sha256':{s:digest((ROOT/s).read_bytes()) for s in source_paths},
        'results':{'fixtures':len(fixtures),'permutations':len(executed),'native_count_mismatches':0,'native_configuration_mismatches':0,'rdkit_unrepaired_count_drifts':sum(f['rdkit_unrepaired_count_drifts'] for f in fixtures),'cases_sha256':digest(json.dumps(executed,sort_keys=True,separators=(',',':')).encode())},
        'fixtures':fixtures,
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'validation/results/rdkit-atom-renumbering-2026-10-14.json')
    args=parser.parse_args()
    previous=Chem.GetUseLegacyStereoPerception()
    try:
        Chem.SetUseLegacyStereoPerception(True)
        record=evaluate()
    finally:
        Chem.SetUseLegacyStereoPerception(previous)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record['results']))


if __name__=='__main__':
    main()
