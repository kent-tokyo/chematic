#!/usr/bin/env python3
"""Inspect declared stereo and CIP availability without generating conformers."""
import argparse
import hashlib
import json
import platform
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sdf', type=Path)
    parser.add_argument('split', type=Path)
    parser.add_argument('--partition', choices=('development', 'holdout'), required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog('rdApp.*')
    split_bytes = args.split.read_bytes()
    split = json.loads(split_bytes)
    source_hash = hashlib.sha256(args.sdf.read_bytes()).hexdigest()
    if split['source_sha256'] != source_hash:
        raise ValueError('source hash mismatch')
    supplier = Chem.SDMolSupplier(str(args.sdf), removeHs=False)
    rows, type_counts = [], Counter()
    for record in split['rows']:
        if record.get('split') != args.partition:
            continue
        row = {'index': record['index'], 'name': record['name'], 'elements': []}
        try:
            reference = supplier[record['index']]
            if reference is None or reference.GetProp('_Name') != record['name']:
                raise ValueError('reference identity mismatch')
            smiles = Chem.MolToSmiles(Chem.RemoveHs(reference), isomericSmiles=True)
            graph = Chem.MolFromSmiles(smiles)
            if graph is None:
                raise ValueError('canonical graph could not be parsed')
            mol = Chem.AddHs(graph)
            declared_atoms = {a.GetIdx() for a in mol.GetAtoms()
                              if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED}
            declared_bonds = {b.GetIdx() for b in mol.GetBonds()
                              if b.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ,
                                                   Chem.BondStereo.STEREOCIS, Chem.BondStereo.STEREOTRANS)}
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
            recognized_atoms, recognized_bonds = set(), set()
            for stereo in Chem.FindPotentialStereo(mol):
                if stereo.specified != Chem.StereoSpecified.Specified:
                    continue
                kind = str(stereo.type)
                item = {'type': kind, 'center': stereo.centeredOn}
                if kind.startswith('Atom_'):
                    recognized_atoms.add(stereo.centeredOn)
                elif kind.startswith('Bond_'):
                    recognized_bonds.add(stereo.centeredOn)
                if kind == 'Atom_Tetrahedral':
                    item['cip_available'] = mol.GetAtomWithIdx(stereo.centeredOn).HasProp('_CIPCode')
                row['elements'].append(item)
                type_counts[kind] += 1
            row['unrecognized_declared_atoms'] = sorted(declared_atoms - recognized_atoms)
            row['unrecognized_declared_bonds'] = sorted(declared_bonds - recognized_bonds)
            row['unrecognized_atom_topologies'] = []
            if row['unrecognized_declared_atoms']:
                ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=True))
                for index in row['unrecognized_declared_atoms']:
                    atom = mol.GetAtomWithIdx(index)
                    neighbor_ranks = [ranks[a.GetIdx()] for a in atom.GetNeighbors()]
                    row['unrecognized_atom_topologies'].append({
                        'index': index, 'degree': len(neighbor_ranks),
                        'distinct_neighbor_ranks': len(set(neighbor_ranks)),
                        'tag': str(atom.GetChiralTag())})
        except Exception as exc:
            row['error'] = repr(exc)
        rows.append(row)
    result = {
        'source_sha256': source_hash,
        'split_sha256': hashlib.sha256(split_bytes).hexdigest(),
        'rdkit_version': rdBase.rdkitVersion,
        'platform': platform.platform(),
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'partition': args.partition,
        'records': len(rows),
        'errors': sum('error' in row for row in rows),
        'specified_types': dict(type_counts),
        'unrecognized_declared_atoms': sum(len(row.get('unrecognized_declared_atoms', [])) for row in rows),
        'unrecognized_declared_bonds': sum(len(row.get('unrecognized_declared_bonds', [])) for row in rows),
        'unrecognized_four_distinct_tetrahedral': sum(
            1 for row in rows for atom in row.get('unrecognized_atom_topologies', [])
            if atom['degree'] == atom['distinct_neighbor_ranks'] == 4
            and atom['tag'] in ('CHI_TETRAHEDRAL_CW', 'CHI_TETRAHEDRAL_CCW')),
        'tetrahedral_without_cip': sum(
            1 for row in rows for item in row['elements']
            if item['type'] == 'Atom_Tetrahedral' and not item['cip_available']),
        'records_with_tetrahedral_without_cip': sum(
            any(item['type'] == 'Atom_Tetrahedral' and not item['cip_available']
                for item in row['elements']) for row in rows),
        'results': rows,
    }
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'results'}))
    if result['errors']:
        raise SystemExit('stereo audit encountered input errors')


if __name__ == '__main__':
    main()
