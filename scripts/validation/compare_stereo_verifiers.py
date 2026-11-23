#!/usr/bin/env python3
"""Compare old and strict declared-stereo checks on exposed development inputs."""
import argparse
import gzip
import hashlib
import json
import platform
import sys
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))
from platinum_development_probe import stereo_violations  # noqa: E402
from platinum_stereo_strict import _raw_tetrahedral_parity, stereo_violations_strict  # noqa: E402


def geometry(mol):
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.numThreads = 1
    params.maxIterations = 20
    if AllChem.EmbedMolecule(mol, params) != 0:
        return None
    conformer = mol.GetConformer()
    coords = [list(conformer.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    return [coords, [[-x, y, z] for x, y, z in coords]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('development', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    RDLogger.DisableLog('rdApp.*')
    source_bytes = gzip.decompress(args.development.read_bytes())
    source = json.loads(source_bytes)
    if hashlib.sha256(source_bytes).hexdigest() != 'e81fd9c1eba4c8fa218308ad1e312939b8c3e7ec306684711f2341449d9448bb':
        raise ValueError('development input fingerprint mismatch')
    if source.get('profile') != 'corrected' or source.get('candidate_limit') != 20:
        raise ValueError('expected the fixed 20-candidate development input')
    old_counts, parity_counts = Counter(), Counter()
    old_names, parity_names = [], []
    for row in source['results']:
        if not row.get('smiles'):
            continue
        mol = Chem.AddHs(Chem.MolFromSmiles(row['smiles']))
        if mol.GetNumAtoms() > 70:
            continue
        specified = [s for s in Chem.FindPotentialStereo(mol)
                     if s.specified == Chem.StereoSpecified.Specified]
        centers = [a.GetIdx() for a in mol.GetAtoms()
                   if a.GetChiralTag() in (Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                                           Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
                   and a.GetDegree() == 4]
        recognized = {s.centeredOn for s in specified if str(s.type) == 'Atom_Tetrahedral'}
        need_old = len(old_names) < 30 and bool(specified)
        need_parity = len(parity_names) < 50 and bool(centers)
        if not (need_old or need_parity):
            continue
        conformers = geometry(mol)
        if conformers is None:
            continue
        ranks = (list(Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=True))
                 if need_parity else None)
        for coords in conformers:
            strict = stereo_violations_strict(mol, [coords])[0]
            if need_old:
                old = stereo_violations(mol, [coords])[0]
                old_atoms, old_bonds = set(old['atoms']), set(old['bonds'])
                old_counts['old_only'] += bool(old_atoms - set(strict['atoms'])
                                               or old_bonds - set(strict['bonds']))
                old_counts['additional_contradictions'] += bool(
                    set(strict['violated_atoms']) - old_atoms
                    or set(strict['violated_bonds']) - old_bonds)
                old_counts['additional_unknowns'] += bool(
                    set(strict['unevaluable_atoms']) - old_atoms
                    or set(strict['unevaluable_bonds']) - old_bonds)
            if need_parity:
                for center in centers:
                    if center not in recognized:
                        parity_counts['rdkit_not_recognized'] += 1
                        continue
                    raw = _raw_tetrahedral_parity(mol, center, coords, ranks)
                    if raw is None:
                        parity_counts['raw_unknown'] += 1
                    elif center in strict['unevaluable_atoms']:
                        parity_counts['rdkit_unknown'] += 1
                    else:
                        parity_counts['compared'] += 1
                        parity_counts['disagreement'] += raw != (center not in strict['violated_atoms'])
        if need_old:
            old_names.append(row['name'])
        if need_parity:
            parity_names.append(row['name'])
        if len(old_names) == 30 and len(parity_names) == 50:
            break
    result = {'development_sha256': hashlib.sha256(source_bytes).hexdigest(),
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'rdkit_version': rdBase.rdkitVersion, 'platform': platform.platform(),
              'selection': 'first eligible exposed development records; at most 70 atoms',
              'seed': 42, 'geometries': 'RDKit ETKDGv3 plus x-mirror',
              'old_comparison': {'molecules': len(old_names), 'names': old_names,
                                 'counts': dict(old_counts)},
              'raw_parity_comparison': {'molecules': len(parity_names), 'names': parity_names,
                                        'counts': dict(parity_counts)}}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('old_comparison', 'raw_parity_comparison')}))
    print(json.dumps({'old': result['old_comparison']['counts'],
                      'raw_parity': result['raw_parity_comparison']['counts']}))
    if (len(old_names) != 30 or len(parity_names) != 50
            or old_counts['old_only'] or parity_counts['disagreement']):
        raise SystemExit('stereo verifier comparison did not pass')


if __name__ == '__main__':
    main()
