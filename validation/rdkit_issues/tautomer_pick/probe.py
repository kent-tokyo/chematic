"""Run with the candidate native package and one explicit RDKit on PYTHONPATH."""
import argparse
import gc
import hashlib
import itertools
import json
from pathlib import Path
import resource
import subprocess
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--output', type=Path)
parser.add_argument('--rdkit-empty', action='store_true')
args = parser.parse_args()

if args.rdkit_empty:
    # An upstream empty-input check is isolated and cannot emit a core dump.
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    from rdkit.Chem.MolStandardize import rdMolStandardize
    try:
        result = rdMolStandardize.TautomerEnumerator().PickCanonical([])
        print(json.dumps({'returned_none': result is None}))
    except Exception as exc:
        print(json.dumps({'exception': type(exc).__name__, 'message': str(exc)}))
    sys.exit(0)

import chematic as c
from rdkit import Chem, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize as S

families = [
    ['C=C(C)O', 'CC(C)=O'],
    ['Oc1ccccn1', 'O=c1cccc[nH]1', 'O=C1CC=CC=N1'],
    ['CC=O', 'C=CO'],
    ['CC(=O)CC(C)=O', 'CC(O)=CC(C)=O', 'C=C(O)CC(C)=O'],
    ['CC(=O)N', 'CC(O)=N'],
    ['CC(=S)N', 'CC(S)=N'],
    ['NC(=N)N', 'N=C(N)N'],
    ['C[C@H](F)Cl', 'C[C@@H](F)Cl'],
    ['C/C=C/C', 'C/C=C\\C'],
    ['[13CH3]CO', 'CC[18OH]'],
    ['[CH3:9][CH2:7][OH:3]', '[OH:3][CH2:7][CH3:9]'],
    ['[NH4+]', 'C[NH3+]'],
    ['', 'C'],
]
pool = ['CCO', 'CCN', 'CCS', 'CCC', 'CC(C)O', 'CC(C)=O', 'C=C(C)O', 'c1ccccc1']
cases = []
for family in families:
    cases.extend(list(p) for p in itertools.permutations(family))
    cases.append(family + [family[0]])
for a in pool:
    for b in pool:
        cases.extend([[a, b], [a, b, a]])
# Includes a non-tautomer mixed set intentionally: the API selects candidates,
# it makes no assertion that they are members of the same tautomer family.
rows = []
te = S.TautomerEnumerator()
for strings in cases:
    rms = [Chem.MolFromSmiles(s) for s in strings]
    assert all(m is not None for m in rms), strings
    scores = [te.ScoreTautomer(m) for m in rms]
    keys = [Chem.MolToSmiles(m) for m in rms]
    expected_index = min(range(len(strings)), key=lambda i: (-scores[i], keys[i], i))
    expected = keys[expected_index]
    oracle = Chem.MolToSmiles(te.PickCanonical(rms))
    assert oracle == expected, (strings, oracle, expected)
    native_scores = [(c.from_smiles(s) if s else c.RWMol().GetMol()).rdkit_tautomer_score() for s in strings]
    assert native_scores == scores, (strings, native_scores, scores)
    for kind in ['list', 'tuple', 'generator']:
        def generated():
            for i, s in enumerate(strings):
                mol = (c.from_smiles(s) if s else c.RWMol().GetMol())
                mol.SetIntProp('source_index', i)
                yield mol
        items = generated()
        if kind == 'list':
            items = list(items)
        elif kind == 'tuple':
            items = tuple(items)
        chosen = c.rdkit_pick_tautomer(items, max_candidates=len(strings))
        gc.collect()
        assert chosen.rdkit_smiles == expected, (strings, kind)
        assert chosen.GetIntProp('source_index') == expected_index, (strings, kind)
    rows.append({'inputs': strings, 'scores': scores, 'keys': keys,
                 'chosen_index': expected_index, 'chosen_smiles': expected})

policy_strings = ['C[C@H](F)Cl', 'CC(F)Cl']
policy_native = c.rdkit_pick_tautomer(map(c.from_smiles, policy_strings)).rdkit_smiles
policy_rdkit = Chem.MolToSmiles(te.PickCanonical([Chem.MolFromSmiles(s) for s in policy_strings]))
assert policy_native == 'CC(F)Cl'

upstream_empty = subprocess.run([sys.executable, str(Path(__file__).resolve()), '--rdkit-empty'],
                                capture_output=True, text=True, timeout=15)
try:
    c.rdkit_pick_tautomer(iter(()))
except ValueError as exc:
    empty_error = str(exc)
else:
    raise AssertionError('empty candidate set did not raise ValueError')
import chematic.chematic as native
native_path = Path(native.__file__)
result = {'rdkit_version': rdBase.rdkitVersion, 'rdkit_path': Chem.__file__,
          'native_path': str(native_path), 'native_sha256': hashlib.sha256(native_path.read_bytes()).hexdigest(),
          'candidate_sets': len(rows), 'selection_comparisons': len(rows) * 3,
          'native_score_comparisons': sum(len(r['inputs']) for r in rows),
          'upstream_empty_input': {'returncode': upstream_empty.returncode,
                                   'stdout': upstream_empty.stdout, 'stderr': upstream_empty.stderr},
          'chematic_empty_error': empty_error,
          'outside_parity_scope_mixed_stereo_tie': {'inputs': policy_strings, 'chematic': policy_native, 'rdkit': policy_rdkit},
          'failures': 0, 'rows': rows}
assert args.output is not None
args.output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k:v for k,v in result.items() if k != 'rows'}, indent=2))
