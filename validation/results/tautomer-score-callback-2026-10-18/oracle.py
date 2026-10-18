"""Pure scoring callbacks, supplied acyclic families and isolated policy cases."""
from pathlib import Path
import json
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize as S
import chematic as C

root = Path(__file__).resolve().parent
families = [
    ['CC(C)=O','C=C(C)O'], ['CC=O','C=CO'],
    ['CCC(C)=O','CC=C(C)O','CCC(O)=C'],
    ['CC(=O)CC(C)=O','CC(O)=CC(C)=O','CC(=O)CC(O)=C'],
    ['NC(N)=O','N=C(N)O'], ['NC(N)=S','N=C(N)S'],
    ['CC(=O)N','CC(O)=N'], ['CC(=O)O','C=C(O)O'],
]
rng = np.random.default_rng(960119)
summary = dict(rdkit_version=rdkit.__version__, seed=960119, trials=128,
               comparisons=0, callback_calls_candidate=0, callback_calls_reference=0,
               index_mismatches=0, graph_mismatches=0)
enumerator = S.TautomerEnumerator()
for trial in range(summary['trials']):
    family = list(families[trial % len(families)])
    rng.shuffle(family)
    values = [int(x) for x in rng.integers(-100,101,size=len(family))]
    for mode in ('default','negative_default','property','constant'):
        for source_kind in ('list','tuple','generator'):
            candidate = [C.from_smiles(s) for s in family]
            reference = [Chem.MolFromSmiles(s) for s in family]
            for i,(a,b) in enumerate(zip(candidate,reference)):
                a.SetIntProp('index',i); b.SetIntProp('index',i)
                a.SetIntProp('custom',values[i]); b.SetIntProp('custom',values[i])
            actual_calls, reference_calls = [],[]
            def actual_score(m):
                actual_calls.append(m.GetIntProp('index'))
                if mode == 'property': return m.GetIntProp('custom')
                if mode == 'constant': return 0
                return m.rdkit_tautomer_score() * (-1 if mode == 'negative_default' else 1)
            def reference_score(m):
                reference_calls.append(m.GetIntProp('index'))
                if mode == 'property': return m.GetIntProp('custom')
                if mode == 'constant': return 0
                return S.TautomerEnumerator.ScoreTautomer(m) * (-1 if mode == 'negative_default' else 1)
            source = candidate if source_kind == 'list' else tuple(candidate) if source_kind == 'tuple' else (m for m in candidate)
            actual = C.rdkit_pick_tautomer_scored(source, actual_score)
            expected = enumerator.PickCanonical(reference, reference_score)
            assert actual_calls == reference_calls == list(range(len(family)))
            assert actual.GetIntProp('index') == expected.GetIntProp('index'), (family,mode,actual.rdkit_smiles,Chem.MolToSmiles(expected))
            assert Chem.MolToSmiles(Chem.MolFromSmiles(actual.rdkit_smiles)) == Chem.MolToSmiles(expected)
            summary['comparisons'] += 1
            summary['callback_calls_candidate'] += len(actual_calls)
            summary['callback_calls_reference'] += len(reference_calls)

# Known deliberate difference: this API keeps the existing lexical policy;
# newer RDKit prefers specified stereo before lexical ties.
inputs = ['CC(F)O','C[C@H](F)O']
actual = C.rdkit_pick_tautomer_scored([C.from_smiles(s) for s in inputs],lambda _:0)
expected = enumerator.PickCanonical([Chem.MolFromSmiles(s) for s in inputs],lambda _:0)
difference = dict(name='specified_vs_unspecified_stereo_tie', inputs=inputs,
                  candidate=actual.rdkit_smiles, reference=Chem.MolToSmiles(expected),
                  counted_as_parity=False)
assert '@' not in actual.rdkit_smiles and '@' in Chem.MolToSmiles(expected)
(root/'separate-diagnostics.json').write_text(json.dumps([difference],indent=2)+'\n')
(root/'oracle-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
print('Separate lexical/stereo policy difference:',difference)
