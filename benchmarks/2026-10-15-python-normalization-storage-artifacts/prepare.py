from pathlib import Path
import json
r=Path('/private/tmp/chematic-descriptor-profile-20261015');old=Path('/private/tmp/chematic-aromatic-components-20261015')
for name in ['run-frontiers.py','run-followups.py','measure.py','deliver.py','freeze.py']:
 txt=(old/name).read_text().replace('chematic-aromatic-components-source-20261015','chematic-descriptor-profile-source-20261015').replace('chematic-aromatic-components-target-20261015','chematic-descriptor-profile-target-20261015').replace('chematic-aromatic-components-20261015','chematic-descriptor-profile-20261015')
 txt=txt.replace('plus-six-accepted-source-candidates','plus-eight-accepted-source-candidates').replace('python-aromatic-component-reuse','python-normalization-storage').replace('aromatic-component-reuse','normalization-storage').replace('aromatic-component-copy','normalization-storage-copy')
 if name=='freeze.py':
  txt=txt.replace('crates/chematic-perception/src/rdkit_parity/typing.rs','crates/chematic-perception/src/aromaticity/model.rs').replace('crates/chematic-perception/src/rdkit_parity/component_typing_tests.rs','crates/chematic-perception/src/aromaticity/normalization_storage_tests.rs').replace('baseline-typing.rs','baseline-model.rs').replace('aromatic-component-source-candidate/v1','normalization-storage-source-candidate/v1')
  txt=txt.replace('Fresh verified six-candidate baseline was built before editing. First native-test link failed with ENOSPC; its log is retained. Only regenerable targets owned by this run build targets were removed between validation/build phases; successful retry logs retained.','Fresh verified eight-candidate baseline was built before editing and CPU-sampled. Only owned completed Cargo targets are reclaimed between native validation, wheel build and timing; all source, wheels and logs retained.')
 (r/name).write_text(txt)
p=json.loads((old/'measurement-plan.json').read_text());p['intended_improvement']='First-use TPSA and Labute ASA aromatic normalization; capacities and indexed hydrogen patch map only';p['immediate_baseline']='All eight previously accepted source candidates, freshly built';(r/'measurement-plan.json').write_text(json.dumps(p,indent=2)+'\n')
print('Measurement plan and pipeline prepared before any candidate timing')
