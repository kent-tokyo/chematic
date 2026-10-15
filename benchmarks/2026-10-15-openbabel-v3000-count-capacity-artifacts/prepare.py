from pathlib import Path
import json,shutil
r=Path('/private/tmp/chematic-openbabel-count-capacity-20261015');old=Path('/private/tmp/chematic-openbabel-zero-z6-20261015')
for n in ['run_formal.py','source_diagnostics.py','measure.py','generate_continuation_corpus.py','deliver.py']:
 t=(old/n).read_text().replace('chematic-openbabel-zero-z6-source-20261015','chematic-openbabel-count-capacity-source-20261015').replace('chematic-openbabel-zero-z6-target-20261015','chematic-openbabel-count-capacity-target-20261015').replace('chematic-openbabel-zero-z6-20261015','chematic-openbabel-count-capacity-20261015').replace('openbabel-v3000-zero-z6','openbabel-v3000-count-capacity').replace('direct-lines-copy','count-capacity-copy');(r/n).write_text(t)
for n in ['continuation-plan.json','continuation-corpus-manifest.json','corpus-continued-chains.molblocks']:shutil.copyfile(old/n,r/n)
p=json.loads((old/'measurement-plan.json').read_text());p.pop('revision_reason',None);p['intended_improvement']='V3000 parsing/roundtrip via one bounded capacity reservation at the ATOM transition';p['baseline_scope']='All nine accepted source candidates, freshly built before this candidate';(r/'measurement-plan.json').write_text(json.dumps(p,indent=2)+'\n')
print('Fixed matrix and adaptive policy recorded before any timing')
