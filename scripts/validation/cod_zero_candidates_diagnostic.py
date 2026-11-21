#!/usr/bin/env python3
"""Inspect stereo-verifier diagnostics for exposed zero-candidate COD inputs."""
import gzip
import json
from pathlib import Path
import chematic
from rdkit import rdBase

with gzip.open('validation/results/cod-linux-development-v4-2026-11-22.json.gz', 'rt') as source:
    rows = json.load(source)['results']
config = chematic.PipelineV2Config.safe(
    force_field='none', stereo_policy='verify_only', ring_torsion_policy='diagnostic_only',
    fail_on_unevaluable_stereo=True, max_attempts=1, total_timeout_ms=20000, embed_seed=42)
results = []
for row in rows:
    if row['status'] != 'ok' or row.get('chematic_count', 0) != 0:
        continue
    item = {'name': row['name'], 'rdkit_count': row.get('rdkit_count'),
            'rdkit_error': row.get('rdkit_error')}
    try:
        chematic.from_smiles(row['smiles']).add_hydrogens().embed_pipeline_v2(config)
        item['pipeline_outcome'] = 'no_strict_failure'
    except Exception as exc:
        diagnostic = getattr(exc, 'diagnostics', {})
        item['pipeline_outcome'] = diagnostic.get('cause', {'error': repr(exc)})
        item['final_stereo'] = diagnostic.get('final_stereo')
    results.append(item)
assert len(results) == 8
out = {'rdkit_version': rdBase.rdkitVersion, 'method': 'pipeline_v2, one seed; not the corrected ETKDG ensemble',
       'results': results}
Path('cod-zero-candidates-diagnostic.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out))
