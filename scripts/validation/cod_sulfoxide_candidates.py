#!/usr/bin/env python3
"""Development-only follow-up on already exposed v4 zero-candidate inputs."""
import json
from pathlib import Path
import sys

import chematic
from rdkit import Chem, rdBase

sys.path.insert(0, str(Path.cwd() / 'scripts/benchmarks'))
from platinum_development_probe import stereo_violations

fixture = json.loads(Path('validation/results/cod-sulfoxide-zero-development-2026-11-22.json').read_text())
rows = []
for case in fixture['results']:
    row = {'name': case['name'], 'v4_chematic_count': case['v4_chematic_count'],
           'v4_rdkit_count': case['v4_rdkit_count']}
    try:
        ref = Chem.AddHs(Chem.MolFromSmiles(case['smiles']))
        mol = chematic.from_smiles(case['smiles']).add_hydrogens()
        conformers = mol.etkdg_embed_multiple_confs_mmff_corrected(
            20, random_seed=42, max_mmff_iters=200, repair_mmff_iters=400)
        checks = stereo_violations(ref, [coords for _, _, coords in conformers])
        row.update(candidate_count=len(conformers),
                   unfinished=sum(status != 0 for status, _, _ in conformers),
                   stereo_violations=sum(bool(c['atoms'] or c['bonds']) for c in checks),
                   stereo_details=checks)
    except Exception as e:
        row['error'] = repr(e)
    rows.append(row)
summary = {'records': len(rows), 'errors': sum('error' in x for x in rows),
           'records_with_candidates': sum(x.get('candidate_count', 0) > 0 for x in rows),
           'candidates': sum(x.get('candidate_count', 0) for x in rows),
           'unfinished': sum(x.get('unfinished', 0) for x in rows),
           'stereo_violations': sum(x.get('stereo_violations', 0) for x in rows)}
output = {'classification': fixture['classification'], 'rdkit': rdBase.rdkitVersion,
          'summary': summary, 'results': rows}
Path('cod-sulfoxide-candidates.json').write_text(json.dumps(output, indent=2) + '\n')
print(json.dumps(summary))
assert not summary['errors'] and not summary['stereo_violations']
