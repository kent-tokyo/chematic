#!/usr/bin/env python3
"""Paired ORCA writer observations with authored expectations, no engine oracle."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.checks.check_interchange_preservation import bounded_run, identity, require, strict_json
from scripts.validation.orca_checked_observer import outcome

WORKER = '''
import json,sys
sys.path.insert(0,sys.argv[1])
import chematic
records=[]
for c in json.load(sys.stdin):
 if c['lane']=='transport':continue
 doc=json.loads(c['input_json']);record={'id':c['id']}
 calls={'checked':chematic.write_orca_input_checked}
 if c['expected']=='accept':calls['legacy']=chematic.write_orca_input
 for name,call in calls.items():
  try:record[name]={'ok':True,'output':call(doc)}
  except (ValueError,TypeError) as e:record[name]={'ok':False,'error':str(e)}
 records.append(record)
print(json.dumps(records,allow_nan=False))
'''


def run(command, data=b''):
    result = bounded_run([str(x) for x in command], data, timeout=30, max_output=1_048_576)
    require(result['resource_status'] == 'completed', 'process resource budget exceeded')
    return {**result, **{k: result[k].decode('utf8', errors='replace') for k in ('stdout', 'stderr')}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python-root', type=Path, required=True)
    parser.add_argument('--wasm-module', type=Path, required=True)
    parser.add_argument('--obabel', type=Path, default=Path('/opt/homebrew/bin/obabel'))
    parser.add_argument('--output', type=Path, default=ROOT/'validation/results/orca-checked-2026-10-24.json')
    args = parser.parse_args()
    fixture = ROOT/'validation/orca-checked-2026-10-24/cases.json'
    cases = strict_json(fixture.read_bytes())['cases']
    raw, records, counts, maximum = {}, [], Counter(), Counter()
    for surface, command in [
        ('python', [sys.executable, '-c', WORKER, args.python_root]),
        ('wasm', ['node', ROOT/'scripts/validation/validate_orca_checked_node.mjs', args.wasm_module]),
    ]:
        result = run(command, json.dumps(cases).encode())
        require(result['exit_code'] == 0, f'{surface} worker failed: {result}')
        raw[surface] = result
        observed = strict_json(result['stdout'].encode())
        selected = [c for c in cases if surface == 'wasm' or c['lane'] != 'transport']
        require(len(selected) == len(observed), 'worker dropped rows')
        for case, row in zip(selected, observed):
            require(case['id'] == row['id'], 'worker reordered rows')
            if case['expected'] == 'refuse':
                require(outcome(case, row['checked'])['status'] == 'refused', 'unsafe checked output')
                counts[f'{surface}.invalid.refused'] += 1
            else:
                for engine in ('checked', 'legacy'):
                    observation = outcome(case, row[engine])
                    if engine == 'checked':
                        require(observation['status'] == 'fields_match', f'{surface} {case["id"]}: {observation}')
                    counts[f'{surface}.{engine}.{case["lane"]}.{observation["status"]}'] += 1
                    if observation.get('metadata_match'): counts[f'{surface}.{engine}.metadata_match'] += 1
                    if case['lane'] == 'common':
                        maximum[f'{surface}.{engine}'] = max(maximum[f'{surface}.{engine}'], observation.get('max_coordinate_error', 0))
                    records.append(dict(surface=surface, engine=engine, id=case['id'], observation=observation))
            records.append(dict(surface=surface, id=case['id'], raw=row))
    for case in [c for c in cases if c['lane'] == 'common']:
        atoms = case['fields']['coords']['atoms']
        xyz = str(len(atoms))+'\noriginal geometry\n'+''.join(f'{a["element"]} {a["x"]!r} {a["y"]!r} {a["z"]!r}\n' for a in atoms)
        result = run([args.obabel, '-ixyz', '-oorcainp', '-ab'], xyz.encode())
        envelope = dict(ok=True, output=result['stdout']) if result['exit_code'] == 0 else dict(ok=False, error=result['stderr'])
        observed = outcome(case, envelope, common_only=True)
        counts[f'openbabel.common.{observed["status"]}'] += 1
        maximum['openbabel'] = max(maximum['openbabel'], observed.get('max_coordinate_error', 0))
        records.append(dict(surface='openbabel', id=case['id'], input_xyz=xyz, process=result, observation=observed))
    base = '38b4c92c430d9b3bfe70e47e97293b9017ced990'
    attestation = {}
    for path, start, end in [
        ('crates/chematic-mol/src/orca.rs', 'pub fn write_orca_input(', '// Output: errors'),
        ('crates/chematic-wasm/src/format_io.rs', 'pub fn write_orca_input_json(', '// ORCA output (.out'),
        ('crates/chematic-py/src/formats/materials.rs', 'pub(super) fn write_orca_input(', '// ORCA output\n'),
    ]:
        previous = subprocess.check_output(['git', 'show', f'{base}:{path}'], cwd=ROOT, text=True)
        current = (ROOT/path).read_text()
        require(start in previous and start in current and end in previous and end in current, 'missing attestation anchor')
        old = previous.split(start, 1)[1].split(end, 1)[0]
        new = current.split(start, 1)[1].split(end, 1)[0]
        require(old == new, 'legacy writer changed from pinned base')
        attestation[path] = hashlib.sha256((start+new).encode()).hexdigest()
    packet = args.output.with_suffix('.packet.json.gz')
    packet.write_bytes(gzip.compress(json.dumps(dict(worker_outputs=raw, records=records), indent=2, allow_nan=False).encode(), mtime=0))
    paths = ['Cargo.lock', 'crates/chematic-mol/src/orca/checked.rs', 'crates/chematic-mol/src/orca/checked_json.rs', 'crates/chematic-mol/src/orca/cartesian.rs', 'crates/chematic-mol/src/orca.rs', 'crates/chematic-mol/src/lib.rs', 'crates/chematic-py/src/formats/orca_checked.rs', 'crates/chematic-py/src/formats/registration.rs', 'crates/chematic-py/src/formats.rs', 'crates/chematic-wasm/src/format_io/orca_checked.rs', 'crates/chematic-wasm/src/format_io.rs', 'scripts/validation/validate_orca_checked.py', 'scripts/validation/validate_orca_checked_node.mjs', 'scripts/validation/orca_checked_observer.py']
    extension = next((args.python_root/'chematic').glob('chematic*.so'))
    report = dict(schema='chematic.orca-checked-validation.v1', utc=datetime.now(timezone.utc).isoformat(),
                  platform=platform.platform(), candidate_base=base, candidate_state='working tree, not a published release',
                  scope='Represented structured xyz input fields and f64 values. Opaque settings are not executable-job validation or source identity.',
                  counts=dict(sorted(counts.items())), common_max_coordinate_error=dict(maximum),
                  fixture=identity(fixture), packet=identity(packet), sources=[identity(ROOT/p) for p in paths],
                  legacy_control=dict(scope='Old writer APIs in the same candidate runtime. Source bodies equal pinned main; not a separately rebuilt main/release binary.', source_commit=base, source_section_sha256=attestation),
                  artifacts=dict(python_extension=identity(extension), python_wrapper=identity(args.python_root/'chematic/__init__.py'), wasm_js=identity(args.wasm_module), wasm_binary=identity(args.wasm_module.with_name('chematic_wasm_bg.wasm')), openbabel=identity(args.obabel.resolve())),
                  versions=dict(python=sys.version, rust=run(['rustc','--version']), node=run(['node','--version']), openbabel=run([args.obabel,'-V'])),
                  comparison_notes=['Common lane: 12 identical ordered element/coordinate arrays, Angstrom, charge 0/spin 1. Open Babel is XYZ-to-ORCA pipeline, not identical internal state.', '12 profile/numeric positives and 38 raw-JSON refusals are candidate controls, not additional Open Babel losses. Python cannot represent duplicate JSON keys and exercises 30 dictionary refusals.', 'No ORCA execution, chemistry validation, timing, release or toolkit-wide superiority claim.'],
                  process_limits=dict(seconds=30, combined_pipe_bytes=1_048_576, memory='not imposed'))
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(counts=report['counts'], common_max_coordinate_error=dict(maximum)), indent=2))


if __name__ == '__main__':
    main()
