#!/usr/bin/env python3
"""Observe all-frame raw coordinates from common HISTORY files.

chematic inspects a structured record; Open Babel exports one CJSON per frame
with -m. Only canonical C/O labels and nonperiodic inputs are compared. This
does not attribute CJSON serialization effects to the HISTORY parser, or score
unit interpretation, masses, charges, velocities, forces or source metadata.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def geometry(elements, flat):
    if not elements or len(flat) != 3*len(elements) or any(not math.isfinite(v) for v in flat):
        raise ValueError('invalid coordinate frame')
    if any(v not in (6,8) for v in elements):
        raise ValueError('comparison is restricted to canonical C/O labels')
    return dict(elements=elements, positions=[flat[i:i+3] for i in range(0,len(flat),3)])


def observe_record(g):
    result = []
    for f in g['frames']:
        result.append(geometry([{'C':6,'O':8}[a['name']] for a in f['atoms']],
                               [v for a in f['atoms'] for v in a['position']]))
    if not result:
        raise ValueError('no frames')
    return result


def observe_cjson(documents):
    result = []
    for g in documents:
        atoms = g['atoms']
        result.append(geometry(atoms['elements']['number'], atoms['coords']['3d']))
    if not result:
        raise ValueError('no CJSON frames')
    return result


def execute(command):
    try:
        p = subprocess.run(command, text=True, encoding='utf-8', capture_output=True, timeout=60, check=False)
    except (OSError, subprocess.TimeoutExpired) as e:
        return dict(status='execution_failed', error=str(e), command=command)
    return dict(status='returned' if p.returncode==0 else 'rejected', command=command,
                returncode=p.returncode, stdout=p.stdout, stderr=p.stderr)


def observe_case(cli, obabel, case):
    expected = observe_record(case['expected'])
    row = dict(id=case['id'], input_sha256=hashlib.sha256(case['input'].encode()).hexdigest(), expected=expected)
    with tempfile.TemporaryDirectory(prefix='chematic-history-comparison-') as directory:
        directory = Path(directory)
        source = directory/'input.HISTORY'
        source.write_text(case['input'], encoding='ascii')
        if hashlib.sha256(source.read_bytes()).hexdigest() != row['input_sha256']:
            raise ValueError('temporary input differs')
        row['chematic'] = execute([str(cli),'dlpoly-history','--inspect','--input',str(source)])
        row['openbabel'] = execute([str(obabel),'-iHISTORY',str(source),'-ocjson','-O',str(directory/'frame.cjson'),'-m'])
        for engine in ('chematic','openbabel'):
            result = row[engine]
            if result['status'] != 'returned':
                continue
            try:
                if engine == 'chematic':
                    observed = observe_record(json.loads(result['stdout']))
                else:
                    paths = list(directory.glob('frame*.cjson'))
                    paths.sort(key=lambda p:int(re.fullmatch(r'frame(\d+)\.cjson',p.name)[1]))
                    receipts = [dict(name=p.name, text=p.read_text(encoding='utf-8')) for p in paths]
                    result['output_files'] = receipts
                    observed = observe_cjson([json.loads(p['text']) for p in receipts])
                result.update(status='observed', geometry=observed,
                              frame_count_preserved=len(observed)==len(expected), exact_all_frame_coordinates=observed==expected)
            except (ValueError, KeyError, IndexError, TypeError) as e:
                result.update(status='unobservable', error=str(e))
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cli', type=Path, required=True)
    parser.add_argument('--obabel', type=Path, default=Path('/opt/homebrew/bin/obabel'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    path = ROOT/'validation/dlpoly_history/cases.json'
    packet = json.loads(path.read_text(encoding='utf-8'))
    cases = [c for c in packet['valid'] if c['compare_coordinates']]
    for c in cases:
        if c['expected']['imcon'] != 0:
            raise ValueError('periodic comparison is outside this observer scope')
    cli, obabel = args.cli.resolve(), args.obabel.resolve()
    rows = [observe_case(cli,obabel,c) for c in cases]
    sources = ['crates/chematic-mol/src/dlpoly_history.rs','crates/chematic-mol/src/dlpoly_history/reader.rs','crates/chematic-mol/src/dlpoly_history/writer.rs','crates/chematic-cli/src/main.rs','crates/chematic-cli/src/dispatch.rs']
    report = dict(schema_version=1, scope='Six original nonperiodic C/O Classic HISTORY inputs; all-frame raw coordinate observation through native inspection versus Open Babel CJSON export with -m; not a pure parser or full semantic-preservation comparison',
                  case_count=len(rows), timeout_seconds=60, input_mode='same seekable temporary file for both tools',
                  command_templates=dict(chematic=['<cli>','dlpoly-history','--inspect','--input','<input.HISTORY>'],openbabel=['<obabel>','-iHISTORY','<input.HISTORY>','-ocjson','-O','<frame.cjson>','-m']),
                  observed_counts={k:sum(c[k]['status']=='observed' for c in rows) for k in ('chematic','openbabel')},
                  exact_all_frame_counts={k:sum(c[k].get('exact_all_frame_coordinates') is True for c in rows) for k in ('chematic','openbabel')},
                  versions=dict(chematic=subprocess.check_output([str(cli),'--version'],text=True).strip(),openbabel=subprocess.check_output([str(obabel),'-V'],text=True).strip()),
                  packet_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
                  binary_sha256=dict(chematic=hashlib.sha256(cli.read_bytes()).hexdigest(),openbabel=hashlib.sha256(obabel.read_bytes()).hexdigest()),cases=rows)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('case_count','observed_counts','exact_all_frame_counts')},indent=2))
    return int(report['exact_all_frame_counts']['chematic'] != len(rows))


if __name__=='__main__':
    raise SystemExit(main())
