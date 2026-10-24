#!/usr/bin/env python3
"""Compare XYZ-to-NWChem Cartesian geometry export, not an nw round trip.

Open Babel advertises nw as write-only. The chematic JSON adapter is verified
against the same indexed XYZ tags/coordinates before either tool is run. The
observer examines geometry tags, coordinates and units only; it does not
validate the surrounding NWChem job deck or infer missing task information.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def xyz_geometry(text):
    rows = text.splitlines()
    count = int(rows[0])
    if count <= 0 or len(rows) != count+2:
        raise ValueError("invalid shared XYZ frame")
    tags, positions = [], []
    for row in rows[2:]:
        w = row.split()
        if len(w) != 4:
            raise ValueError("invalid shared XYZ atom")
        tags.append(w[0])
        positions.append([float(s) for s in w[1:]])
    if not all(math.isfinite(n) for v in positions for n in v):
        raise ValueError("nonfinite shared XYZ geometry")
    return dict(tags=tags, positions_angstrom=positions)


def observe(text):
    blocks, rows, header = [], None, None
    for raw in text.splitlines():
        w = raw.split('#', 1)[0].split()
        if not w:
            continue
        if w[0].lower() == "geometry":
            if rows is not None:
                raise ValueError("nested geometry")
            rows, header = [], w
        elif rows is not None and w[0].lower() == "end":
            if len(w) != 1:
                raise ValueError("invalid END")
            blocks.append((header, rows))
            rows = None
        elif rows is not None:
            if len(w) != 4:
                raise ValueError("expected simple Cartesian export")
            rows.append(w)
    if rows is not None or len(blocks) != 1 or not blocks[0][1]:
        raise ValueError("missing, multiple, or incomplete geometry blocks")
    header, rows = blocks[0]
    lower = [s.lower() for s in header]
    unit = "angstroms"
    if "units" in lower:
        unit = lower[lower.index("units")+1]
    if unit not in ("an", "angstroms"):
        raise ValueError("shared XYZ export must declare/default to angstroms")
    positions = [[float(s.replace('D', 'e').replace('d', 'e')) for s in w[1:]] for w in rows]
    if not all(math.isfinite(n) for v in positions for n in v):
        raise ValueError("nonfinite export")
    return dict(tags=[w[0] for w in rows], positions_angstrom=positions)


def run(command, text, expected):
    try:
        p = subprocess.run(command, input=text, text=True, encoding="utf-8", capture_output=True, timeout=30, check=False)
    except (OSError, subprocess.TimeoutExpired) as e:
        return dict(status="execution_failed", error=str(e))
    receipt = dict(returncode=p.returncode, stdout=p.stdout, stderr=p.stderr,
                   command_input_sha256=hashlib.sha256(text.encode()).hexdigest())
    if p.returncode:
        return dict(status="rejected", **receipt)
    try:
        record = observe(p.stdout)
    except (ValueError, IndexError) as e:
        return dict(status="unobservable", error=str(e), **receipt)
    exact = record == expected
    error = None
    if record["tags"] == expected["tags"]:
        error = max(abs(a-b) for va,vb in zip(record["positions_angstrom"], expected["positions_angstrom"]) for a,b in zip(va,vb))
    return dict(status="observed", geometry=record, exact_indexed_geometry=exact,
                maximum_absolute_component_error_angstrom=error, **receipt)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--obabel", type=Path, default=Path('/opt/homebrew/bin/obabel'))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    path = ROOT/'validation/nwchem_geometry/cases.json'
    packet = json.loads(path.read_text(encoding='utf-8'))
    commands = dict(chematic=[str(args.cli.resolve()), 'nwchem-geometry', '--from-json'],
                    openbabel=[str(args.obabel.resolve()), '-ixyz', '-onw', '-ab', '-l1'])
    rows = []
    for c in packet['xyz_exports']:
        expected = xyz_geometry(c['xyz_input'])
        adapter = dict(tags=[a['tag'] for a in c['expected']['centers']],
                       positions_angstrom=[a['position'] for a in c['expected']['centers']])
        if adapter != expected or c['expected']['units'] != 'angstroms':
            raise ValueError('JSON adapter differs from shared XYZ geometry')
        row = dict(id=c['id'], xyz_sha256=hashlib.sha256(c['xyz_input'].encode()).hexdigest(), expected=expected)
        for engine, command in commands.items():
            text = json.dumps(c['expected'], allow_nan=False) if engine == 'chematic' else c['xyz_input']
            row[engine] = run(command, text, expected)
        rows.append(row)
    sources = ['crates/chematic-mol/src/nwchem.rs', 'crates/chematic-mol/src/nwchem/reader.rs', 'crates/chematic-mol/src/nwchem/writer.rs', 'crates/chematic-cli/src/main.rs', 'crates/chematic-cli/src/dispatch.rs']
    report = dict(schema_version=1, scope='Six original shared indexed XYZ geometries exported to NWChem; direct JSON adapter verified; tags and angstrom coordinates only; nw read/round-trip capability not compared',
                  case_count=len(rows), timeout_seconds=30, commands=commands,
                  observed_counts={k:sum(c[k]['status']=='observed' for c in rows) for k in commands},
                  exact_geometry_counts={k:sum(c[k].get('exact_indexed_geometry') is True for c in rows) for k in commands},
                  versions={k:subprocess.check_output([v[0], '--version' if k=='chematic' else '-V'], text=True).strip() for k,v in commands.items()},
                  packet_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
                  binary_sha256={k:hashlib.sha256(Path(v[0]).read_bytes()).hexdigest() for k,v in commands.items()}, cases=rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('case_count','observed_counts','exact_geometry_counts')}, indent=2))
    return int(report['exact_geometry_counts']['chematic'] != len(rows))


if __name__ == '__main__':
    raise SystemExit(main())
