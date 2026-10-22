#!/usr/bin/env python3
"""Authored CML → CML chemistry-field development packet, not a holdout.

python3 -m scripts.validation.check_cml_character_references --probe PROBE
--obabel OBABEL --output NEW.json. PROBE is file_io_semantic_probe.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import platform
import tempfile
import xml.etree.ElementTree as ET

from scripts.check_interchange_preservation import bounded_run, identity, require

SOURCE = '''<molecule xmlns="http://www.xml-cml.org/schema">
<atomArray>
<atom id="a1" elementType="N" formalCharge="1" isotopeNumber="15" hydrogenCount="4" x2="1.25" y2="-2.5"/>
<atom id="h1" elementType="H" isotopeNumber="2" hydrogenCount="0" x2="2.25" y2="-2.5"/>
</atomArray>
<bondArray>
<bond atomRefs2="a1 h1" order="1"/>
</bondArray>
</molecule>
'''


def packet():
    changes = [
        ("literal-control", None, None, False),
        ("decimal-charge", 'formalCharge="1"', 'formalCharge="&#49;"', False),
        ("hex-isotope", 'isotopeNumber="15"', 'isotopeNumber="&#x31;&#x35;"', False),
        ("decimal-element", 'elementType="N"', 'elementType="&#78;"', False),
        ("decimal-coordinate", 'x2="1.25"', 'x2="&#49;.25"', False),
        ("hex-atom-id", 'id="a1"', 'id="a&#x31;"', False),
        ("hex-bond-order", 'order="1"', 'order="&#x31;"', False),
        ("invalid-charge", 'formalCharge="1"', 'formalCharge="oops"', True),
        ("invalid-isotope", 'isotopeNumber="15"', 'isotopeNumber="oops"', True),
        ("invalid-coordinate", 'x2="1.25"', 'x2="oops"', True),
        ("isotope-alias-conflict", 'isotopeNumber="15"', 'isotope="14" isotopeNumber="15"', True),
        ("duplicate-decoded-id", '<atom id="h1"', '<atom id="a&#49;" elementType="C"/>\n<atom id="h1"', True),
        ("invalid-xml-character", 'formalCharge="1"', 'formalCharge="&#0;"', True),
        ("escaped-reference-is-literal", 'formalCharge="1"', 'formalCharge="&amp;#49;"', True),
    ]
    cases = []
    for case, before, after, refuse in changes:
        text = SOURCE if before is None else SOURCE.replace(before, after)
        cases.append({"id": case, "source": text, "expected_refusal": refuse,
                      "sha256": hashlib.sha256(text.encode()).hexdigest()})
    return cases


def snapshot(text):
    # ElementTree supplies independent XML 1.0 character-reference decoding.
    root = ET.fromstring(text)
    atoms, edges, ids = [], [], {}
    for node in root.iter():
        if node.tag.split('}')[-1] != 'atom':
            continue
        a = node.attrib
        require(a['id'] not in ids, 'ambiguous decoded atom identity')
        ids[a['id']] = len(atoms)
        atoms.append({'element': a['elementType'], 'charge': int(a.get('formalCharge', 0)),
                      'isotope': int(a.get('isotopeNumber', a.get('isotope', 0))),
                      'hydrogens': int(a['hydrogenCount']) if 'hydrogenCount' in a else None,
                      'coordinates': [float(a.get('x2', 0)), float(a.get('y2', 0))]})
    for node in root.iter():
        if node.tag.split('}')[-1] == 'bond':
            a, b = node.attrib['atomRefs2'].split()
            order = {'S': 1, 'D': 2, 'T': 3}.get(node.attrib['order'])
            edges.append([*sorted((ids[a], ids[b])), order or int(node.attrib['order'])])
    return {'atoms': atoms, 'bonds': sorted(edges)}


def failed_process(engine, run, expected_refusal):
    # A panic, signal or runtime/launch failure is never a semantic rejection.
    code, stderr = run['exit_code'], run['stderr']
    if isinstance(stderr, bytes): stderr = stderr.decode('utf8', errors='replace')
    if code not in (0, 1): return 'process_failure'
    rejected = ((engine == 'chematic' and code == 1 and stderr.startswith('Error:'))
                or (engine == 'openbabel' and '0 molecules converted' in stderr))
    return 'rejected' if rejected and expected_refusal else 'read_write_failure'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--obabel', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), 'output already exists')
    version = bounded_run([str(args.obabel), '-V'])
    require(version['exit_code'] == 0 and b'Open Babel 3.2.1' in version['stdout'], 'expected Open Babel 3.2.1')
    cases, rows = packet(), []
    expected = snapshot(SOURCE)
    with tempfile.TemporaryDirectory(prefix='chematic-cml-ref-') as directory:
        for case in cases:
            require(len(case['source'].encode()) < 16384, 'packet input too large')
            path, output = Path(directory)/'input.cml', Path(directory)/'output.cml'
            path.write_text(case['source'])
            for engine, command in [
                ('chematic', [str(args.probe), '--format', 'cml', '--path', str(path), '--rewrite-output', str(output)]),
                ('openbabel', [str(args.obabel), '-icml', '-ocml']),
            ]:
                output.unlink(missing_ok=True)
                run = bounded_run(command, case['source'].encode(), timeout=10, max_output=262144)
                row = {'case': case['id'], 'engine': engine, 'command': command,
                       **{k: v.decode('utf8', errors='replace') if isinstance(v, bytes) else v for k,v in run.items()}}
                row['terminal'] = 'process_failure'
                if run['resource_status'] != 'completed':
                    row['terminal'] = run['resource_status']
                elif run['exit_code'] != 0 or (engine == 'openbabel' and not run['stdout'].strip()):
                    row['terminal'] = failed_process(engine, run, case['expected_refusal'])
                else:
                    try:
                        require(engine != 'chematic' or output.stat().st_size <= 262144, 'output too large')
                        text = output.read_text() if engine == 'chematic' else row['stdout']
                        row['output'] = text
                        row['snapshot'] = snapshot(text)
                        row['fields'] = {k: row['snapshot'][k] == v for k,v in expected.items()}
                        row['terminal'] = ('output_on_ambiguous_or_invalid_source' if case['expected_refusal']
                                           else 'packet_fields_preserved' if all(row['fields'].values()) else 'packet_field_loss')
                    except (ValueError, KeyError, IndexError, ET.ParseError) as error:
                        row['terminal'], row['error'] = 'unresolved_output', str(error)
                rows.append(row)
    result = {'schema':'chematic.cml-character-reference-development.v1',
              'scope':'Original CML→CML development packet: indexed elements, charges, isotopes, total H counts, 2D coordinates and explicit bonds only. No full XML/CML conformance, source metadata, stereo, 3D, independent holdout or toolkit-wide superiority claim. Rejections are not successful coverage.',
              'provenance':'Original minimal ammonium/deuterium structures, 2026-10-22, recipe packet(); no third-party or sealed corpus input.',
              'license':'MIT OR Apache-2.0','redistribution':'allowed: original inputs',
              'minimization':'Two atoms isolate XML lexical changes; an extra duplicate ID tests ambiguity. Positions are supplied angstrom coordinates, not optimized geometries.',
              'oracle':{'name':'Python ElementTree plus authored typed-field expectations','version':platform.python_version(),'script':identity(__file__),'expected':expected,'tolerance':'exact integers and supplied dyadic 2D floats; no tolerance'},
              'engines':{'chematic':identity(args.probe),'openbabel':identity(args.obabel)},
              'openbabel_version':version['stdout'].decode().strip(), 'cases':cases,'rows':rows,
              'counts':{e:dict(Counter(r['terminal'] for r in rows if r['engine']==e)) for e in ['chematic','openbabel']}}
    with args.output.open('x') as output:
        output.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['counts'],sort_keys=True))


if __name__ == '__main__':
    main()
