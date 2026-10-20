#!/usr/bin/env python3
"""Bounded RDKit #9359 check: unique-center Bold/Hash/Dash inference.

Provide a fresh chematic extension and RDKit 2026.09.1 on PYTHONPATH.
--wasm-module checks real Node/WASM. --baseline-pythonpath diagnoses an
older source runtime in its own process, not a released wheel.
"""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import random
import subprocess
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / 'validation/rdkit_issues/cdxml/nondirectional-stereo.cdxml'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def cases():
    result = []
    for display in ['Bold', 'Hash', 'Dash']:
        for reverse in [False, True]:
            for angle in [0, 30, 90, 170, 240, 315]:
                for mirror in [False, True]:
                    for seed in range(3):
                        tree = ET.fromstring(FIXTURE.read_text())
                        fragment = tree.find('.//fragment')
                        nodes = list(fragment.findall('n'))
                        s, c = math.sin(math.radians(angle)), math.cos(math.radians(angle))
                        for node in nodes:
                            x, y = map(float, node.get('p').split())
                            if mirror:
                                y = -y
                            node.set('p', f'{(x*c-y*s)*2+12:.12f} {(x*s+y*c)*2-5:.12f}')
                            fragment.remove(node)
                        random.Random(seed).shuffle(nodes)
                        for i, node in enumerate(nodes):
                            fragment.insert(i, node)
                        bond = fragment.find("b[@id='101']")
                        bond.set('Display', display)
                        if reverse:
                            bond.set('B', '2'); bond.set('E', '1')
                        text = ET.tostring(tree, encoding='unicode')
                        # A unique eligible carbon endpoint is the intended narrow end.
                        # Use RDKit's unambiguous directional parser as independent control.
                        bond.set('B', '1'); bond.set('E', '2')
                        bond.set('Display', 'WedgeBegin' if display == 'Bold' else 'WedgedHashBegin')
                        result.append(dict(display=display,reverse=reverse,angle=angle,mirror=mirror,seed=seed,
                                           input=text,reference=ET.tostring(tree,encoding='unicode')))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--wasm-module', type=Path)
    p.add_argument('--baseline-pythonpath', type=Path)
    args = p.parse_args()
    import chematic
    from rdkit import Chem, rdBase
    assert rdBase.rdkitVersion == '2026.09.1'
    rows = cases()
    def canonical(smiles):
        m = Chem.MolFromSmiles(smiles)
        assert m is not None
        return Chem.MolToSmiles(m, isomericSmiles=True)
    for row in rows:
        reference = Chem.MolsFromCDXML(row['reference'])[0]
        expected = Chem.MolToSmiles(reference, isomericSmiles=True)
        expected_codes = sorted(a.GetProp('_CIPCode') for a in reference.GetAtoms() if a.HasProp('_CIPCode'))
        assert len(expected_codes) == 1
        mol = chematic.from_cdxml(row['input'], infer_nondirectional_stereo=True)
        assert canonical(mol.smiles) == expected, row
        codes = sorted(v['descriptor'] for v in mol.cip_stereo())
        assert codes == expected_codes
        default = chematic.from_cdxml(row['input'])
        upstream = Chem.MolsFromCDXML(row['input'])[0]
        assert not default.cip_stereo()
        assert not any(a.HasProp('_CIPCode') for a in upstream.GetAtoms())
        assert canonical(default.smiles) == Chem.MolToSmiles(upstream)
        row.update(expected=expected,expected_codes=expected_codes,candidate=mol.smiles,
                   candidate_codes=codes,default=default.smiles,upstream_default=Chem.MolToSmiles(upstream))
    directional = []
    for row in rows:
        if row['display'] != 'Bold':
            continue
        for display in ['WedgeBegin', 'WedgeEnd', 'WedgedHashBegin', 'WedgedHashEnd']:
            tree = ET.fromstring(row['input'])
            tree.find(".//b[@id='101']").set('Display', display)
            text = ET.tostring(tree, encoding='unicode')
            ref = Chem.MolsFromCDXML(text)[0]
            expected = Chem.MolToSmiles(ref)
            expected_codes = sorted(a.GetProp('_CIPCode') for a in ref.GetAtoms() if a.HasProp('_CIPCode'))
            default = chematic.from_cdxml(text)
            enabled = chematic.from_cdxml(text, infer_nondirectional_stereo=True)
            assert canonical(default.smiles) == canonical(enabled.smiles) == expected
            assert sorted(v['descriptor'] for v in default.cip_stereo()) == expected_codes
            directional.append(dict(input=text,display=display,reverse=row['reverse'],angle=row['angle'],
                                    mirror=row['mirror'],seed=row['seed'],expected=expected,
                                    expected_codes=expected_codes,candidate=enabled.smiles,default=default.smiles))
    all_rows = rows + directional
    wasm_hashes = None
    if args.wasm_module:
        js = """const fs=require('fs'),w=require(process.argv[1]);const rows=JSON.parse(fs.readFileSync(0,'utf8'));const out=rows.map(r=>{const m=w.mol_from_cdxml_with_options(r.input,true),d=w.mol_from_cdxml(r.input);const v={smiles:m.canonical_smiles(),codes:JSON.parse(w.cip_assignments_json(m)).map(v=>v.cipCode).sort(),default:d.canonical_smiles(),batch:JSON.parse(w.cdxml_to_smiles_json_with_options(r.input,true))};m.free();d.free();return v});console.log(JSON.stringify(out));"""
        actual = json.loads(subprocess.check_output(['node','-e',js,str(args.wasm_module.resolve())],input=json.dumps(all_rows),text=True,timeout=60))
        assert len(actual) == len(all_rows)
        for row, got in zip(all_rows, actual):
            assert canonical(got['smiles']) == row['expected']
            assert got['codes'] == row['expected_codes']
            assert canonical(got['default']) == canonical(row['default'])
            assert len(got['batch']) == 1 and canonical(got['batch'][0]) == row['expected']
            row['wasm'] = got
        wasm_hashes = {f.name:sha(f) for f in args.wasm_module.parent.iterdir() if f.is_file()}
    baseline = None
    if args.baseline_pythonpath:
        code = """import chematic,json,sys
data=json.load(sys.stdin);s=data['fixture'];m=chematic.from_cdxml(s)
try: chematic.from_cdxml(s,infer_nondirectional_stereo=True);has=True
except TypeError: has=False
print(json.dumps({'option_available':has,'default':m.smiles,'directional':[chematic.from_cdxml(v['input']).smiles for v in data['controls']]}))
"""
        baseline = json.loads(subprocess.check_output([sys.executable,'-c',code],input=json.dumps({'fixture':FIXTURE.read_text(),'controls':directional}),text=True,env=dict(os.environ,PYTHONPATH=str(args.baseline_pythonpath)),timeout=30))
        baseline['directional_mismatch_count'] = sum(canonical(got) != row['expected'] for got,row in zip(baseline['directional'],directional))
        baseline.update(scope='Retained pre-fix source runtime, not a released package',extension_sha256=sha(args.baseline_pythonpath/'chematic/chematic.so'))
    out = ROOT / 'validation/results/rdkit-cdxml-stereo-options-2026-10-16.json'
    out.parent.mkdir(parents=True,exist_ok=True)
    packet = out.with_suffix('.packet.json.gz')
    with packet.open('wb') as f:
        with gzip.GzipFile(fileobj=f,mode='wb',mtime=0) as z:
            z.write((json.dumps({'nondirectional':rows,'directional':directional,'baseline':baseline},separators=(',',':'))+'\n').encode())
    sources = ['crates/chematic-mol/src/cdxml.rs','crates/chematic-py/src/formats.rs','crates/chematic-wasm/src/mol_io.rs']
    result = dict(schema='chematic.cdxml-stereo-options-validation.v1',utc=datetime.now(timezone.utc).isoformat(),
                  issue='https://github.com/rdkit/rdkit/issues/9359',source_base=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
                  source_sha256={f:sha(ROOT/f) for f in sources},runner_sha256=sha(__file__),fixture_sha256=sha(FIXTURE),
                  rdkit_version=rdBase.rdkitVersion,rdkit_runtime_sha256=sha(rdBase.__file__),
                  python_extension_sha256=sha(importlib.import_module('chematic.chematic').__file__),wasm_runtime_sha256=wasm_hashes,
                  corpus_cases=len(rows),python_configurations_matched=len(rows),python_default_controls_matched=len(rows),
                  wasm_configurations_matched=len(rows) if wasm_hashes else None,wasm_batch_configurations_matched=len(rows) if wasm_hashes else None,
                  wasm_default_controls_matched=len(rows) if wasm_hashes else None,upstream_nondirectional_unassigned=len(rows),
                  python_directional_controls_matched=len(directional),wasm_directional_controls_matched=len(directional) if wasm_hashes else None,
                  candidate_failures=0,baseline=None if baseline is None else {k:v for k,v in baseline.items() if k!='directional'},packet=packet.name,packet_sha256=sha(packet),
                  scope='Unique degree-4 carbon / terminal F drawings; CIP and full isomeric graph comparisons against explicit directional controls. Transformations, B/E reversal and node-order changes. Heuristic inference, not general ChemDraw compatibility or source preservation; source candidate only.')
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k.endswith('matched') or k in ['corpus_cases','candidate_failures']}))


if __name__ == '__main__':
    main()
