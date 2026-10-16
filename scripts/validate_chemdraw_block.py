#!/usr/bin/env python3
"""Native binary CDX/UTF-8 CDXML dispatch against independent RDKit XML reads."""
import argparse
import gzip
import hashlib
import json
import platform
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

import chematic
import chematic.chematic as extension
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdChemDraw

ROOT=Path(__file__).resolve().parents[1]
FIXTURE=ROOT/'validation/rdkit_issues/chemdraw-block/cases.json'

def sha(b):return hashlib.sha256(b).hexdigest()

def graph(m):
    assert m is not None
    return dict(atoms=[(a.GetAtomicNum(),a.GetFormalCharge(),a.GetIsotope(),a.GetTotalNumHs()) for a in m.GetAtoms()],bonds=sorted((min(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),max(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),b.GetBondTypeAsDouble()) for b in m.GetBonds()),smiles=Chem.MolToSmiles(m))

def graph_block(block):return graph(Chem.MolFromMolBlock(block,sanitize=True,removeHs=False))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--wasm-module',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    assert rdBase.rdkitVersion=='2026.09.1';RDLogger.DisableLog('rdApp.*')
    cases=json.loads(FIXTURE.read_text(encoding='utf-8'))['cases']
    jobs=[]
    for c in cases:
        for mode in ('AUTO','CDX'):
            jobs.append(dict(id=c['id'],input_hex=c['cdx_hex'],format=mode,allow_metadata_loss=c['allow_metadata_loss']))
        for mode in ('AUTO','CDXML'):
            jobs.append(dict(id=c['id'],input_hex=c['cdxml'].encode().hex(),format=mode,allow_metadata_loss=False))
    node=subprocess.run(['node',str(ROOT/'scripts/validate_chemdraw_block_node.mjs'),str(a.wasm_module.resolve())],input=json.dumps(jobs),text=True,capture_output=True,check=True)
    node_rows=json.loads(node.stdout);assert len(node_rows)==len(jobs)
    by_id={c['id']:c for c in cases};raw=[]
    for job,js in zip(jobs,node_rows):
        c=by_id[job['id']];source=bytes.fromhex(job['input_hex'])
        f=chematic.from_chemdraw_block(source,format=job['format'],allow_metadata_loss=job['allow_metadata_loss'])
        reference=[graph(m) for m in Chem.MolsFromCDXML(c['cdxml'],sanitize=True,removeHs=False)]
        native=[graph_block(m.to_mol_v3000([])) for m in f['molecules']]
        wasm=[graph_block(b) for b in js['blocks']]
        assert native==reference,(job,native,reference)
        assert wasm==reference,(job,wasm,reference)
        expected_coords=[[[float(v) for v in n.attrib['p'].split()] for n in frag.findall('n')] for frag in ET.fromstring(c['cdxml']).iter('fragment')]
        assert f['coordinates']==expected_coords==js['report']['coordinates']
        assert f['source_bytes']==source==bytes(js['report']['source_bytes'])
        for key in ('format','projected_cdxml','omissions'):assert f[key]==js['report'][key]
        raw.append(dict(job,rdkit=reference,python=native,python_coordinates=f['coordinates'],python_projection=f['projected_cdxml'],python_omissions=f['omissions'],wasm=js))
    probe=bytes.fromhex(cases[0]['cdx_hex'])
    block_result=rdChemDraw.MolsFromChemDrawBlock(probe)
    assert len(block_result)==0
    with tempfile.TemporaryDirectory(prefix='chematic-chemdraw-oracle-') as directory:
        file=Path(directory)/'ethanol.cdx';file.write_bytes(probe)
        try:rdChemDraw.MolsFromChemDrawFile(str(file))
        except RuntimeError as e:file_error=str(e)
        else:raise AssertionError('expected pinned RDKit XML dispatch failure')
    try:rdChemDraw.MolsFromChemDrawBlock(probe,format=rdChemDraw.CDXFormat.CDX)
    except TypeError as e:explicit_error=str(e)
    else:raise AssertionError('expected missing explicit format argument')
    reproduction=dict(url='https://github.com/rdkit/rdkit/issues/9604',input_sha256=sha(probe),binary_block_molecule_count=0,binary_file_error=file_error,explicit_format_error=explicit_error,native_binary_molecule_count=1)
    a.output.mkdir(parents=True,exist_ok=True)
    packed=gzip.compress(json.dumps(dict(comparisons=raw,rdkit_issue=reproduction),sort_keys=True).encode(),mtime=0);(a.output/'comparisons.json.gz').write_bytes(packed)
    sources=['crates/chematic-mol/src/chemdraw.rs','crates/chematic-mol/src/chemdraw/binary.rs','crates/chematic-mol/src/chemdraw/properties.rs','crates/chematic-mol/src/chemdraw/xml.rs','crates/chematic-py/src/formats/chemdraw.rs','crates/chematic-wasm/src/chemdraw.rs','scripts/generate_chemdraw_fixtures.py','scripts/validate_chemdraw_block.py','scripts/validate_chemdraw_block_node.mjs','validation/rdkit_issues/chemdraw-block/cases.json']
    result=dict(schema='chematic.chemdraw-validation.v1',source_base=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),source_sha256={s:sha((ROOT/s).read_bytes()) for s in sources},runtime=dict(rdkit=rdBase.rdkitVersion,python=platform.python_version(),platform=platform.platform(),python_extension_sha256=sha(Path(extension.__file__).read_bytes()),wasm_sha256=sha(a.wasm_module.with_name('chematic_wasm_bg.wasm').read_bytes())),results=dict(authored_binary_recipes=len(cases),python_graph_comparisons=len(jobs),wasm_graph_comparisons=len(jobs),source_bytes_and_coordinate_checks=len(jobs),failures=0),issue=reproduction,packet=dict(path='comparisons.json.gz',sha256=sha(packed)),scope='bounded concrete CDX profile and explicit dispatch; CDXML graph policy retained. No full ChemDraw, reaction/query, marked CDX stereo, vendor round trip, speed or published-artifact claim')
    (a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
