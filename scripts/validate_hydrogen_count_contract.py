#!/usr/bin/env python3
"""Check hydrogen accounting against pinned RDKit; diagnose Open Babel reports.

Provide fresh chematic and RDKit 2026.09.1 on PYTHONPATH. Optional
--wasm-module checks total and per-atom attached H in real Node/WASM;
--baseline-pythonpath probes an older source runtime in a separate process.
"""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import importlib
import json
import os
from pathlib import Path
import random
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / 'validation/openbabel_issues/hydrogen-counts/cases.tsv'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def formula_h(formula):
    # Open Babel may spell hydrogen isotopes D/T; aggregate nuclei, not masses.
    return sum(int(n or 1) for symbol,n in re.findall(r'([A-Z][a-z]?|\*)(\d*)',formula) if symbol in ['H','D','T'])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--wasm-module',type=Path)
    p.add_argument('--baseline-pythonpath',type=Path)
    p.add_argument('--obabel',default='/opt/homebrew/bin/obabel')
    args = p.parse_args()
    import chematic
    from rdkit import Chem, rdBase
    assert rdBase.rdkitVersion == '2026.09.1'
    def reference(smiles):
        params=Chem.SmilesParserParams();params.removeHs=False;params.sanitize=False
        m=Chem.MolFromSmiles(smiles,params)
        assert m is not None,smiles
        m.UpdatePropertyCache(strict=False)
        return m
    originals=[]
    for line in FIXTURE.read_text().splitlines():
        if not line.startswith('#'):
            smiles,h=line.split('\t')
            originals.append(dict(smiles=smiles,expected=int(h),kind='original'))
    rows=list(originals)
    bases=['C','CC','CCC','CCCC','C1CCCCC1','CCO','CC(=O)O','N','CN','[NH4+]','C[NH3+]','O','[OH-]','c1ccccc1','c1ccncc1','[nH]1cccc1','CC(=O)N','CC(C)O','[2H]O[3H]','N[C@@H](C)O']
    for base in bases:
        params=Chem.SmilesParserParams();params.removeHs=False
        m=Chem.MolFromSmiles(base,params)
        assert m is not None
        for style in ['compact','bracket','explicit_nodes']:
            source=Chem.AddHs(m) if style=='explicit_nodes' else m
            for seed in range(3):
                order=list(range(source.GetNumAtoms()));random.Random(seed).shuffle(order)
                permuted=Chem.RenumberAtoms(source,order)
                smiles=Chem.MolToSmiles(permuted,canonical=False,allHsExplicit=style=='bracket')
                rows.append(dict(smiles=smiles,kind='generated',base=base,style=style,seed=seed))
    components=['[OH2]','[H][O][H]','[NH4+]','[Cl-]','[2H]','[H-]','[CH4]','[He]']
    for a in components:
        for b in components:
            rows.append(dict(smiles=a+'.'+b,kind='component_pair'))
    for row in rows:
        ref=reference(row['smiles'])
        expected=sum(int(a.GetAtomicNum()==1)+a.GetTotalNumHs() for a in ref.GetAtoms())
        if 'expected' in row:
            assert row['expected']==expected,row
        row['expected']=expected
        row['attached']=[a.GetTotalNumHs(includeNeighbors=True) for a in ref.GetAtoms()]
        m=chematic.from_smiles(row['smiles'])
        assert m.num_hydrogens==expected,row
        assert chematic.from_smiles(m.smiles).num_hydrogens==expected,row
        assert formula_h(m.formula)==expected,row
        row.update(candidate=m.num_hydrogens,candidate_smiles=m.smiles,candidate_formula=m.formula)
    wasm_hashes=None
    if args.wasm_module:
        code="""const fs=require('fs'),w=require(process.argv[1]);const rows=JSON.parse(fs.readFileSync(0,'utf8'));const out=rows.map(r=>{const m=w.parse_smiles(r.smiles),v={total:w.num_hydrogens(m),attached:Array.from({length:m.atom_count()},(_,i)=>JSON.parse(w.get_atom_info(m,i)).totalHydrogens)};m.free();return v});console.log(JSON.stringify(out));"""
        actual=json.loads(subprocess.check_output(['node','-e',code,str(args.wasm_module.resolve())],input=json.dumps(rows),text=True,timeout=60))
        assert len(actual)==len(rows)
        for row,got in zip(rows,actual):
            assert got['total']==row['expected'],row
            assert got['attached']==row['attached'],(row,got)
            row['wasm']=got
        wasm_hashes={f.name:sha(f) for f in args.wasm_module.parent.iterdir() if f.is_file()}
    baseline=None
    if args.baseline_pythonpath:
        code="""import chematic,json,sys
rows=json.load(sys.stdin);out=[]
for r in rows:
 m=chematic.from_smiles(r['smiles']);out.append({'total':m.num_hydrogens,'formula':m.formula})
print(json.dumps(out))
"""
        old=json.loads(subprocess.check_output([sys.executable,'-c',code],input=json.dumps(rows),text=True,env=dict(os.environ,PYTHONPATH=str(args.baseline_pythonpath)),timeout=30))
        assert len(old)==len(rows)
        formula_mismatches=sum(formula_h(v['formula'])!=r['expected'] for v,r in zip(old,rows))
        assert formula_mismatches==0
        baseline=dict(scope='Retained pre-fix source runtime, not a released package',extension_sha256=sha(args.baseline_pythonpath/'chematic/chematic.so'),rows=old,mismatch_count=sum(v['total']!=r['expected'] for v,r in zip(old,rows)),formula_h_mismatch_count=formula_mismatches)
    ob=[]
    for row in originals:
        r=subprocess.run([args.obabel,'-ismi','-ocan','--append','formula'],input=row['smiles']+'\n',text=True,capture_output=True,timeout=10)
        record=dict(input=row['smiles'],expected=row['expected'],returncode=r.returncode,stdout=r.stdout,stderr=r.stderr)
        line=r.stdout.rstrip('\r\n')
        if r.returncode==0 and '\t' in line:
            smiles,formula=line.split('\t',1)
            record['formula']=formula
            if formula:
                record['formula_h']=formula_h(formula)
            else:
                record['formula_status']='empty formula returned'
            parsed=reference(smiles)
            record['output_smiles_h']=sum(int(a.GetAtomicNum()==1)+a.GetTotalNumHs() for a in parsed.GetAtoms())
        ob.append(record)
    out=ROOT/'validation/results/openbabel-hydrogen-counts-2026-10-16.json'
    out.parent.mkdir(parents=True,exist_ok=True)
    packet=out.with_suffix('.packet.json.gz')
    with packet.open('wb') as f:
        with gzip.GzipFile(fileobj=f,mode='wb',mtime=0) as z:
            z.write((json.dumps(dict(cases=rows,baseline=baseline,openbabel=ob),separators=(',',':'))+'\n').encode())
    compared=[r for r in ob if 'formula_h' in r]
    serialized=[r for r in ob if 'output_smiles_h' in r]
    sources=['crates/chematic-chem/src/descriptors/composition.rs','crates/chematic-wasm/src/mol_depict.rs','crates/chematic-wasm/src/mol_descriptors.rs']
    result=dict(schema='chematic.hydrogen-count-validation.v1',utc=datetime.now(timezone.utc).isoformat(),
                issues=['https://github.com/openbabel/openbabel/issues/3008','https://github.com/openbabel/openbabel/issues/2624'],
                source_base=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),source_sha256={f:sha(ROOT/f) for f in sources},
                runner_sha256=sha(__file__),fixture_sha256=sha(FIXTURE),rdkit_version=rdBase.rdkitVersion,rdkit_runtime_sha256=sha(rdBase.__file__),
                python_extension_sha256=sha(importlib.import_module('chematic.chematic').__file__),wasm_runtime_sha256=wasm_hashes,
                cases=len(rows),distinct_smiles=len({r['smiles'] for r in rows}),python_matches=len(rows),python_smiles_roundtrip_matches=len(rows),
                python_formula_h_controls=len(rows),wasm_total_matches=len(rows) if wasm_hashes else None,wasm_atom_info_matches=sum(len(r['attached']) for r in rows) if wasm_hashes else None,candidate_failures=0,
                baseline=None if baseline is None else {k:v for k,v in baseline.items() if k!='rows'},
                openbabel_version=subprocess.check_output([args.obabel,'-V'],text=True).strip(),openbabel_binary_sha256=sha(args.obabel),
                openbabel=dict(formula_compared=len(compared),formula_missing=len(ob)-len(compared),
                               missing_formula_inputs=[r['input'] for r in ob if 'formula_h' not in r],
                               smiles_compared=len(serialized),formula_h_mismatches=sum(r['formula_h']!=r['expected'] for r in compared),
                               smiles_h_mismatches=sum(r['output_smiles_h']!=r['expected'] for r in serialized)),
                packet=packet.name,packet_sha256=sha(packet),scope='Hydrogen nuclei and attached-H accounting only. Native formulas were already correct on these cases. Unsanitized RDKit graph/cache references retain unusual explicit-H inputs; this is not chemical-validity, mass or general engine-parity evidence. Native H transformation internals and SMILES serialization are unchanged; source candidate only.')
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['cases','distinct_smiles','python_matches','wasm_total_matches','wasm_atom_info_matches','candidate_failures','openbabel']}))


if __name__=='__main__':
    main()
