#!/usr/bin/env python3
"""Original ViewMol geometry packet: actual bindings/CLI and Open Babel rewrite.

Observations concern local same-format workflows, not full ViewMol support,
isolated parser/writer accuracy, released packages, speed or broad superiority.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
ROOT=Path(__file__).resolve().parents[2]
PACKET=ROOT/'validation/viewmol_preservation/cases.json'
NODE=r'''
const fs=require('node:fs'),w=require(process.argv[1]),p=JSON.parse(fs.readFileSync(0,'utf8'));
const out={valid:[],refused:[],projection_refused:[]};
const status=e=>!(e instanceof WebAssembly.RuntimeError)&&String(e).startsWith('ViewMol line ')?'refused':'error';
for(const c of p.valid){
 try{let s=w.parse_viewmol_json(c.text);const states=[JSON.parse(s)],projections_json=[w.viewmol_positions_angstrom_json(s)];
  for(let i=0;i<3;i++){s=w.parse_viewmol_json(w.write_viewmol_from_json(s));states.push(JSON.parse(s));projections_json.push(w.viewmol_positions_angstrom_json(s));}
  out.valid.push({id:c.id,states,projections_json});
 }catch(e){out.valid.push({id:c.id,error:String(e)});}
}
for(const c of p.refused){try{out.refused.push({id:c.id,status:'accepted',output:w.parse_viewmol_json(c.text)});}catch(e){out.refused.push({id:c.id,status:status(e),error:String(e)});}}
for(const c of p.projection_refused){
 const s=w.parse_viewmol_json(w.write_viewmol_from_json(w.parse_viewmol_json(c.text))),r={id:c.id,state:JSON.parse(s)};
 try{r.output=w.viewmol_positions_angstrom_json(s);r.status='accepted';}catch(e){r.status=status(e);r.error=String(e);}out.projection_refused.push(r);
}
process.stdout.write(JSON.stringify(out));
'''
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def run(args,text=None):
 r=subprocess.run(args,input=text,text=True,capture_output=True,timeout=30)
 return {'argv':[str(a) for a in args],'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr}
def cli_refusal(r):return r['returncode']==2 and not r['stdout'] and r['stderr'].startswith('chematic: ViewMol line ')
def format_refusal(e):return isinstance(e,ValueError) and str(e).startswith('ViewMol line ')
def bits(rows):return [[float(v).hex() for v in row] for row in rows]
def projected(record):
 factor=float(record['factor_to_angstrom'])
 result=[[float(v)*factor for v in a['position']] for a in record['atoms']]
 if any(not math.isfinite(v) for row in result for v in row):raise ValueError('nonfinite projected external coordinate')
 return result

def observe(text):
 """Independent extraction of the external writer's data groups."""
 lines=[line.removesuffix("\r") for line in text.split("\n")]
 groups=[(i,line.split()) for i,line in enumerate(lines) if line.lstrip().startswith('$')]
 coords=[(i,f) for i,f in groups if f and f[0]=='$coord']
 titles=[i for i,f in groups if f==['$title']]
 if len(coords)!=1 or len(titles)>1 or not groups or groups[-1][1]!=['$end']:
  raise ValueError('missing/duplicate external groups')
 if any(f[0] not in ('$coord','$title','$end') for _,f in groups):raise ValueError('unexpected external group')
 index,header=coords[0]
 if len(header)!=2:raise ValueError('missing external scale')
 end=next(i for i,_ in groups if i>index)
 rows=[line.split() for line in lines[index+1:end] if line.strip()]
 if not rows or any(len(row) not in (4,5) for row in rows):raise ValueError('unexpected external atom section')
 title=None
 if titles:
  i=titles[0]
  if i+1>=len(lines) or lines[i+1].lstrip().startswith('$'):raise ValueError('missing external title')
  title=lines[i+1]
 record={'title':title,'factor_to_angstrom':header[1],
         'atoms':[{'element':r[3],'position':r[:3],'fixed_axes':r[4] if len(r)==5 else None} for r in rows]}
 projected(record)
 return record

def comparisons(actual,expected):
 result={'full_token_record_exact':actual==expected,
         'raw_scale_coordinates_exact':actual['factor_to_angstrom']==expected['factor_to_angstrom'] and [a['position'] for a in actual['atoms']]==[a['position'] for a in expected['atoms']],
         'angstrom_positions_bits_exact':bits(projected(actual))==bits(projected(expected)),
         'element_order_exact':[a['element'] for a in actual['atoms']]==[a['element'] for a in expected['atoms']],
         'fixed_axes_exact':[sorted(a['fixed_axes'] or '') for a in actual['atoms']]==[sorted(a['fixed_axes'] or '') for a in expected['atoms']],
         'title_exact':actual['title']==expected['title']}
 result['declared_state_exact']=all(result[k] for k in ['angstrom_positions_bits_exact','element_order_exact','fixed_axes_exact','title_exact'])
 return result

def python_lane(p,m):
 out={'valid':[],'refused':[],'projection_refused':[]}
 for c in p['valid']:
  r=m.parse_viewmol(c['text']);states=[r];projections=[bits(m.viewmol_positions_angstrom(r))]
  for _ in range(3):r=m.parse_viewmol(m.write_viewmol(r));states.append(r);projections.append(bits(m.viewmol_positions_angstrom(r)))
  out['valid'].append({'id':c['id'],'states':states,'projection_bits':projections})
 for c in p['refused']:
  try:out['refused'].append({'id':c['id'],'status':'accepted','output':m.parse_viewmol(c['text'])})
  except Exception as e:out['refused'].append({'id':c['id'],'status':'refused' if format_refusal(e) else 'error','error':repr(e)})
 for c in p['projection_refused']:
  r=m.parse_viewmol(m.write_viewmol(m.parse_viewmol(c['text'])));row={'id':c['id'],'state':r}
  try:row.update(status='accepted',output=m.viewmol_positions_angstrom(r))
  except Exception as e:row.update(status='refused' if format_refusal(e) else 'error',error=repr(e))
  out['projection_refused'].append(row)
 return out

def cli_lane(p,cli):
 out={'valid':[],'refused':[],'projection_refused':[]}
 def checked(text,args):
  r=run([cli,'viewmol',*args],text)
  if r['returncode']:raise RuntimeError(r)
  return r['stdout']
 for c in p['valid']:
  r=json.loads(checked(c['text'],['--inspect']));states=[r];projections=[bits(json.loads(checked(json.dumps(r),['--from-json','--angstrom'])))]
  for _ in range(3):
   text=checked(json.dumps(r),['--from-json']);r=json.loads(checked(text,['--inspect']));states.append(r);projections.append(bits(json.loads(checked(json.dumps(r),['--from-json','--angstrom']))))
  out['valid'].append({'id':c['id'],'states':states,'projection_bits':projections})
 for c in p['refused']:
  r=run([cli,'viewmol'],c['text']);out['refused'].append({'id':c['id'],'status':'refused' if cli_refusal(r) else 'error','receipt':r})
 for c in p['projection_refused']:
  text=checked(c['text'],[]);state=json.loads(checked(text,['--inspect']));r=run([cli,'viewmol','--angstrom'],text)
  out['projection_refused'].append({'id':c['id'],'state':state,'status':'refused' if cli_refusal(r) else 'error','receipt':r})
 return out

def summarize(lane,p):
 expected={c['id']:c for c in p['valid']};unsafe={c['id']:c['expected'] for c in p['projection_refused']}
 exact=projection=0
 for r in lane['valid']:
  c=expected[r['id']]
  exact+=len(r.get('states',[]))==4 and all(s==c['expected'] for s in r['states'])
  values=r.get('projection_bits')
  if values is None:values=[bits(json.loads(text)) for text in r.get('projections_json',[])]
  projection+=len(values)==4 and all(v==bits(c['expected_angstrom']) for v in values)
 return {'valid_total':len(p['valid']),'exact_all_four_states':exact,'projected_bits_all_four':projection,
         'refusal_total':len(p['refused']),'typed_refused':sum(r['status']=='refused' for r in lane['refused']),
         'unsafe_projection_total':len(unsafe),'raw_projection_records_exact':sum(r['state']==unsafe[r['id']] for r in lane['projection_refused']),
         'unsafe_projection_refused':sum(r['status']=='refused' for r in lane['projection_refused'])}

def openbabel_lane(p,binary):
 rows=[]
 with tempfile.TemporaryDirectory(prefix='chematic-viewmol-') as directory:
  for c in p['valid']:
   source=Path(directory)/(c['id']+'.vmol');output=Path(directory)/(c['id']+'.out.vmol');source.write_text(c['text'])
   receipt=run([binary,'-ivmol',source,'-ovmol','-O',output,'-ab']);text=output.read_text() if output.exists() else ''
   row={'id':c['id'],'receipt':receipt,'output':text,'source_has_fixed_axes':any(a['fixed_axes'] for a in c['expected']['atoms'])}
   try:
    if receipt['returncode']:raise ValueError('nonzero external exit')
    actual=observe(text);row.update(status='observed',state=actual,comparisons=comparisons(actual,c['expected']))
   except (ValueError,IndexError,KeyError) as e:row.update(status='unobservable',error=str(e))
   rows.append(row)
 return rows

def unit_probe(binary):
 text='$title\nunit probe\n$coord 2.0\n.125 .25 .5 C\n$end\n'
 with tempfile.TemporaryDirectory(prefix='chematic-viewmol-unit-') as directory:
  source=Path(directory)/'units.vmol';source.write_text(text);receipt=run([binary,'-ivmol',source,'-oxyz','-ab'])
  try:positions=[float(v) for v in receipt['stdout'].splitlines()[2].split()[1:]]
  except (ValueError,IndexError):positions=None
  return {'source':text,'receipt':receipt,'expected_angstrom':[.25,.5,1.0],'observed':positions,'exact':receipt['returncode']==0 and positions==[.25,.5,1.]}

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 for name in ['python-path','wasm-module','cli','output']:parser.add_argument('--'+name,type=Path,required=True)
 parser.add_argument('--openbabel',default='/opt/homebrew/bin/obabel');a=parser.parse_args();p=json.loads(PACKET.read_text())
 sys.path.insert(0,str(a.python_path.resolve()));m=importlib.import_module('chematic')
 node=run(['node','-e',NODE,str(a.wasm_module.resolve())],json.dumps(p))
 if node['returncode']:raise RuntimeError(node)
 lanes={'python':python_lane(p,m),'wasm_node':json.loads(node['stdout']),'cli':cli_lane(p,str(a.cli.resolve()))};summary={k:summarize(v,p) for k,v in lanes.items()}
 external=openbabel_lane(p,a.openbabel);observable=[r for r in external if r['status']=='observed']
 keys=['declared_state_exact','full_token_record_exact','raw_scale_coordinates_exact','angstrom_positions_bits_exact','element_order_exact','fixed_axes_exact','title_exact']
 ob={'total':len(external),'observed':len(observable),'fixed_axis_inputs':sum(r['source_has_fixed_axes'] for r in external),'fixed_axis_inputs_observed':sum(r['source_has_fixed_axes'] for r in observable)}
 ob.update({k:{'exact':sum(r['comparisons'][k] for r in observable),'total':len(observable)} for k in keys})
 emitted={'valid':[dict(c,text=m.write_viewmol(c['expected'])) for c in p['valid']]};rewrites=openbabel_lane(emitted,a.openbabel)
 sources=['crates/chematic-mol/src/viewmol.rs','crates/chematic-mol/src/viewmol/reader.rs','crates/chematic-mol/src/viewmol/writer.rs','crates/chematic-py/src/formats/viewmol.rs','crates/chematic-wasm/src/viewmol.rs','crates/chematic-cli/src/viewmol.rs','scripts/validation/check_viewmol_preservation.py']
 report={'schema_version':1,'scope':__doc__,'candidate_base':'7704ee2030dda8f9ca54a07eb3fc22e94623fb49','packet_sha256':sha(PACKET),
         'source_sha256':{s:sha(ROOT/s) for s in sources},'runtime_sha256':{'python_extension':sha(m.chematic.__file__),'wasm':sha(a.wasm_module.with_name('chematic_wasm_bg.wasm')),'cli':sha(a.cli)},
         'versions':{'python':sys.version,'node':run(['node','--version']),'openbabel':run([a.openbabel,'-V'])},
         'openbabel_inventory':run([a.openbabel,'-L','vmol']),'unit_probe':unit_probe(a.openbabel),
         'summary':{'candidate':summary,'openbabel':ob},'candidate_receipts':lanes,'openbabel_receipts':external,
         'candidate_text_rewrite':{'scope':'External same-format output observability only; fixed-axis input may be outside the installed OB subset.','observed':sum(r['status']=='observed' for r in rewrites),'total':len(rewrites),'receipts':rewrites}}
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n');print(json.dumps(report['summary'],indent=2))
 good=all(r['exact_all_four_states']==r['valid_total'] and r['projected_bits_all_four']==r['valid_total'] and r['typed_refused']==r['refusal_total'] and r['raw_projection_records_exact']==r['unsafe_projection_total'] and r['unsafe_projection_refused']==r['unsafe_projection_total'] for r in summary.values())
 return 0 if good and report['unit_probe']['exact'] else 1
if __name__=='__main__':raise SystemExit(main())
