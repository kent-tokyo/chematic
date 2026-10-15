from pathlib import Path
import sys,json,hashlib,subprocess,importlib.util,platform
from rdkit import Chem
from rdkit.Geometry import Point3D
import rdkit
r=Path('/private/tmp/chematic-openbabel-ring-queue-20261015');s=Path('/private/tmp/chematic-openbabel-ring-queue-source-20261015');sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();plan=json.loads((r/'large-neighborhood-plan.json').read_text());assert plan['declared_before_generation_and_timing'];assert rdkit.__version__=='2026.09.1'
m=json.loads((r/'source-manifest.json').read_text());assert all(sha(s/n)==h for n,h in m['frozen_build_source_sha256'].items())
path=r/'corpus-large-neighborhood.molblocks';original=[];records=[];frames=[]
for depth in plan['recipe']['depths']:
 rw=Chem.RWMol();rw.AddAtom(Chem.Atom(6));rw.AddAtom(Chem.Atom(6));rw.AddBond(0,1,Chem.BondType.DOUBLE);front=[0,1]
 for _ in range(depth):
  nxt=[]
  for parent in front:
   for _ in range(2):
    idx=rw.AddAtom(Chem.Atom(6));rw.AddBond(parent,idx,Chem.BondType.SINGLE);nxt.append(idx)
  front=nxt
 molecule=rw.GetMol();Chem.SanitizeMol(molecule);n=molecule.GetNumAtoms();assert n==2+4*(2**depth-1)
 for variant in plan['recipe']['variants']:
  shift=variant*n//8;order=list(range(shift,n))+list(range(shift))
  if variant%2:order.reverse()
  mol=Chem.RenumberAtoms(molecule,order);name=f'branched-depth-{depth}-variant-{variant}';mol.SetProp('_Name',name);conf=Chem.Conformer(n);conf.Set3D(False)
  for i in range(n):conf.SetAtomPosition(i,Point3D(i/2,(i%7)/4,0))
  mol.AddConformer(conf);text=Chem.MolToMolBlock(mol,forceV3000=True);frames.append(text.rstrip('\n')+'\n$$$$\n');expected=Chem.MolToSmiles(mol,isomericSmiles=True)
  central=mol.GetBonds()[0];assert central.GetBondType()==Chem.BondType.DOUBLE;start=central.GetBeginAtomIdx();goal=central.GetEndAtomIdx();seen={start};front=[start];found=False
  for _ in range(6):
   nxt=[]
   for a in front:
    for bond in mol.GetAtomWithIdx(a).GetBonds():
     if bond.GetIdx()==central.GetIdx():continue
     nb=bond.GetOtherAtomIdx(a)
     if nb==goal:found=True
     if nb not in seen:seen.add(nb);nxt.append(nb)
   front=nxt
  assert not found
  records.append({'name':name,'depth':depth,'variant':variant,'atoms':n,'production_ring_search_visited_atoms':len(seen),'uses_hash_overflow':len(seen)>32,'expected_canonical_smiles':expected});original.append(mol)
assert len(frames)==64 and any(v['uses_hash_overflow'] for v in records) and any(not v['uses_hash_overflow'] for v in records)
path.write_text(''.join(frames));(r/'large-neighborhood-records.json').write_text(json.dumps(records,indent=2)+'\n')
manifest={'id':'v3000-branched-ring-search64-v1','provenance':{'source':'Original project-authored valid acyclic carbon-tree recipe; no third-party molecular rows','generated':'2026-10-15','recipe':'large-neighborhood-plan.json; large_neighborhood.py','transformation':plan['recipe'],'producer':{'rdkit':rdkit.__version__,'python':platform.python_version()}},'license':'MIT OR Apache-2.0','redistribution':{'status':'allowed','rationale':'New synthetic project-authored graphs and recipe under project dual license'},'minimization':'Eight fixed depths cross the 32-visited-atom storage boundary and reach 1022 atoms; eight declared atom permutations check order sensitivity. All 64 records retained. No biological, physical or toolkit-wide representativeness is claimed.','oracle':{'native_source_version':'1.1.4 plus thirteen accepted source candidates versus frozen ring-queue prototype','native_binaries_sha256':{arm:m['binaries'][arm]['sha256'] for arm in ['baseline','candidate']},'observer':'RDKit '+rdkit.__version__,'settings':'SanitizeMol; canonical isomeric MolToSmiles; strict sanitized MolFromMolBlock(removeHs=False); exact coordinates and names; native snapshot/dump repeats=1 expected=64'},'expected':{'records':64,'status':'Every source snapshot/rewrite byte-identical; every source canonical and rewritten graph equals generated RDKit graph; names and exact 2D coordinates preserved','coordinate_units':'angstrom','absolute_tolerance':0,'scope':'Separate post-hoc exploratory source-only diagnostic; never Open Babel ranking or replacement of original matrices'},'sha256':sha(path),'plan_sha256':sha(r/'large-neighborhood-plan.json'),'recipe_sha256':sha(r/'large_neighborhood.py')};(r/'large-neighborhood-corpus-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
a=r/'large-neighborhood.artifacts';a.mkdir(exist_ok=True);streams={}
for arm in ['baseline','candidate']:
 for kind,mode in [('molblocks','dump'),('snapshot','snapshot')]:
  data=subprocess.run([str(r/(arm+'-v3000')),str(path),mode,'1','64'],capture_output=True,check=True,timeout=180).stdout;p=a/(arm+'.'+kind);p.write_bytes(data);streams[arm,kind]=data
for kind in ['molblocks','snapshot']:assert streams['baseline',kind]==streams['candidate',kind],kind
snap=json.loads(streams['candidate','snapshot']);assert len(snap)==64
rewrites=streams['candidate','molblocks'].decode().split('$$$$\n');assert rewrites[-1]=='' and len(rewrites)==65
for expected,mol,row,block in zip(records,original,snap,rewrites[:-1],strict=True):
 parsed=Chem.MolFromSmiles(row['canonical_smiles']);assert parsed is not None and Chem.MolToSmiles(parsed,isomericSmiles=True)==expected['expected_canonical_smiles'],expected['name']
 parsed=Chem.MolFromMolBlock(block,sanitize=True,removeHs=False,strictParsing=True);assert parsed is not None and Chem.MolToSmiles(parsed,isomericSmiles=True)==expected['expected_canonical_smiles'],expected['name'];assert parsed.GetProp('_Name')==expected['name']
 for i in range(mol.GetNumAtoms()):
  p=mol.GetConformer().GetAtomPosition(i);q=parsed.GetConformer().GetAtomPosition(i);assert (p.x,p.y,p.z)==(q.x,q.y,q.z)
gate={'all_64_source_snapshots_and_rewrites_byte_identical':True,'all_64_original_graphs_names_and_exact_coordinates_preserved':True,'sha256':{p.name:sha(p) for p in a.iterdir()},'corpus_manifest_sha256':sha(r/'large-neighborhood-corpus-manifest.json'),'post_hoc_exploratory_not_openbabel_comparison':True};(r/'large-neighborhood-equivalence.json').write_text(json.dumps(gate,indent=2)+'\n');print('All 64 generated records sanitized and passed exact source/RDKit gates; timing starts now',flush=True)
spec=importlib.util.spec_from_file_location('diagnostic',r/'source_diagnostics.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
try:module.run_profile('large-neighborhood21',path,64,gate)
finally:
 if module.reserve.exists():module.reserve.unlink()
