from pathlib import Path
import gzip,json,hashlib,subprocess
R=Path('/private/tmp/chematic-openbabel-writer-20261017'); O=Path(__file__).resolve().parent
profiles=[]
def save(kind,name,count,data):
 suffix='sdf' if kind=='sdf' else 'molblocks'; path=O/f'{kind}-{name}.{suffix}'
 path.write_bytes(data); profiles.append({'kind':kind,'profile':name,'records':count,'path':path.name,'sha256':hashlib.sha256(data).hexdigest()})
D=R/'validation/benchmark_corpora'
for name,filename in [('0','openbabel-sdf-speed-1000-v1.sdf.gz'),('5000','openbabel-sdf-speed-slice5000-1000-v1.sdf.gz')]:
 save('sdf',name,1000,gzip.decompress((D/filename).read_bytes()))
for name,filename in [('0','openbabel-v3000-speed-slice0-1000-v1.molblocks.gz'),('5000','openbabel-v3000-speed-slice5000-1000-v1.molblocks.gz'),('chains','openbabel-v3000-speed-chains64-v1.molblocks.gz')]:
 save('v3000',name,1000 if name in ['0','5000'] else 64,gzip.decompress((D/filename).read_bytes()))
chain=(O/'v3000-chains.molblocks').read_text()
continued=[]; in_atoms=False
for line in chain.splitlines():
 if line=='M  V30 BEGIN ATOM': in_atoms=True
 elif line=='M  V30 END ATOM': in_atoms=False
 if in_atoms and line.startswith('M  V30 ') and line.split()[2].isdigit():
  t=line.split(); line='M  V30 '+' '.join(t[2:4])+' -\nM  V30 '+' '.join(t[4:])
 continued.append(line)
save('v3000','continued',64,('\n'.join(continued)+'\n').encode())
for variant in ['spiral','offset']:
 for branched in [False,True]:
  result=[]
  for k in range(1,65):
   n=16*k; total=n+int(branched)
   lines=[f'chain-{k}', '  chematic          3D','', '  0  0  0     0  0            999 V3000','M  V30 BEGIN CTAB',f'M  V30 COUNTS {total} {total-1} 0 0 0','M  V30 BEGIN ATOM']
   for i in range(total):
    if variant=='spiral': x=i*.1251;y=((i*37)%101-50)*.1234;z=((i*73)%97-48)*.2345
    else: x=(i-512)*1.2345;y=((i*19)%53-26)*.3456;z=-((i*43)%79-39)*.4567
    element='F' if i==n else 'C'
    lines.append(f'M  V30 {i+1} {element} {x:.4f} {y:.4f} {z:.4f} 0')
   lines+=['M  V30 END ATOM','M  V30 BEGIN BOND']
   for i in range(1,n): lines.append(f'M  V30 {i} 1 {i} {i+1}')
   if branched: lines.append(f'M  V30 {n} 1 2 {n+1}')
   lines+=['M  V30 END BOND','M  V30 END CTAB','M  END','$$$$']
   result+=lines
  save('v3000_3d',('branched-' if branched else '')+variant,64,('\n'.join(result)+'\n').encode())
plan={'base_commit':json.loads((O/'source-proof.json').read_text())['base_commit'],'candidate':'Decimal pairs and stack-buffer padding in shared fixed4 coordinate formatter','profiles':profiles,'primary':['sdf-0','sdf-5000'],'blocks':21,'repeats':3,'adoption':'Exact source snapshots and rewrites, frozen helper and standard formatter match decimal pairs, random f64 bits and directed boundaries; both primary write paired CI lower>1; all 30 source medians>=0.99; disclose negative intervals','source_valid_epochs':1,'adaptive_reruns':False,'build':'Rust opt3 full LTO 1 CGU 1 job, incremental false; separate empty per-arm targets and separate source roots','comparator':'Open Babel installed native C++ API compiled -O3, RDKit installed observer version recorded','timed_scope':'Native batch parse/prepared write/roundtrip. Inputs loaded and V3000 framed before timer. Diagnostic readers and retained full batches. Startup and plugin loading excluded. Whole-process RSS includes preparation and warmup.','scope':'One macOS host, source build. Synthetic 16..1024-carbon chains, fluorine attached to carbon 2 for branched primary. Arbitrary finite nonplanar coordinates, not physical conformers or holdout. No universal or published-package claim.','coordinates':{'spiral':'x=i*.1251;y=((i*37)%101-50)*.1234;z=((i*73)%97-48)*.2345','offset':'x=(i-512)*1.2345;y=((i*19)%53-26)*.3456;z=-((i*43)%79-39)*.4567'}}
(O/'measurement-plan.json').write_text(json.dumps(plan,indent=2)+'\n')
subprocess.run(['clang++','-O3','-std=c++17','-I/opt/homebrew/include/openbabel3',str(R/'tools/openbabel_v3000_throughput.cpp'),'-L/opt/homebrew/lib','-lopenbabel','-o',str(O/'openbabel-v3000')],check=True)
subprocess.run(['clang++','-O3','-std=c++17','-I/opt/homebrew/include/openbabel3',str(R/'tools/openbabel_sdf_throughput.cpp'),'-L/opt/homebrew/lib','-lopenbabel','-o',str(O/'openbabel-sdf')],check=True)
print('Frozen 10 profiles before timing')
