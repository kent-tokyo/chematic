import os,sys,json,time,hashlib
from pathlib import Path
import chematic
source=Path('/private/tmp/chematic-basis-roots-source-20261015')
offset=int(sys.argv[1]);op=sys.argv[2];path=source/'validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi';data=path.read_bytes();rows=[l.split()[0] for l in data.decode().splitlines() if l.strip()][offset:offset+1000];assert len(rows)==1000
print(json.dumps({'ready':True,'pid':os.getpid(),'offset':offset,'op':op,'rows':len(rows),'version':chematic.__version__,'corpus_sha256':hashlib.sha256(data).hexdigest()}),flush=True)
stop=time.monotonic()+8;count=0;value=0.0
while time.monotonic()<stop:
 molecules=[chematic.from_smiles(s) for s in rows]
 for m in molecules:value+=getattr(m,op)
 count+=len(molecules)
 del molecules
print(json.dumps({'complete':True,'evaluations':count,'checksum':value.hex(),'errors':0}),flush=True)
