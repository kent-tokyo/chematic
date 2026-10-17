from pathlib import Path
import sys,time,json
import chematic
R=Path('/private/tmp/chematic-python-asa-preparation-20261017')
offset=int(sys.argv[1]); rows=[line.split()[0] for line in (R/'validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi').read_text().splitlines() if line.strip()][offset:offset+1000]
assert len(rows)==1000
started=time.monotonic();loops=0;checksum=0.0
while time.monotonic()-started<8:
 molecules=[chematic.from_smiles(s) for s in rows]
 checksum+=sum(m.labute_asa for m in molecules);loops+=1
print(json.dumps({'scope':'Diagnostic only, mixed parse/getter/free process; attribution must use getter ancestry, not global counts. No formal timing samples.','offset':offset,'loops':loops,'checksum':checksum}))
