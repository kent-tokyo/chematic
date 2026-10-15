from pathlib import Path
import subprocess,json,hashlib
r=Path('/private/tmp/chematic-huckel-fifo-20261015');s=Path('/private/tmp/chematic-huckel-fifo-source-20261015');py='/private/tmp/chematic-cycle-speed-20261014/env/bin/python';m=json.loads((r/'source-manifest.json').read_text())
def verify():
 assert all(hashlib.sha256((s/n).read_bytes()).hexdigest()==h for n,h in m['candidate_sha256'].items())
for offset in [5000,0]:
 verify()
 cmd=[py,str(s/'scripts/bench_python_speed_frontier.py'),'--corpus',str(s/'validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi'),'--offset',str(offset),'--limit','1000','--blocks','20','--baseline-wheel',str(next((r/'baseline-wheel').glob('*.whl'))),'--baseline-path',str(r/'baseline'),'--candidate-wheel',str(next((r/'candidate-wheel').glob('*.whl'))),'--candidate-path',str(r/'candidate'),'--source-manifest',str(r/'source-manifest.json'),'--base-revision','cac480fb23cd8adeec6e5e9bed1f18677f87b392-plus-ten-accepted-source-candidates','--output',str(r/f'frontier-{offset}.json')]
 with (r/f'frontier-{offset}.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
 verify();print(offset,'complete',flush=True)
