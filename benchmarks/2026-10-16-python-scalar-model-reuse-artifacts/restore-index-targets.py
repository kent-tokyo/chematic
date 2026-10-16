from pathlib import Path
from urllib.parse import unquote
import json, re, subprocess
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
OUT=Path(__file__).resolve().parent
targets={'benchmarks/README.md','benchmarks/archive-index.md'}
for path in sorted(targets):
 text=subprocess.check_output(['git','show','HEAD:'+path],cwd=ROOT,text=True)
 for raw in re.findall(r'\[[^\]]+\]\(([^)]+)\)',text):
  target=raw.strip().split(None,1)[0].split('#',1)[0]
  if target and not target.startswith(('https:','http:','mailto:')):
   p=(ROOT/'benchmarks'/unquote(target)).resolve()
   assert p.is_relative_to(ROOT),target
   if not p.is_file(): targets.add(p.relative_to(ROOT).as_posix())
tree=set(subprocess.check_output(['git','ls-tree','-r','--name-only','HEAD'],cwd=ROOT,text=True).splitlines())
assert targets<=tree,targets-tree
paths=sorted(targets)
data=subprocess.check_output(['git','cat-file','--batch'],input=''.join('HEAD:'+p+'\n' for p in paths).encode(),cwd=ROOT)
at=0;written=0
for path in paths:
 end=data.index(b'\n',at);header=data[at:end].decode().split();size=int(header[2]);at=end+1
 raw=data[at:at+size];at+=size+1
 target=ROOT/path;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw);written+=size
assert at==len(data)
(OUT/'restored-index-targets.json').write_text(json.dumps({'tracked_paths':paths,'bytes':written,'reason':'Restore only index targets in sparse checkout; all unchanged Git content'},indent=2)+'\n')
print('Restored',len(paths),'unchanged index targets, bytes',written)
