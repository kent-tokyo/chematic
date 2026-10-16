from pathlib import Path
import sys,errno,time
r=Path('/private/tmp/chematic-aromatic-proof-evidence-20261016');s=Path('/private/tmp/chematic-aromatic-proof-20261016');reserve=r/'python-checkpoint-reserve.bin'
with reserve.open('wb') as f:
 for _ in range(64):f.write(b'x'*1024*1024)
original_text=Path.write_text;original_bytes=Path.write_bytes

def guarded(method,path,*args,**kwargs):
 for attempt in range(60):
  try:return method(path,*args,**kwargs)
  except OSError as error:
   if error.errno!=errno.ENOSPC or not path.is_relative_to(r):raise
   if reserve.exists():reserve.unlink();print('Released owned 64MiB reserve; current observation remains in memory',flush=True)
   if attempt==59:raise
   time.sleep(1)
Path.write_text=lambda path,*args,**kwargs:guarded(original_text,path,*args,**kwargs)
Path.write_bytes=lambda path,*args,**kwargs:guarded(original_bytes,path,*args,**kwargs)
sys.path.insert(0,str(s/'scripts'));import bench_python_speed_frontier as bench
try:bench.main()
finally:
 if reserve.exists():reserve.unlink()
