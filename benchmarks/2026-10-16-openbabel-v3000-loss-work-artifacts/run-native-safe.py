"""Unchanged native runner; retry only identical in-memory output on ENOSPC."""
from pathlib import Path
import errno, runpy, sys, time
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
EVIDENCE=Path(__file__).resolve().parent
reserve=EVIDENCE/'checkpoint-reserve.tmp'
try:
    with reserve.open('wb') as output:
        for _ in range(64):output.write(b'x'*1024*1024)
except BaseException:
    reserve.unlink(missing_ok=True)
    raise
original_text=Path.write_text
original_bytes=Path.write_bytes

def guarded(method,path,*args,**kwargs):
    for attempt in range(60):
        try:return method(path,*args,**kwargs)
        except OSError as error:
            if error.errno!=errno.ENOSPC or not path.is_relative_to(EVIDENCE):raise
            if reserve.exists():
                reserve.unlink()
                print('Released owned 64MiB reserve; same in-memory write only',flush=True)
            if attempt==59:raise
            time.sleep(1)
Path.write_text=lambda path,*args,**kwargs:guarded(original_text,path,*args,**kwargs)
Path.write_bytes=lambda path,*args,**kwargs:guarded(original_bytes,path,*args,**kwargs)
runner=ROOT/'scripts'/sys.argv[1]
sys.argv=[str(runner)]+sys.argv[2:]
sys.path.insert(0,str(runner.parent))
try:runpy.run_path(str(runner),run_name='__main__')
finally:reserve.unlink(missing_ok=True)
