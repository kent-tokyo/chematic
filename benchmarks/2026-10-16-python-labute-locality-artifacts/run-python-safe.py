from pathlib import Path
import errno, sys, time
OUT=Path(__file__).resolve().parent
ROOT=Path('/private/tmp/chematic-labute-locality-20261016')
reserve=OUT/'enospc-checkpoint-reserve.bin'
try:
    with reserve.open('wb') as output:
        for _ in range(64):output.write(b'x'*1024*1024)
except BaseException:
    reserve.unlink(missing_ok=True)
    raise
original=Path.write_text

def guarded(path,*args,**kwargs):
    for attempt in range(60):
        try:return original(path,*args,**kwargs)
        except OSError as error:
            if error.errno!=errno.ENOSPC or not path.is_relative_to(OUT):raise
            if reserve.exists():
                reserve.unlink()
                print('Released owned 64MiB reserve; retrying the same in-memory checkpoint only',flush=True)
            if attempt==59:raise
            time.sleep(1)
Path.write_text=guarded
sys.path.insert(0,str(ROOT/'scripts'))
import bench_python_speed_frontier as bench
try:bench.main()
finally:reserve.unlink(missing_ok=True)
