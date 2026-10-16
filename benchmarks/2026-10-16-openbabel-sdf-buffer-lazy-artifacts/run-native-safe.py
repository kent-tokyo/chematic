"""Run the unchanged native runner; only retry owned checkpoint writes."""
from pathlib import Path
import errno
import runpy
import sys

ROOT = Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
EVIDENCE = Path(__file__).resolve().parent
reserve = EVIDENCE / 'checkpoint-reserve.tmp'
with reserve.open('wb') as stream:
    stream.write(bytes(64 * 1024 * 1024))
original = Path.write_text

def safe_write(path, data, *args, **kwargs):
    try:
        return original(path, data, *args, **kwargs)
    except OSError as error:
        if error.errno != errno.ENOSPC or path.parent != EVIDENCE or path.suffix != '.partial':
            raise
        reserve.unlink(missing_ok=True)
        return original(path, data, *args, **kwargs)

Path.write_text = safe_write
runner = ROOT / 'scripts' / sys.argv[1]
sys.argv = [str(runner)] + sys.argv[2:]
sys.path.insert(0, str(runner.parent))
try:
    runpy.run_path(str(runner), run_name='__main__')
finally:
    reserve.unlink(missing_ok=True)
