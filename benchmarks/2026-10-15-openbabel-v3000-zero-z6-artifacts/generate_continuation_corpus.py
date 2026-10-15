"""Reframe every atom/bond row of the unchanged synthetic chain profile."""
from pathlib import Path
import argparse,hashlib
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--input',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();lines=[]
for raw in args.input.read_text(encoding='utf-8').splitlines(keepends=True):
    payload=raw.removesuffix('\n').removeprefix('M  V30 ')
    if raw.startswith('M  V30 ') and payload.split()[0].isdigit():
        first,rest=payload.split(' ',1)
        lines.extend([f'M  V30 {first} -\n',f'M  V30 {rest}\n'])
    else:
        lines.append(raw)
args.output.write_text(''.join(lines),encoding='utf-8',newline='\n');print(hashlib.sha256(args.output.read_bytes()).hexdigest())
