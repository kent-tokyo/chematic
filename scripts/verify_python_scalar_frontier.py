#!/usr/bin/env python3
"""Retain every selected descriptor/fingerprint value or refusal on fixed corpora."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'validation/cosmolkit_comparison'))
from bench_corpus import table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--engine', choices=('chematic','rdkit','cosmolkit'), required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--ops', default='tpsa,labute_asa,num_rings')
    args = parser.parse_args()
    ops = args.ops.split(',')
    if set(ops) - {'tpsa','labute_asa','num_rings','rotatable_bonds','mqn','fp_layered','qed'}:
        parser.error('unsupported operation')
    if args.engine == 'rdkit':
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')
    parse, functions = table(args.engine)
    corpora = [ROOT / 'validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi',
               ROOT / 'validation/cosmolkit_comparison/corpora/stress-v1.smi']
    with args.output.open('w') as out:
        for corpus in corpora:
            data = corpus.read_bytes()
            rows = [line.split()[0] for line in data.decode().splitlines() if line.strip()]
            out.write(json.dumps({'corpus':corpus.name,'sha256':hashlib.sha256(data).hexdigest(),
                                  'rows':len(rows),'engine':args.engine})+'\n')
            for index, smiles in enumerate(rows):
                row = {'corpus':corpus.name,'index':index,'values':{}}
                try:
                    molecule = parse(smiles)
                    if molecule is None:raise ValueError('parse failed')
                    for op in ops:
                        try:
                            value = functions[op](molecule)
                            if op == 'fp_layered':
                                from bench_python_speed_frontier import normal
                                value = normal(value, args.engine, op)
                            if op in ('mqn','fp_layered'):
                                if (not isinstance(value,list) or any(type(v) is not int for v in value)
                                        or (op == 'mqn' and len(value) != 42)):
                                    raise ValueError('noninteger vector')
                                row['values'][op] = {'value':value}
                            else:
                                if type(value) not in (int,float):raise ValueError('nonnumeric scalar')
                                row['values'][op] = ({'value':value} if math.isfinite(value)
                                                     else {'unavailable':'nonfinite'})
                        except Exception as error:
                            row['values'][op] = {'error':type(error).__name__+': '+str(error)}
                except Exception as error:
                    row['parse_error'] = type(error).__name__+': '+str(error)
                out.write(json.dumps(row,sort_keys=True,allow_nan=False)+'\n')
            print(corpus.name,len(rows),'done',flush=True)


if __name__ == '__main__':main()
