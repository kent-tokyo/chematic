#!/usr/bin/env python3
"""Conservative, time-bounded invalid graph overlap audit by scaffold."""
import argparse
import hashlib
import json
import multiprocessing as mp
from collections import Counter
from pathlib import Path
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks'))
from platinum_invalid_exclusion_audit import relaxed_graph_query


def check(query, targets, connection):
    try:
        connection.send([i for i, target in enumerate(targets) if target.HasSubstructMatch(query)])
    finally:
        connection.close()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('sdf', type=Path)
    p.add_argument('split', type=Path)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--seconds', type=float, default=0.5)
    args = p.parse_args()
    RDLogger.DisableLog('rdApp.*')
    split_bytes = args.split.read_bytes(); data = json.loads(split_bytes)
    source_hash = hashlib.sha256(args.sdf.read_bytes()).hexdigest()
    if source_hash != data['source_sha256']:
        raise ValueError('source mismatch')
    targets=[]
    for entry in data['invalid_exclusions']:
        m=Chem.MolFromSmiles(entry['smiles'], sanitize=False)
        if m is None:
            raise ValueError('unrepresentable exclusion')
        m.UpdatePropertyCache(strict=False); Chem.GetSymmSSSR(m); targets.append(m)
    supplier = Chem.SDMolSupplier(str(args.sdf), removeHs=False)
    if len(supplier)!=len(data['rows']):
        raise ValueError('record count mismatch')
    cache={}; overlaps=[]; indeterminate=[]
    ctx=mp.get_context('fork')
    for row,mol in zip(data['rows'], supplier):
        if row.get('split')!='holdout':continue
        key=row['scaffold_sha256']
        if key not in cache:
            try:
                heavy=Chem.RemoveHs(mol); graph=Chem.MolFromSmiles(Chem.MolToSmiles(heavy,isomericSmiles=True))
                scaf=MurckoScaffold.GetScaffoldForMol(graph)
                if scaf.GetNumAtoms()==0:
                    cache[key]=list(range(len(targets)))
                else:
                    query=relaxed_graph_query(scaf); parent,child=ctx.Pipe(duplex=False)
                    process=ctx.Process(target=check,args=(query,targets,child));process.start();child.close()
                    if parent.poll(args.seconds):
                        try:cache[key]=parent.recv()
                        except EOFError:cache[key]=None
                    else:cache[key]=None
                    if process.is_alive():process.terminate()
                    process.join(timeout=1)
                    if process.is_alive():process.kill();process.join()
                    parent.close()
            except Exception:
                cache[key]=None
        if cache[key] is None:
            indeterminate.append({'index':row['index'],'scaffold_sha256':key})
        elif cache[key]:
            overlaps.append({'index':row['index'],'split':'holdout','invalid_indices':cache[key], 'scaffold_sha256':key})
        if len(cache)%200==0:print('scaffolds',len(cache),flush=True)
    out={'source_sha256':source_hash, 'split_sha256':hashlib.sha256(split_bytes).hexdigest(),
         'rdkit_version':rdBase.rdkitVersion,'invalid_exclusion_count':len(targets),
         'overlap_counts':dict(Counter(x['split'] for x in overlaps)), 'overlaps':overlaps,
         'indeterminate':indeterminate,'eligible_scaffolds':len(cache),'seconds_per_scaffold':args.seconds,
         'rule':'exclude overlap and every indeterminate scaffold; atom-number substructure ignoring bond order',
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print('overlap',len(overlaps),'indeterminate',len(indeterminate),flush=True)


if __name__=='__main__':main()
