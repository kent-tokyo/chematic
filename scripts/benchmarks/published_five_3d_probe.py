#!/usr/bin/env python3
"""Compare five published wheels on identical exposed molecules, one worker at a time."""
import argparse
import hashlib
import json
import multiprocessing as mp
import platform
import resource
import sys
import time
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def worker(conn, engine, package_dir, smiles):
    try:
        if package_dir:
            sys.path.insert(0, package_dir)
        from platinum_development_probe import graph_matches
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        if engine.startswith('chematic'):
            import chematic
            if chematic.__version__ != engine.removeprefix('chematic-'):
                raise ValueError('wrong chematic version loaded')
            loaded = {'version': chematic.__version__, 'extension_sha256': digest(chematic.chematic.__file__)}
            other = chematic.from_smiles(smiles).add_hydrogens()
            if not graph_matches(other, mol):
                raise ValueError('indexed graph mismatch')
        elif engine == 'cosmolkit':
            import cosmolkit as ck
            if ck.__version__ != '0.6.3': raise ValueError('wrong COS version loaded')
            loaded = {'version': ck.__version__, 'package_file': ck.__file__}
            other = ck.Molecule.from_smiles(smiles).with_hydrogens()
            atoms = [a.atomic_number for a in other.atoms]
            bonds = sorted((min(b.begin,b.end),max(b.begin,b.end),b.order_name.upper()) for b in other.bonds)
            rbonds = sorted((min(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),max(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),str(b.GetBondType()).upper()) for b in mol.GetBonds())
            if atoms != [a.GetAtomicNum() for a in mol.GetAtoms()] or bonds != rbonds:
                raise ValueError('indexed graph mismatch')
        if engine == 'rdkit': loaded = {'version': rdBase.rdkitVersion, 'extension_sha256': digest(rdBase.__file__)}
        conn.send(('ready', loaded))
        if conn.recv() != 'start':
            return
        for attempt in range(20):
            seed = 42 * (attempt + 1)
            candidate = None
            if engine.startswith('chematic'):
                result = other.rdkit_embed_multiple_confs_mmff(1, random_seed=seed, max_mmff_iters=200)
                if result:
                    candidate = (result[0][0], result[0][2])
            elif engine == 'rdkit':
                trial = Chem.Mol(mol)
                params = AllChem.ETKDGv3(); params.randomSeed = seed; params.numThreads = 1
                if AllChem.EmbedMolecule(trial, params) == 0:
                    status = AllChem.MMFFOptimizeMolecule(trial, maxIters=200)
                    candidate = (status, [list(trial.GetConformer().GetAtomPosition(i)) for i in range(trial.GetNumAtoms())])
            else:
                other.embed_3d_conformers_(1, random_seed=seed, num_threads=1)
                result = other.with_mmff_optimized_conformers(num_threads=1, max_iterations=200)
                if result.molecule.conformers_3d:
                    candidate = (result.conformer_results[0].status_code, result.molecule.conformers_3d[0].coordinates)
            conn.send(('candidate', (candidate, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)))
        conn.send(('done', None))
    except BaseException as exc:
        try: conn.send(('error', f'{type(exc).__name__}: {exc}'))
        except (EOFError, OSError): pass
    finally:
        conn.close()


def run(engine, directory, smiles, seconds):
    parent, child = mp.get_context('spawn').Pipe()
    process = mp.get_context('spawn').Process(target=worker,args=(child,engine,directory,smiles))
    process.start(); child.close()
    out = {'status':'error','error':None,'attempted':0,'candidates':[], 'arrival_seconds':[], 'worker_peak_rss':None}
    start = None
    try:
        if not parent.poll(30):
            out['error']='startup timeout'; return out
        kind,value=parent.recv()
        if kind != 'ready':
            out['error']=value; return out
        out['loaded'] = value
        start=time.monotonic(); parent.send('start'); deadline=start+seconds
        while True:
            remaining=deadline-time.monotonic()
            if remaining <= 0 or not parent.poll(remaining):
                out['status']='timeout'; break
            kind,value=parent.recv(); arrival=time.monotonic()-start
            if arrival > seconds:
                out['status']='timeout'; break
            if kind == 'candidate':
                candidate,out['worker_peak_rss']=value; out['attempted']+=1
                if candidate is not None:
                    out['candidates'].append(candidate); out['arrival_seconds'].append(arrival)
            elif kind == 'done':
                out['status']='ok'; break
            else:
                out['error']=value; break
    except (EOFError, OSError) as exc:
        out['error']=repr(exc)
    finally:
        if process.is_alive(): process.terminate()
        process.join(timeout=1)
        if process.is_alive(): process.kill(); process.join()
        parent.close()
        if start is not None: out['elapsed_with_cleanup']=time.monotonic()-start
    return out


def main():
    from cod_robust_scoring import annotate
    from platinum_geometry import prepare_geometry_checks
    p=argparse.ArgumentParser(); p.add_argument('--sdf',required=True); p.add_argument('--split',required=True); p.add_argument('--output',required=True); p.add_argument('--limit',type=int,default=100); p.add_argument('--seconds',type=float,default=10.)
    args=p.parse_args(); split=json.loads(Path(args.split).read_text())
    if digest(args.sdf)!=split['source_sha256']: raise ValueError('SDF hash mismatch')
    rows=sorted((r for r in split['rows'] if r.get('split')=='holdout'), key=lambda r:hashlib.sha256(r['name'].encode()).digest())[:args.limit]
    engines={'chematic-1.2.2':'/private/tmp/chematic-goal-122','chematic-1.2.3':'/private/tmp/chematic-goal-123','chematic-1.2.4':'/private/tmp/chematic-goal-124-published','rdkit':None,'cosmolkit':'/private/tmp/chematic-bench-cos063'}
    wheels={p.name:digest(p) for p in sorted(Path('/private/tmp/chematic-goal-wheels').glob('*.whl'))}
    metadata={'type':'metadata','sdf_sha256':digest(args.sdf),'split_sha256':digest(args.split),'script_sha256':digest(__file__),'scorer_hashes':{n:digest(Path(__file__).with_name(n)) for n in ('cod_robust_scoring.py','platinum_development_probe.py','platinum_stereo_strict.py','platinum_geometry.py')},'wheels':wheels,'rdkit':rdBase.rdkitVersion,'platform':platform.platform(),'count':20,'seconds':args.seconds,'timed_prefix_seconds':0.2,'total':len(rows),'selection':'first 100 SHA-256 name-sorted v3 holdout records; previously exposed','seed_schedule':'42 * one-based attempt','timing_scope':'generation and minimization after per-engine preparation and worker startup; cleanup included in elapsed_with_cleanup','rss_scope':'worker high-water mark at last completed attempt, bytes on macOS; import included; timeouts may be undercounted','fixed_acceptance':'complete 20 attempts before deadline','engines':list(engines)}
    supplier=Chem.SDMolSupplier(args.sdf,removeHs=False)
    with open(args.output,'x') as f:
        f.write(json.dumps(metadata)+'\n'); f.flush()
        for num,row in enumerate(rows):
            result={'type':'result','index':row['index'],'name':row['name'],'engines':{}}
            try:
                ref=supplier[row['index']]
                if ref is None or ref.GetProp('_Name')!=row['name']: raise ValueError('input identity mismatch')
                ref=Chem.RemoveHs(ref); smiles=Chem.MolToSmiles(ref,isomericSmiles=True); mol=Chem.AddHs(Chem.MolFromSmiles(smiles)); checks=prepare_geometry_checks(mol)
                result['geometry_reference'] = {'bonds_supported':len(checks[0]), 'bonds_total':mol.GetNumBonds(), 'angles_supported':len(checks[1]), 'angles_total':sum(a.GetDegree() * (a.GetDegree() - 1) // 2 for a in mol.GetAtoms())}
                names=list(engines); order=names[num%len(names):]+names[:num%len(names)]; result['order']=order
                for name in order:
                    out=run(name,engines[name],smiles,args.seconds)
                    try: out=annotate(out,mol,ref,checks)
                    except Exception as exc: out['scoring_error']=repr(exc)
                    result['engines'][name]=out
            except Exception as exc: result['input_error']=repr(exc)
            f.write(json.dumps(result,allow_nan=False)+'\n'); f.flush(); print(num+1,'/',len(rows),flush=True)
        f.write(json.dumps({'type':'complete','rows':len(rows)})+'\n')

if __name__ == '__main__': main()
