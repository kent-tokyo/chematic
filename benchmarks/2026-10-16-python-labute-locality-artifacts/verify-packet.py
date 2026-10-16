from pathlib import Path
import gzip,hashlib,json,re,subprocess,sys,tempfile
ROOT=Path('/private/tmp/chematic-labute-locality-20261016')
OUT=Path(__file__).resolve().parent
PACKET=ROOT/'benchmarks/2026-10-16-python-labute-locality-artifacts'
sys.path.insert(0,str(ROOT/'scripts'));import bench_python_speed_frontier as bench
manifest=json.loads((OUT/'source-manifest.json').read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert all(sha(ROOT/p)==v for p,v in manifest['candidate_sha256'].items())
changed=manifest['changed']+manifest['added']
with tempfile.TemporaryDirectory(prefix='chematic-labute-patch-proof-',dir='/private/tmp') as d:
    root=Path(d)
    subprocess.run(['git','init','-q',str(root)],check=True)
    for path in manifest['changed']:
        original=subprocess.check_output(['git','show',manifest['base']+':'+path],cwd=ROOT)
        target=root/path;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(original)
        assert sha(target)==manifest['baseline_sha256'][path]
    subprocess.run(['git','apply',str(OUT/'candidate.patch')],cwd=root,check=True)
    assert all(sha(root/path)==manifest['candidate_sha256'][path] for path in changed)
checks=0;strict=0;unranked=[]
for offset in [0,5000]:
    data=json.loads((OUT/f'frontier-{offset}.json').read_text())
    assert data['complete'] and len(data['blocks'])==20
    assert data['base_revision']==manifest['base']
    assert data['source_manifest']==manifest
    assert data['corpus']['sha256']==sha(ROOT/'validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi')
    for op in bench.OPS:
        for mode in bench.MODES:
            key=op+'/'+mode
            for seed,arm in enumerate(['baseline','rdkit','cosmolkit']):
                value=data['summary'][key][arm]
                allowed=bench.comparable(data['agreement'][op],arm)
                expected=bench.paired_median_ratio_ci([b[arm]['operations'][op][mode] for b in data['blocks']],[b['candidate']['operations'][op][mode] for b in data['blocks']],seed=seed) if allowed else None
                assert value['speedup']==expected,(offset,key,arm)
                win=bool(expected and expected['ci95'][0]>1 and expected['winning_blocks']==20)
                assert value['demonstrated_win']==win
                checks+=1;strict+=int(win)
                if not allowed:unranked.append([offset,key,arm])
    for arm in ['baseline','candidate','rdkit','cosmolkit']:
        raw=gzip.decompress((OUT/f'observer-{offset}-{arm}.json.gz').read_bytes())
        observer=json.loads(raw)
        assert all(b[arm]['output_sha256']==observer['output_sha256'] for b in data['blocks'])
parity=json.loads((OUT/'parity.json').read_text());assert parity['outputs']['baseline']==parity['outputs']['candidate']
for arm,values in parity['outputs'].items():
    for kind,value in values.items():assert hashlib.sha256(gzip.decompress((OUT/f'{arm}-{kind}.jsonl.gz').read_bytes())).hexdigest()==value['sha256']
restart=json.loads((OUT/'administrative-restart-plan.json').read_text());interrupted=json.loads((OUT/'interrupted-frontier-5000.json').read_text());assert not interrupted['complete'] and len(interrupted['blocks'])==restart['saved_complete_blocks']==7
proof={'source_files_verified':len(manifest['candidate_sha256']),'patch_reconstructed_from_base':manifest['base'],'changed_and_added_files_verified':changed,'all_matrix_comparisons_recomputed':checks,'strict_wins_including_source':strict,'all_unranked_lanes':unranked,'all_20_hashes_match_fresh_observers':True,'full_source_dump_hashes_verified':True,'interrupted_epoch_retained_separately':True,'universal_superiority':False}
(OUT/'packet-audit.json').write_text(json.dumps(proof,indent=2)+'\n')
(PACKET/'packet-audit.json').write_bytes((OUT/'packet-audit.json').read_bytes())
files={p.name:sha(p) for p in sorted(PACKET.iterdir()) if p.is_file() and p.name!='packet.sha256.json'}
(PACKET/'packet.sha256.json').write_text(json.dumps(files,indent=2)+'\n')
assert all(sha(PACKET/p)==value for p,value in files.items())
print(json.dumps(proof,indent=2))
