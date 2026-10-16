from pathlib import Path
import gzip,hashlib,json,subprocess,sys,tempfile
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
import bench_openbabel_sdf_throughput as sdf
import bench_openbabel_v3000_throughput as v3
from rdkit import Chem
manifest=json.loads((OUT/'source-manifest.json').read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert all(sha(ROOT/p)==value for p,value in manifest['source_files_sha256'].items())
assert all(sha(Path(p))==value for p,value in manifest['native_libraries_sha256'].items())
seal=json.loads((OUT/'binary-seal.json').read_text());assert all(sha(OUT/p)==value for p,value in seal['files_sha256'].items())
with tempfile.TemporaryDirectory(prefix='chematic-v3000-loss-replay-',dir='/private/tmp') as d:
    root=Path(d);subprocess.run(['git','init','-q',str(root)],check=True)
    for path in manifest['changed']:
        target=root/path;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(subprocess.check_output(['git','show',manifest['base_revision']+':'+path],cwd=ROOT))
    subprocess.run(['git','apply',str(OUT/'candidate.patch')],cwd=root,check=True)
    assert all(sha(root/path)==manifest['source_files_sha256'][path] for path in manifest['changed']+manifest['added'])
profiles=json.loads((OUT/'measurement-plan.json').read_text())['profiles'];proof={};comparisons=0
for item in profiles:
    kind=item['kind'];profile=item['profile'];key=kind+'-'+profile
    data=json.loads((OUT/f'paired-{key}.json').read_text());assert data['complete'] and len(data['blocks'])==21 and data['source_manifest']==manifest
    suffix='sdf' if kind=='sdf' else 'molblocks';corpus=OUT/f'{key}.{suffix}';assert sha(corpus)==data['corpus']['sha256']
    reference=sdf.observe(corpus.read_bytes());artifacts=OUT/f'paired-{key}.artifacts'
    outputs={}
    for arm in ['baseline','candidate','openbabel']:
        raw=(artifacts/f'{arm}.{suffix}').read_bytes();outputs[arm]=raw;assert hashlib.sha256(raw).hexdigest()==data['rewrite_sha256'][arm]
        observed=sdf.semantic_check(reference,sdf.observe(raw)) if kind=='sdf' else v3.rewrite_check(reference,raw)
        assert observed==data['semantics'][arm]
        if arm!='openbabel':
            snap=gzip.decompress((OUT/f'preflight-{key}-{arm}-snapshot.gz').read_bytes())
            assert hashlib.sha256(snap).hexdigest()==data['snapshot_sha256'][arm]
            strings=[row['canonical_smiles'] for row in json.loads(snap)]
        elif kind=='v3000':strings=[line.split()[0] for line in (artifacts/'openbabel.smiles').read_text().splitlines()]
        else:
            raw_smiles=subprocess.check_output([str(OUT/'openbabel-sdf'),str(corpus),'smiles','1',str(item['records'])])
            (OUT/f'post-timing-{key}-openbabel-smiles.gz').write_bytes(gzip.compress(raw_smiles,mtime=0))
            strings=[line.split()[0] for line in raw_smiles.decode().splitlines()]
        assert len(strings)==item['records']
        rejected=[];bad=[]
        for i,(row,string) in enumerate(zip(reference,strings,strict=True)):
            mol=Chem.MolFromSmiles(string)
            if mol is None:rejected.append(i)
            elif Chem.MolToSmiles(mol)!=row['graph']:bad.append(i)
        expected=data['semantics'][arm+'_parse'];assert expected['passed']==(not rejected and not bad) and expected['mismatches']==len(bad)
        if kind=='v3000':assert expected['mismatch_rows']==bad and expected['rejected_rows']==rejected
        else:assert not rejected and expected['examples']==bad[:8]
    assert outputs['baseline']==outputs['candidate']
    for mode,values in data['summary'].items():
        for seed,arm in enumerate(['baseline','openbabel']):
            allowed=sdf.eligible(data['semantics'],data['source_parse_snapshots_equal'],mode,arm)
            interval=sdf.paired_median_ratio_ci([b[arm][mode]['elapsed_ns'] for b in data['blocks']],[b['candidate'][mode]['elapsed_ns'] for b in data['blocks']],seed=seed) if allowed else None
            assert values[arm]['semantic_prerequisite']==allowed and values[arm]['speedup']==interval
            assert values[arm]['strict_speed_win']==bool(interval and interval['ci95'][0]>1 and interval['winning_blocks']==21)
            comparisons+=1
    proof[key]={'complete_blocks':21,'source_snapshots_and_rewrites_exact':True,'saved_rewrites_and_parse_semantics_revalidated':True}
result={'source_file_hashes_verified':len(manifest['source_files_sha256']),'native_library_entries_verified':len(manifest['native_libraries_sha256']),'patch_reconstructed_from':manifest['base_revision'],'all_comparisons_recomputed':comparisons,'profiles':proof,'post_timing_sdf_smiles_observers_not_timing_samples':True,'universal_superiority':False}
(OUT/'packet-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
