"""Check all five packets, source patches and frozen adoption decisions."""
from pathlib import Path
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
NAMES = [
    '2026-10-16-openbabel-sdf-buffered',
    '2026-10-16-openbabel-sdf-buffered-main',
    '2026-10-16-openbabel-sdf-buffer-lazy',
    '2026-10-16-openbabel-sdf-direct-fields',
    '2026-10-16-openbabel-sdf-loss-work',
]

def sha(data):
    return hashlib.sha256(data).hexdigest()

def main():
    audited = []
    for index, name in enumerate(NAMES):
        folder = ROOT/'benchmarks'/(name+'-artifacts')
        manifest = json.loads((folder/'source-manifest.json').read_text())
        decision = json.loads((folder/'decision.json').read_text())
        reconstruction = json.loads((folder/'reconstruction.json').read_text())
        patch = gzip.decompress((folder/'candidate.patch.gz').read_bytes())
        assert sha(patch) == reconstruction['patch_sha256']
        with tempfile.TemporaryDirectory(prefix='chematic-sdf-packet-') as temp:
            root = Path(temp)
            for n in reconstruction['changed_source_files']:
                old = subprocess.run(['git','show',manifest['base_revision']+':'+n],cwd=ROOT,capture_output=True)
                if old.returncode == 0:
                    p=root/n
                    p.parent.mkdir(parents=True,exist_ok=True)
                    p.write_bytes(old.stdout)
            subprocess.run(['git','apply','--unsafe-paths','-'],cwd=root,input=patch,check=True)
            for n in reconstruction['changed_source_files']:
                assert sha((root/n).read_bytes()) == manifest['source_files_sha256'][n], n
        negative = {}
        below_floor = {}
        primary_pass = True
        for profile in ['sdf-0','sdf-5000','v3000-0','v3000-5000','v3000-chains']:
            d=json.loads(gzip.decompress((folder/f'paired-{profile}.json.gz').read_bytes()))
            assert d['complete'] and len(d['blocks'])==21 and len(d['summary'])==3
            assert d['source_parse_snapshots_equal']
            assert d['rewrite_sha256']['baseline']==d['rewrite_sha256']['candidate']
            negative[profile]={}
            below_floor[profile]={}
            for mode in ['parse','write','roundtrip']:
                source=d['summary'][mode]['baseline']
                assert source['semantic_prerequisite']
                ratio=source['speedup']
                assert ratio['blocks']==21
                if ratio['median'] < .99:
                    below_floor[profile][mode]=ratio
                if ratio['ci95'][1] < 1:
                    negative[profile][mode]=ratio
                if profile.startswith('sdf-') and (mode!='write' or index==4):
                    primary_pass &= ratio['ci95'][0]>1
                for engine in ['baseline','openbabel']:
                    entry=d['summary'][mode][engine]
                    if entry['strict_speed_win']:
                        assert entry['semantic_prerequisite']
                        assert entry['speedup']['ci95'][0]>1
                        assert entry['speedup']['winning_blocks']==21
            assert below_floor[profile]==decision['below_floor'][profile]
            assert negative[profile]==decision['all_negative_intervals'][profile]
        accepted=primary_pass and not any(below_floor.values())
        assert accepted==decision['local_adoption_gate_passed']
        assert accepted==(index in [0,4])
        audited.append({'name':name,'base_revision':manifest['base_revision'],
                        'patch_reconstruction':True,'complete_profiles':5,'source_lanes':15,
                        'blocks_per_profile':21,'adoption_gate':accepted,
                        'negative_intervals':sum(len(v) for v in negative.values())})
    final=json.loads((ROOT/'benchmarks'/(NAMES[-1]+'-artifacts')/'source-manifest.json').read_text())
    exact=all(sha((ROOT/n).read_bytes())==h for n,h in final['source_files_sha256'].items())
    target=json.loads((ROOT/'benchmarks'/(NAMES[-1]+'-artifacts')/'reconstruction.json').read_text())['changed_source_files']
    assert all(sha((ROOT/n).read_bytes())==final['source_files_sha256'][n] for n in target)
    if '--require-exact-source' in sys.argv:
        assert exact, 'check out measured source checkpoint for whole-source identity'
    print(json.dumps({'passed':True,'current_optimization_files_match':True,'current_whole_source_matches_measured':exact,
                      'no_timing_pooled_or_substituted':True,'packets':audited},indent=2))

if __name__=='__main__':
    main()
