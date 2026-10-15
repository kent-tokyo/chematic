from pathlib import Path
import json,hashlib,shutil,os,ctypes
r=Path('/private/tmp/chematic-basis-roots-20261015');src=r/'delivery';dst=Path('/Users/k_tanabe/.codex/worktrees/relevant-cycle-speed/chematic');m=json.loads((r/'source-manifest.json').read_text());index=json.loads((r/'delivery-sha256.json').read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();adopt=json.loads((r/'decision.json').read_text())['adopt_source'];changed=(m['changed_documentation_files']+m['changed_production_files']) if adopt else [];new=m['added_test_files'] if adopt else []
lib=ctypes.CDLL('/usr/lib/libSystem.B.dylib',use_errno=True);clone=lib.clonefile;clone.argtypes=[ctypes.c_char_p,ctypes.c_char_p,ctypes.c_int]
def copy_guarded(origin,temporary):
 assert not temporary.is_symlink()
 if temporary.exists():temporary.unlink()
 if clone(os.fsencode(origin),os.fsencode(temporary),0)!=0:shutil.copyfile(origin,temporary)
for n,h in m['baseline_source']['sha256'].items():
 expected={h,index[n]} if n in changed else {h}
 assert sha(dst/n) in expected,'Attached baseline changed: '+n
for n in index:
 assert n in changed+new or n=='benchmarks/2026-10-15-python-basis-roots.md' or n.startswith('benchmarks/2026-10-15-python-basis-roots-artifacts/'),n
 assert sha(src/n)==index[n],n
 target=dst/n
 assert not target.is_symlink(),n
 if target.exists() and n not in changed:assert sha(target)==index[n],'New destination differs: '+n
for n,h in index.items():
 target=dst/n;target.parent.mkdir(parents=True,exist_ok=True);tmp=target.with_name(target.name+'.basis-roots-copy');mode=(target if target.exists() else src/n).stat().st_mode & 0o777;copy_guarded(src/n,tmp);tmp.chmod(mode);assert sha(tmp)==h,n;os.replace(tmp,target);assert sha(target)==h,n
for n,h in m['baseline_source']['sha256'].items():
 if n not in changed:assert sha(dst/n)==h,'Unrelated attached file changed: '+n
(r/'delivery-verification.json').write_text(json.dumps({'delivered_sha256':index,'unchanged_baseline_files':len(m['baseline_source']['sha256'])-len(changed),'applied_exact_measured_patch_on_previous_source_candidate':adopt,'baseline_source_unchanged':not adopt,'no_shared_git_mutation':True,'adopted_source':adopt,'merged_or_published':False},indent=2)+'\n')
print('Delivered',len(index),'files; all',len(m['baseline_source']['sha256'])-len(changed),'unrelated frozen source/input files unchanged; no Git mutation')
