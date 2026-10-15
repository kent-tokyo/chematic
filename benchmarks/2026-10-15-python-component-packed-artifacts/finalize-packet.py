from pathlib import Path
import json,hashlib,shutil,importlib.util
r=Path('/private/tmp/chematic-component-packed-20261015');d=r/'delivery';a=d/'benchmarks/2026-10-15-python-component-packed-artifacts';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
p=r/'package.py';t=p.read_text().replace('Every operation, mode, block, error and mismatch is retained.','Every successfully checkpointed operation, mode, block, error and mismatch is retained; the unavailable observation is disclosed below.').replace("if p.is_file()};(a/'artifact-sha256.json')", "if p.is_file() and p.name!='artifact-sha256.json'};(a/'artifact-sha256.json')");p.write_text(t)
report=d/'benchmarks/2026-10-15-python-component-packed.md';t=report.read_text().replace('Every operation, mode, block, error and mismatch is retained.','Every successfully checkpointed operation, mode, block, error and mismatch is retained; the unavailable observation is disclosed below.');report.write_text(t)
for name in ['package.py','deliver.py','package.log','finalize-packet.py']:shutil.copyfile(r/name,a/name)
s=Path('/private/tmp/chematic-component-packed-source-20261015');m=json.loads((r/'source-manifest.json').read_text());assert all(sha(s/n)==h for n,h in m['candidate_sha256'].items())
for n in m['changed_production_files']+m['changed_documentation_files']+m['added_test_files']:assert sha(d/n)==m['candidate_sha256'][n]
index={p.relative_to(a).as_posix():sha(p) for p in sorted(a.rglob('*')) if p.is_file() and p.name!='artifact-sha256.json'};(a/'artifact-sha256.json').write_text(json.dumps(index,indent=2)+'\n');assert all(sha(a/n)==h for n,h in index.items())
outer={p.relative_to(d).as_posix():sha(p) for p in sorted(d.rglob('*')) if p.is_file()};(r/'delivery-sha256.json').write_text(json.dumps(outer,indent=2)+'\n');assert all(sha(d/n)==h for n,h in outer.items())
spec=importlib.util.spec_from_file_location('names',s/'scripts/check_product_names.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);assert not module.scan_paths(d,[report,a/'README.md',d/'CHANGELOG.md'])
print('Final packet hashes verified:',len(index),'non-self-indexed artifacts;',len(outer),'delivery files; measured source exact')
