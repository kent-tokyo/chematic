import hashlib, json, subprocess, sys
from pathlib import Path
OUT = Path("validation/results/mmff94-correctness-2026-10-11/optimizer-recovery")
ARGS = ['scripts/diagnostics/mmff94_same_explicit_h_energy.py', '--output', 'validation/results/mmff94-correctness-2026-10-11/optimizer-recovery/candidate.jsonl', '--summary', 'validation/results/mmff94-correctness-2026-10-11/optimizer-recovery/candidate-summary.json', '--rdkit-term-breakdown', '--expected-schematic', '1.0.20', '--expected-rdkit', '2026.03.6', '--schematic-artifact', 'crates/chematic-py/python/chematic/chematic.cpython-313-darwin.so', '--rdkit-artifact', '.venv-v1gate/lib/python3.13/site-packages/rdkit/rdBase.so', '--schematic-build-command', 'VIRTUAL_ENV=/private/tmp/chematic-mmff-development-20261011 .venv-v1gate/bin/maturin develop --release --offline -j 2 -m crates/chematic-py/Cargo.toml -- -C strip=none -C debuginfo=1', '--gradient-input-index', '99', '--gradient-input-index', '94', '--gradient-input-index', '91', '--gradient-input-index', '192', '--gradient-input-index', '222', '--gradient-input-index', '231']
extension=Path("crates/chematic-py/python/chematic/chematic.cpython-313-darwin.so")
before=hashlib.sha256(extension.read_bytes()).hexdigest()
result=subprocess.run([sys.executable,*ARGS],capture_output=True,text=True)
(OUT/"summary.stdout.json").write_text(result.stdout)
(OUT/"energy.stderr.log").write_text(result.stderr)
assert hashlib.sha256(extension.read_bytes()).hexdigest()==before
summary=json.loads((OUT/"candidate-summary.json").read_text())
assert result.returncode==2
assert {k:v for k,v in summary["gate"].items() if not v}=={"source_tree_clean":False}
print(json.dumps({"extension_sha256":before,"gate":summary["gate"]},indent=2))
