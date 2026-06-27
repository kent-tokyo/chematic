# RDKit Issue-Inspired Regression Corpus

SMILES fixtures drawn from RDKit GitHub issues, used as chematic regression tests.

**Goal:** not to fix RDKit bugs, but to verify chematic doesn't reproduce the same failure modes.

## Categories

| Directory | Coverage |
|---|---|
| `stereo/` | E/Z bond fragment extraction, atropisomer stereo degradation, canonical stereo parity |
| `canonicalization/` | Idempotence (canonical² == canonical), aromatic↔Kekulé roundtrip, charged heteroaromatics |
| `fragments/` | BRICS fragment extraction near E/Z bonds — never panics |

## Fixture format

```
# RDKit #NNNN: short description
# Expected chematic behavior: <preserve | drop_with_warning | parse_and_warn | never_panic>
SMILES  name
```

## Expected behaviors

- `preserve` — stereo / metadata must survive the round trip
- `drop_with_warning` — lossy but explicit; warning returned, no crash
- `never_panic` — correctness not guaranteed, but must not panic or infinite-loop
- `idempotent` — canonical(canonical(smi)) == canonical(smi)

## Running

```bash
# Quick idempotence check across all fixtures
python3 -c "
import chematic, pathlib
for f in pathlib.Path('validation/rdkit_issues').rglob('*.smi'):
    for line in f.read_text().splitlines():
        if line.startswith('#') or not line.strip(): continue
        smi = line.split()[0]
        try:
            mol = chematic.from_smiles(smi)
            c1 = mol.smiles
            c2 = chematic.from_smiles(c1).smiles
            assert c1 == c2, f'NOT idempotent: {f.name}: {smi} -> {c1} -> {c2}'
        except Exception as e:
            print(f'ERROR {f.name}: {smi}: {e}')
print('corpus check done')
"

# Rust tests
cargo test -p chematic-smiles --test canonical_robustness --quiet
```
