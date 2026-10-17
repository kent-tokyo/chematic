# Open Babel 2D MOL identity regression: cumulated ring bonds

On the fixed 1,216-identity stereo stress selection, the CheMatic CLI built
from main `a5c45ee6410e80d8043083247e01e06c4863ea67` lost three molecules
in both V2000 and V3000. The preceding [main audit](main-stereo-extended-2026-10-17.md)
had zero Open Babel-only matches. The loss followed the CLI's new CTAB
aromaticity normalization: its ring model counted a carbon with two
endocyclic double bonds as an ordinary 1π aromatic center, changing a
nonaromatic cumulated ring into aromatic radical atoms. This patch excludes
that carbon from the ring's aromatic circuit while retaining normalization
of genuinely aromatic CTABs.

RDKit 2026.09.1 wrote and reread every 2D MOL input. Rows enter the comparison
only when the requested MOL version and source identity survive that gate.
CheMatic and Open Babel 3.2.1 then receive identical MOL bytes. RDKit reads
their output SMILES and compares identity with the source, including specified
stereochemistry. The measured host was macOS arm64.

| MOL input | Comparable | CheMatic before | CheMatic after | Open Babel | Source exclusions |
| --- | ---: | ---: | ---: | ---: | ---: |
| V2000 | 1,157 | 1,154 | **1,157** | 974 | 59 |
| V3000 | 1,192 | 1,189 | **1,192** | 936 | 24 |

The only changed rows are corpus lines 1476, 1643, and 1647 in each format:
six failures became strict source-identity matches. All 2,432 selected
format attempts retained the same source-gate status and identity; there
were no CheMatic regressions. The CLI's existing aromatic CTAB test also
passes, covering pyrrole and charged aromatic rings in both formats.

The [before packet](results/openbabel-stereo-2d-stress-1216-regression-2026-10-18.json.gz)
has compressed SHA-256
`e5b3850ab5860dd971dbe861fff4bb7a5e7808a87cec1ef7d4d6338819e6292e`.
The [after packet](results/openbabel-stereo-2d-stress-1216-cumulated-ring-2026-10-18.json.gz)
has compressed SHA-256
`c49ca2e270baed322c19f5e96b161e424ebcf9a293c7b63e4744d61e8d696487`.
Both preserve per-row outputs, input MOL hashes, exclusions, classifications,
and tool metadata. The tested after-patch CLI SHA-256 is
`3ad08f26a1a9aeecd98c11596710d5a65fed8304db29b9069fa144770d064803`;
Open Babel's executable SHA-256 is
`733b2c925d60b17a5c9d8ce9ec841e918cc20c7634301e23f9008a187d1d9bc5`.
The corpus SHA-256 is
`879166ab5481cac805ecd2d77f37d628383ec124085b90e2ee1659d55ce76277`.

## Reproduce

With RDKit 2026.09.1 and Open Babel 3.2.1 installed, build the CLI and run:

```sh
cargo build -p chematic-cli
python3 scripts/compare_openbabel_stereo_formats.py \
  --cli target/debug/chematic --obabel /opt/homebrew/bin/obabel \
  --corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --stress-first-unique 1216 --output /tmp/openbabel-2d-1216.json
```

These results establish the CheMatic CLI's lead on this paired 2D MOL-to-SMILES
workflow and restore the result measured before the regression. RDKit source
identity is an operational oracle, not independent chemical truth. This does
not establish toolkit-wide accuracy superiority, published-package behavior,
or results for other input formats and chemical classes.
