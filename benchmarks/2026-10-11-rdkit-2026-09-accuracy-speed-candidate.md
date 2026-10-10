# RDKit 2026.09.1 accuracy and speed candidate

This record covers source commit `ecad2162` on macOS arm64 with Python 3.13.6.
It is candidate evidence, not a claim about a published chematic package.

## Scope

- Corpus: `validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi`
- Corpus SHA-256: `f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`
- Reference: RDKit 2026.09.1
- Secondary comparator: COSMolKit 0.5.0rc22
- Accuracy: all 10,000 rows for TPSA, Labute ASA, and ring count
- Speed: first 1,000 rows, 21 alternating-order blocks, medians and paired
  bootstrap 95% confidence intervals

## Accuracy

| Operation | Result against RDKit 2026.09.1 |
|---|---:|
| TPSA | 10,000 / 10,000 exact |
| Ring count | 10,000 / 10,000 exact |
| Labute ASA | 8,863 / 10,000 bit-identical; 10,000 / 10,000 within `1e-9` |

RDKit 2026.09 changed its default ring backend. Before this change, chematic
matched 9,994 / 10,000 ring counts; the six residuals were large symmetric
macrocycles. The candidate routes those ring families to relevant-cycle
enumeration while retaining the faster legacy-equivalent path for compact
ring systems.

The checked summary is
[`../validation/results/rdkit-2026-09-selected-accuracy-candidate.json`](../validation/results/rdkit-2026-09-selected-accuracy-candidate.json).

## Speed against RDKit

Ratios above 1 mean chematic is faster.

| Operation | chematic median | RDKit median | Direct ratio (95% CI) | Winning blocks |
|---|---:|---:|---:|---:|
| TPSA | 1.451 ms | 2.191 ms | 1.501× (1.493–1.535) | 21 / 21 |
| Chiral Morgan R2 | 14.022 ms | 15.753 ms | 1.121× (1.116–1.130) | 21 / 21 |
| Labute ASA | 1.377 ms | 1.068 ms | 0.773× (0.764–0.780) | 0 / 21 |
| Ring count | 3.333 ms | 0.077 ms | 0.023× (0.022–0.024) | 0 / 21 |

RDKit computes ring information while parsing, whereas chematic keeps it lazy.
The direct ring-count lane therefore compares an RDKit cache lookup with
chematic's first computation. The parse-inclusive lane below is the relevant
end-to-end comparison; the direct deficit remains open rather than being
hidden by moving work into parsing.

## Parse-inclusive pipelines

All four pipelines beat RDKit in all 21 blocks:

| Pipeline | Median speedup (95% CI) |
|---|---:|
| Parse + TPSA | 11.906× (11.870–11.984) |
| Parse + Labute ASA | 11.870× (11.833–11.922) |
| Parse + ring count | 8.147× (8.109–8.171) |
| Parse + chiral Morgan R2 | 4.000× (3.972–4.009) |

The raw timing record is
[`2026-10-11-rdkit-2026-09-accuracy-speed-candidate.json`](2026-10-11-rdkit-2026-09-accuracy-speed-candidate.json).

## Decision

Adopt the ring-accuracy correction and the output-neutral no-closure bond-order
hot path. Do not claim complete speed superiority: direct Labute ASA and the
first direct ring-count call remain slower. The next speed work should target
those two lanes without moving hidden work into every parse.
