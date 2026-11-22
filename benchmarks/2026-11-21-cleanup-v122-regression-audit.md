# Cleanup v1.2.2 regression audit

The v1.2.2 source comparison against RDKit 2026.09.1 matches 4,022 of 4,069
comparable stress strings and 9,996 of 10,000 broad strings. Compared with
the October 24 baseline, 47 previously exact rows differ in the stress
corpus and two in the broad corpus.
[Stress and sanitization data](../validation/results/rdkit-cleanup-sanitized-stress-2026-11-21.json) ·
[Broad data](../validation/results/rdkit-cleanup-sanitized-10k-2026-11-21.json)

## Published wheel and candidate

| Check | Published v1.2.2 wheel | Local candidate `7fc29a3e` |
|---|---|---|
| The 47 selected rows | 31 graph differences from RDKit. The other 16 agree after sanitization. Row 2025 loses an input hydrogen in both libraries. | 45 outputs are unchanged. Rows 2025 and 2203 retain their input metal hydrogen. All 47 produce outputs, conserve input formula and retain the same string on repeat Cleanup. |
| Rooted SMILES for the 31 graph differences | One chematic output per molecule across 797 spellings. RDKit gives multiple outputs on 24 molecules. | The candidate retains one output per molecule and the same result on repeat Cleanup. |

[Published atom audit](../validation/results/cleanup-v1.2.2-47-atom-audit-2026-11-21.json) ·
[Candidate atom audit](../validation/results/cleanup-47-metal-ligand-stable-audit-2026-11-21.json) ·
[Rooted inputs](../validation/results/cleanup-rooted-smiles-31-v1.2.2-2026-11-21.json)

The candidate source matches 4,000/4,069 stress strings and 9,994/10,000
broad strings. Those denominators include seven and two refusals respectively.
Reionize matches 4,072/4,072 and 10,000/10,000. The separate 575-row metal
audit finds identical outcomes for compact and explicit-H pairs, with seven
paired refusals. Row 2796 still changes the hydrogen total under the existing
isolated-potassium/phenol policy.
[Candidate stress](../validation/results/rdkit-cleanup-metal-ligand-stable-stress-2026-11-21.json) ·
[Candidate broad](../validation/results/rdkit-cleanup-metal-ligand-stable-10k-2026-11-21.json) ·
[Metal audit](../validation/results/cleanup-metal-stress-metal-ligand-stable-2026-11-21.json)

The audit measures output agreement, formula conservation and repeatability on
fixed inputs. Adjudication of the remaining chemical differences is open.
The [full dated experiment log](https://github.com/kent-tokyo/chematic/blob/0be29688d88c13baebc323418a5f4143cd31836d/benchmarks/2026-11-21-cleanup-v122-regression-audit.md)
retains all intermediate candidates, reverted experiments and raw records.
