# Cleanup v1.2.2 regression audit

The fixed stress corpus, RDKit 2026.09.1 and v1.2.2 source dump reproduce 4,022 exact Cleanup strings among 4,069 compatible rows. The 47 previously exact rows have the following classification. The published macOS arm64 Python 3.13 wheel was checked against its archive before inspecting its in-memory results. The wheel output matches the source record on all 47 rows.

| Recorded group | Chemical difference | Same after sanitization | Refused |
|---|---:|---:|---:|
| Acid/base | 15 | 0 | 0 |
| Metal charge | 2 | 0 | 0 |
| Metal disconnection | 14 | 16 | 0 |
| Total | 31 | 16 | 0 |

For the 31 chemical differences, element, formal charge, attached hydrogen and connectivity do not match under graph isomorphism. For the other 16, isotope and bond graphs match, and sanitizing a copy of RDKit's Cleanup product yields the same canonical isomeric SMILES as reparsing the chematic output. This result describes the sanitized products. The original RDKit product's radical state differs from the state inferred on parsing the chematic SMILES, and is retained in the row record. The graph and repeat-Cleanup checks are in [cleanup-v1.2.2-47-atom-audit-2026-11-21.json](../validation/results/cleanup-v1.2.2-47-atom-audit-2026-11-21.json).

All 47 chematic outputs have the same recorded compatibility SMILES after a second Cleanup. RDKit's compatibility SMILES changes after a second Cleanup on 39 of the rows. The eight chemical differences whose RDKit output stays the same on a second Cleanup are indices 65, 83, 101, 809, 911, 1664, 2076 and 2855. Stability alone does not establish the correct chemical result.

For nine of the 31 chemical differences, the second RDKit Cleanup string equals canonical SMILES obtained by reparsing the chematic output. The row record retains this separate check; it does not establish that the first in-memory RDKit product is the same molecule.

For each of those eight inputs, [cleanup-rooted-smiles-v1.2.2-2026-11-21.json](../validation/results/cleanup-rooted-smiles-v1.2.2-2026-11-21.json) records deterministic SMILES rooted at each atom. The two libraries parse every variant to the same canonical input. RDKit Cleanup returns between two and nine distinct strings per input, depending on atom order. The published chematic wheel returns one string per input. Choosing a single input ordering to imitate one RDKit result would change this observed representation invariance.

An experiment using RDKit's first proposed acid/base transfer fixed 24 stress-string differences and both prior regressions in the 10,000-row corpus. Four existing tests failed, covering acid strength, Kekulé and aromatic forms, a zinc-cut conjugate base, and ferrocene ions. The experiment was reverted. Input-order and query-embedding-order experiments introduced new mismatches and were also reverted. Their rows, counts and source hashes are in [cleanup-order-experiments-2026-11-21.json](../validation/results/cleanup-order-experiments-2026-11-21.json) and its linked reports. The remaining chemical differences still need a fix or an explicit supported-policy decision.

The updated source gate preserves exact-string counts and separately records equivalence after sanitizing RDKit's product. On the stress corpus, 16 of 47 mismatches pass that additional check. On the 10,000-row corpus, one of four mismatches passes. The exact-string counts remain 4,022/4,069 and 9,996/10,000, with respectively 47 and two previously exact rows differing from the 2026-10-24 source comparison. Raw Reionize matches 4,072/4,072 and 10,000/10,000. Results and corpus hashes are in [rdkit-cleanup-sanitized-stress-2026-11-21.json](../validation/results/rdkit-cleanup-sanitized-stress-2026-11-21.json) and [rdkit-cleanup-sanitized-10k-2026-11-21.json](../validation/results/rdkit-cleanup-sanitized-10k-2026-11-21.json).
