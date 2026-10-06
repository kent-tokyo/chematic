# v1.0.36 published reruns, cache, write_mol, RDKit.js 2026.09.1 reactions

Follows [batch 10](2026-10-06-734-754-followups-batch10.md). v1.0.36 is the
first release with the batch 6–10 fixes; this record reruns their evidence
on the published artifacts and closes four follow-ups. One Linux x86-64
host, RDKit 2026.03.6 (CPython 3.11), Node 22. Source lane: code revision
`e4e2aa73` (the evidence commit that follows changes documents and result
files only).

## Published v1.0.36 artifacts

| Artifact | Hash |
|---|---|
| PyPI `chematic-1.0.36-cp39-cp39-manylinux_2_17_x86_64.manylinux2014_x86_64.whl` | sha256 `5d897d8b…4860cf` |
| PyPI `chematic-1.0.36.tar.gz` | sha256 `8fa030ab…475673` |
| npm `@kent-tokyo/chematic@1.0.36` | `sha512-8abA6x6A…hqZD0w==` |
| crates.io `chematic` 1.0.36 | checksum `0a074112…983d61f` |

The Linux wheel is a CPython 3.9 build and RDKit 2026.03.6 ships no
CPython 3.9 wheel, so the published wheel ran in its own interpreter and
RDKit in another: `scripts/chematic_reaction_worker.py` (driven by
`biotransformer_rule_corpus.py --chematic-python`) for reactions and
`scripts/chematic_chemistry_dump.py` with `compare_chemistry_dump_rdkit.py`
for the rest.

| Check | Published v1.0.36 | Source before this record |
|---|---|---|
| BioTransformer corpus (974 rules × 400, RDKit profile) | implicit-H 1,051/1,051, explicit-H 6,604/6,652 same sets, 14 typed re-sanitize rows; **all 233 non-exact rows identical to source** ([summary](../validation/results/v1.0.36-published-pypi-linux-cp39-biotransformer-corpus.json), [rows](../validation/results/v1.0.36-published-pypi-linux-cp39-biotransformer-corpus-rows.jsonl)) | same |
| 83 reaction fixtures, checked provenance: PyPI wheel / npm / crate | 80 graph/origin/map + 3 jointly invalid on each ([pypi](../validation/results/v1.0.36-published-pypi-linux-cp39-reaction-checked-provenance-83.json), [npm](../validation/results/v1.0.36-published-npm-reaction-checked-provenance-83.json), [crate](../validation/results/v1.0.36-published-crate-reaction-checked-provenance-83.json)) | same |
| xsmarts-autoconf lane (`6ce94c9`), sdist built for CPython 3.11 | 75/83, the same 8 expected differences, 0 unexpected ([json](../validation/results/v1.0.36-published-sdist-xsmarts-autoconf-lane.json)) | same |
| SMARTS opt-in profile, 310,000 pinned cells (archived oracle) | 309,982 exact, 18 typed refusals, 0 unexpected ([json](../validation/results/v1.0.36-published-pypi-linux-cp39-smarts-optin-310k.json)) | same |
| Accurate CIP vs `rdCIPLabeler` | exposed 10k 4,394 agree + 1 abstention; ChEMBL 5k 4,186 agree (20 marked Kekulé-dependent) | same |
| Hybridization vs `GetHybridization` | 220,015/220,015 and 138,655/138,655 atoms | same |
| MMFF94 heavy atoms typed differently | 81 (10k) and 5 (5k); hydrogens 0 | same |
| MOL writer, RDKit reads the block as the input | 1,686/1,687 and 1,669/1,670; `strict=True` refuses 1 row each | same |
| WASM Node tests (`crates/chematic-wasm/tests`, 20 files) on the npm tarball | all pass, the Kekulé-dependent CIP JSON included | — |

[Chemistry summary](../validation/results/v1.0.36-published-pypi-linux-cp39-chemistry-vs-rdkit.json).
CIP (both engines, E/Z included) and MOL blocks of the published wheel and
the source after this record are identical on all 15,000 rows. macOS and
Windows wheels were not run on this host.

## Release channels and Pages

The v1.0.36 Pages run failed on a relative link out of `docs/` that MkDocs
cannot resolve; [#765](https://github.com/kent-tokyo/chematic/pull/765)
fixed it and [#766](https://github.com/kent-tokyo/chematic/pull/766)
recorded the release channels. Independently re-measured here: PyPI 17
wheels and 1 sdist, the npm integrity above, all 20 crates on
static.crates.io, a non-draft Release with its metadata asset, every tag
workflow successful, and the Pages run after #765 successful.
`scripts/check_release_docs_consistency.py` (run in CI) now fails on any
link under `docs/` that leaves the MkDocs tree, so the next one fails CI
rather than the Pages deploy; `mkdocs build --strict` passes.

## SMIRKS template cache

`PreparedReaction::shared` (used by every SMIRKS-text entry point) kept up
to 2,048 prepared templates and cleared them all when full, so a working
set was prepared again after every 2,048 new templates. It is now 16
shards per reading, each evicting its least recently used template; a hit
no longer allocates its key.

| Load | Clear-all (v1.0.36) | Sharded LRU |
|---|---:|---:|
| 300 templates used every round + 5,000 one-off, 100 rounds: preparations | 5,900 | **5,300** (each once) |
| 40,000 hits on one template (time) | 4.5–6.2 ms | **3.0 ms** |
| 8 threads × 1,000 shared templates × 5 rounds: preparations | — | 1,000–1,100 (test bound) |

`crates/chematic-rxn/tests/prepared_cache.rs` checks the first and third
rows. The host has two cores, so no contention figure is claimed.

## write_mol

| `to_mol_block`, exposed 10k (Python, best of 5) | All rows | 1,687 stereo rows |
|---|---:|---:|
| v1.0.34 (no coordinates; stereo not written) | 18.1 ms | 5.8 ms |
| v1.0.35 | 85.3 ms | 68.3 ms |
| v1.0.36 | 78.7–80.6 ms | 60.4–63.5 ms |
| source | **70.3–72.3 ms** | **51.8–52.0 ms** |

Since v1.0.35 the writer lays out and wedges every stereo molecule, so it
writes 1,686 of 1,687 of them back for RDKit instead of none; the other
8,313 rows cost what they did in v1.0.34. The source skips the CIP ranking
and small-ring search for double bonds without `/`/`\` markers (E/Z
perception everywhere, not only here), reuses the molecule's memoized SSSR
in the layout, hashes layout sets with FxHash and walks bonds without a
per-step allocation: 12% fewer instructions writing 600 stereo molecules,
byte-identical on 15,000 blocks. The rest is the layout's ring perception
and the canonical order that fixes the Kekulé form.

## RDKit.js 2026.09.1 reactions

RDKit 2026.09.1 has no Python wheel, but RDKit.js 2026.09.1 is published
with `get_rxn` / `run_reactants`. `scripts/rdkitjs_reaction_corpus.mjs`
(one rule per process, `rdkitjs_reaction_corpus_driver.py`) ran the 974
single-reactant BioTransformer rules on the same 400 reactants, both
hydrogen modes, with RDKit.js 2026.03.6 and 2026.09.1, sanitizing each raw
product by reading its MOL block back:

| | Rows |
|---|---:|
| (rule, reactant, hydrogen mode) rows with raw products, or aborted | 8,706 |
| same sanitized product sets in RDKit.js 2026.03.6 and 2026.09.1 | **8,627** |
| differ, and 2026.09.1 now equals the Python 2026.03.6 oracle | 70 |
| differ: E/Z perceived from generated coordinates on double bonds the input leaves unmarked | 8 |
| aborted the WebAssembly instance in both builds (BTMR0654 on one explicit-H reactant) | 1 |

[Differential](../validation/results/rdkitjs-2026.03.6-vs-2026.09.1-biotransformer-reactions.json).
Every one of the 78 changes is stereo or hydrogen spelling read back from
the MOL block, not a different product: 2026.09.1 reads back 70 rows (a
peptide's histidine centre, H-pinned radicals) as the Python harness
writes them, and its new 2D layout moves 8 unmarked double bonds to the
other side. The read-back itself is why RDKit.js agrees with the Python
2026.03.6 oracle on 8,448 (2026.03.6) and 8,518 (2026.09.1) of the 8,705
rows rather than all of them. No reaction-semantics change between the two
releases shows on this corpus, so chematic's 2026.03.6 profile needs no
re-baselining here; the Python 2026.09.1 wheel, when published, remains
the stronger check.
