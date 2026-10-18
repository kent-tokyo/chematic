# RDKit workflow follow-up: reaction, Cleanup, and R-group

Date: 2026-10-18. Base commit: `f99c7b8ab33098a2e1fd9fd568aee370bf79ab38` plus this branch's source changes. Oracle: RDKit `2026.09.1`. These are source-candidate measurements, not published-package results or general chemical accuracy claims.

## Cleanup

The source-built `cleanup_rdkit_dump` example and RDKit were rerun on all 4,072 stress-corpus rows. The two parsers agreed on the input for all rows; RDKit Cleanup accepted 4,069. Exact output agreement rose from **3,885 to 3,891/4,069 (95.63%)**. Six formerly different rows now agree, with **zero regressions** on the same scored set. The seven prior CheMatic output errors are now zero; one of those seven still produces a valid but different structure. The change preserves aromatic N-oxide `N+-O-` bonds instead of converting them to invalid `N=O` valence. [Row-level result](2026-10-18-rdkit-cleanup-n-oxide-followup.json).

The remaining **178** differences include metal handling, acid/base priority, and other representation changes. Exact agreement with RDKit is a compatibility measure, not independent chemical truth.

## Labelled R-group

Terminal mapped wildcards can now match an implicit hydrogen when the corresponding substituent is absent. A hydrogen R-group is emitted only if another row establishes that column. Symmetric partial substitution prefers the lower available label for the real substituent. The fallback is limited to single-bond terminal wildcards at distinct core anchors and retains the 1,024-match cap. The source-built Python extension matches RDKit core and R-group graphs with matching labels on **8/8** small paired cases, up from the previous 4/4 mapped-core-atom cases. [Case-level result](2026-10-18-rdkit-rgroup-hydrogen-followup.json).

## Reaction

The same rebuilt extension was rerun against the pinned 83-case checked-product gate: **80** exact graph, atom-origin, and template-map outcomes, plus **three** inputs invalid in both engines. The separate opposite-enantiomer case still refuses the mismatched reactant in CheMatic (**1/1**), whereas RDKit's raw runner emits a product set. [83-case result](2026-10-18-rdkit-reaction-83-followup.json) · [stereo safety result](2026-10-18-rdkit-reaction-stereo-safety-followup.json). Reaction source was unchanged. These fixtures do not measure yield, selectivity, or general SMIRKS precision and recall.

## Reproduce

```sh
cargo build -p chematic-chem --example cleanup_rdkit_dump --offline
target/debug/examples/cleanup_rdkit_dump \
  < validation/cosmolkit_comparison/corpora/stress-v1.smi \
  > /tmp/chematic-cleanup-followup.jsonl
python scripts/rdkit_cleanup_source_gate.py \
  --dump /tmp/chematic-cleanup-followup.jsonl \
  --previous-source benchmarks/2026-10-18-rdkit-cleanup-source-gate.json \
  --output /tmp/chematic-cleanup-followup-gate.json
python scripts/rdkit_rgroup_hydrogen_gate.py \
  --output /tmp/chematic-rgroup-followup-gate.json
python scripts/reaction_python_checked_provenance_gate.py \
  --expected-rdkit 2026.09.1 --output /tmp/chematic-reaction-followup.json
python scripts/reaction_product_stereo_safety_gate.py \
  --output /tmp/chematic-reaction-stereo-followup.json
```

The Python commands require the source-built CheMatic extension and RDKit `2026.09.1` in the same interpreter. The machine-readable results record source and binary hashes where applicable.
