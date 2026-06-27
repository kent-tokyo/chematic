# RDKit Issue Lessons

This document explains why chematic maintains a `validation/rdkit_issues/` corpus
and what design principles we derived from studying RDKit's GitHub issues.

## Why this corpus exists

RDKit is the reference implementation for cheminformatics. Its issue tracker is a
living archive of edge cases that real molecules can hit. Rather than copying RDKit's
fixes, we use these issues as:

1. **Regression targets** — verify chematic doesn't reproduce the same failure mode
2. **Design guidance** — understand where leniency causes downstream confusion
3. **Differentiation signal** — areas where a Rust-native, Result-based API can do better

## Key lessons by category

### Canonical SMILES (RDKit #8759, #8775)

**Issue:** `MolToSmiles(MolFromSmiles(MolToSmiles(mol)))` can produce a different string
than `MolToSmiles(mol)` for certain stereocenters and fused ring systems.

**chematic principle:** `canonical_smiles(parse(canonical_smiles(mol))) == canonical_smiles(mol)`
is tested in `crates/chematic-smiles/tests/canonical_robustness.rs` for all SMILES in
`validation/rdkit_issues/stereo/canonical_idempotence.smi`.

### E/Z fragment extraction (RDKit #9368)

**Issue:** `MolFragmentToSmiles()` raises a C++ pre-condition violation when the
fragment boundary is directly adjacent to an E/Z double bond.

**chematic principle:** `brics_fragments()` and `brics_bonds()` must never panic.
Stereo may be dropped silently when it cannot be preserved through a cut — that is
acceptable. Verified in `rdkit_9368_ez_fragment_no_panic` test.

### Atropisomer stereo degradation (RDKit #9338)

**Issue:** `TautomerEnumerator.Canonicalize()` raises on atropisomer-like bond stereo
that cannot survive tautomer enumeration.

**chematic principle:** If stereo cannot be preserved through standardization or
tautomerization, it must be cleared explicitly (not panicked). `parse_smiles_report()`
exposes this as `W002_DROPPED_STEREO`.

### Structured warnings vs stderr (RDKit #2683)

**Issue:** C++-level warnings from large SMILES batches (46M+ molecules) flood stderr
and cannot be captured by Python's `warnings` module.

**chematic principle:** `parse_smiles_report(smiles)` returns `(mol, warnings)` as
structured data instead of writing to stderr. Warnings carry a code (`W001_`, `W002_`,
`W003_`) for programmatic filtering.

### Similarity metric naming (RDKit #8317)

**Issue:** RDKit's `AllBit`, `Asymmetric`, `BraunBlanquet` etc. don't match the names
in textbooks. `AllBit` == Rand, `Asymmetric` == Simpson/Overlap.

**chematic principle:** Functions are named by the standard literature term.
Formula is documented inline or in the API reference.

## Corpus structure

```
validation/rdkit_issues/
  stereo/
    canonical_idempotence.smi     — RDKit #8759 stereo cases
    ez_fragment_extraction.smi    — RDKit #9368 BRICS fragment cases
    atropisomer_invalid_stereo.smi — RDKit #9338 bond stereo degradation
  canonicalization/
    aromatic_kekule_roundtrip.smi  — aromatic ↔ Kekulé stability
    charged_heteroaromatic.smi     — N+/O- near aromatic rings
  fragments/
    ez_near_fragment_bond.smi      — BRICS cuts near E/Z bonds
```

## What we deliberately don't do

- We do not port RDKit's bug fixes verbatim — different architecture.
- We do not claim "RDKit compatibility" as a goal — we aim for chemically correct behavior
  per IUPAC / Daylight / OpenSMILES specs.
- We do not implement every RDKit feature — see `docs/limitations.md` for scope.
