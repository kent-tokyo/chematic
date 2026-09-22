# RDKit accuracy plan

Updated 2026-09-16. The objective is operation-level equivalence on declared
profiles and, where an independent gold standard exists, evidence that chematic
is equivalent to or better than the pinned comparator. RDKit agreement alone is
not chemical truth.

## Evaluation contract

Every reported lane fixes:

- chematic revision/package and binding;
- RDKit version and API/options;
- input corpus, split, hash, and acquisition time;
- operation, preprocessing, tolerance, and ordering rules;
- valid, refused, failed, and timed-out denominators;
- raw output location/hash and summary generator;
- hardware/runtime for performance measurements.

Native and RDKit-compatible profiles are separate. Cross-profile recall is a
diagnostic and must not be reported as compatibility. Canonical-SMILES string
equality is not a substitute for semantic identity.

## Frozen data policy

Development data may be inspected and used for fixes. A sealed holdout may be
consumed only after the candidate is frozen and the source is attested unused.
Once evaluated, it becomes exposed and cannot be reused as the sole evidence for
the next candidate.

The current precondition is recorded in
`validation/results/sealed-cohort-preflight-trust-eval-candidate-20260916.json`:
2,000 development rows and 8,000 sealed rows. No score has yet been calculated.

## Acceptance thresholds

| Operation | Required result on declared valid scope |
|---|---|
| Integer descriptors | exact equality |
| Floating descriptors | absolute error ≤1e-6 unless a narrower published contract applies |
| Fingerprint bits/counts | exact configured output and provenance |
| Retrieval | identical ordered IDs, recall 1.0, unrounded score error ≤1e-12 |
| Canonical/identity | semantic identity, idempotency, no false merge, no information loss |
| Stereo/CIP | zero wrong confident labels; abstentions and unsupported rows reported separately |
| Equivalence against independent gold | paired 95% difference interval inside ±0.1 percentage point, with coverage safeguards |
| Superiority | positive lower confidence bound and no material coverage/safety regression |

The thresholds apply only after parser and support-domain accounting. Dropping
difficult rows cannot improve the compatibility score.

## A0 — Evidence integrity

Complete for the RDKit 2025.09.3 core-eight profile. Candidate `5e9211a6` was
frozen before source acquisition, and its overlap-audited split passed
2,000/2,000 development plus 8,000/8,000 one-time sealed rows for molecular
weight, HBA, HBD, TPSA, LogP, molar refractivity, Fsp3, and aromatic-ring count.
All rows were accounted for; mismatches, parse failures, and unsupported values
were zero. Source, split, binary, evaluator, attestation, and raw-result hashes
are preserved in `validation/results/a0-core-eight-sealed-acceptance-20260922.json`.

Reopen A0 only when the oracle, profile, field set, or tolerance changes. A new
adoption must use a new annotated freeze, post-freeze/unused source, overlap
audit, immutable raw result, and one-time sealed evaluation.

## A1 — Perception and descriptors

Exposed regression lanes show high or exact agreement for core descriptors,
rotatable bonds, and several ring/heteroatom families. Potential stereocenter
count retains known residuals and must not be generalized from the exposed set.

Remaining exit:

- retain the completed core-eight unused result and run unused rows only for
  newly declared descriptor/perception families;
- retain native-profile behavior while testing the opt-in compatibility profile;
- classify every mismatch by perception, normalization, formula, or unsupported
  chemistry;
- rerun affected Rust, Python, and Node/WASM bindings.

## A2 — CIP, canonicalization, and identity

Existing gates cover atom-order invariance, 300 development structures, and
5,115 SMILES spelling variants. Remaining known boundaries include coupled
aromatic E/Z carrier systems and phosphorus-CIP representation instability.

Remaining exit:

- resolve or explicitly refuse all known wrong-confident cases;
- add atom/bond permutations and V2000/V3000 round trips;
- verify canonical idempotency and stable-key no-false-merge behavior;
- adjudicate disputed cases independently of one RDKit spelling.

## A3 — Fingerprints and retrieval

The RDKit-compatible Morgan lane has exact exposed-corpus and cross-binding
evidence for configured folded bits, sparse counts, bitInfo, and top-k search.
Native ECFP4 remains a separate fingerprint definition.

Remaining exit:

- run exact fingerprint and retrieval checks on sealed inputs;
- preserve raw environment/provenance evidence;
- rerun every affected binding after algorithm changes;
- report preprocessing failures and tie ordering explicitly.

## A4 — Workflows and interchange

This package covers SMARTS/substructure, standardization, reaction application,
and V3000. The ordinary-V3000 baseline is complete for its declared external
reader versions; coordination, haptic, polymer expansion, and ENDPTS/ATTACH
semantics are not.

Remaining exit:

- classify and resolve the current SMARTS residual set without partial-output
  promotion;
- measure mapped embedding and reaction-product precision/recall;
- test standardization identity and stereo preservation;
- add typed V3000 capabilities only with external-reader round trips.

## A5 — Independent adjudication

RDKit parity cannot establish superiority over RDKit. Prepare an adjudication
packet containing predeclared metrics, blinded cases, raw outputs, disagreements,
and reproduction instructions. Completion requires absolute gold labels and a
non-maintainer review. Until then, use “RDKit-compatible on the declared lane,”
not “more accurate than RDKit.”

## A6 — 3D and force fields

Treat 3D as a separate experimental profile. Evaluate MMFF94/UFF atom typing,
charges, parameters, energy, analytic gradients, convergence, timeouts, stereo
retention, and conformer quality on identical coordinates and seeds. Bounded
typing or finite energy alone does not establish force-field or ETKDG parity.

## Promotion workflow

1. Fix and test on development data.
2. Freeze an annotated candidate and dependency lockfile.
3. Verify the unused-data attestation and split hashes.
4. Run the sealed lane once and preserve raw output.
5. Generate operation-level reports and classify failures.
6. Obtain independent review where A5 claims are involved.
7. Update the compatibility dashboard and release notes.

Any tuning after step 4 requires a new candidate and a new unused holdout.
