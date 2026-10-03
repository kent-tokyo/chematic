# Published v1.0.30 output-audit decision

The P0.1 task was to rerun, account for, and compare **published artifacts**.
That measurement audit is complete on the exposed corpus. It is not a claim
that v1.0.30 has full RDKit compatibility. The
[artifact record](2026-10-02-v1.0.30-published-artifact-gates.md) and
`python3 scripts/check_v1030_artifact_packet.py` pin the wheel, npm tarball,
crates.io graph, RDKit 2026.03.6, input and raw-output hashes, and row totals.
The source-only formula fix and older timing records are excluded.

| Result | Adoption policy |
|---|---|
| 10,000 chemistry inputs and 310,000 SMARTS cells per applicable artifact | Count every input/cell. The 200 different SMARTS **match sets** are compatibility failures, not matches or typed refusals. Of these, 43 also change the hit/no-hit Boolean; 157 have the same Boolean but different matched atom sets. All are named in the archived classification: 194 symmetrized-ring and six organometallic ring-semantics cells across 64 rows. Native SSSR remains the default; a future RDKit-style ring profile must be opt-in and newly measured. |
| CIP and compatible Morgan on the same 10,000 inputs | CIP is 9,995 exact plus five **typed abstentions**, not 100% labelled parity. Morgan is 9,999 exact plus one typed refusal, not 10,000 exact. Preserve both denominators and reasons. |
| Canonical SMILES | Literal RDKit spelling agrees on 64/10,000, while semantic round-trip succeeds on 10,000/10,000. Treat semantic identity and exact spelling as separate gates. |
| 63-operation matrix, v1.0.29 → v1.0.30 | Published Python and Rust each cover all 63 operations and 210,410 output rows per release. Only HBA and its Lipinski bundle change, on the same 1,359 rows. The 10,000-row chemistry stream does not include HBA; its byte identity is not evidence of HBA parity. HBA's independent 5,000/5,000 RDKit result is separate. |
| npm 63-operation slice | 59 equivalent public adapters exist. Four operations have **no equivalent exported API**: N/O-only RDKit TPSA, RDKit-compatible Atom Pair, RDKit-compatible Pattern fingerprint, and MMFF94 minimization of explicit ETKDG coordinates. These are API gaps, not zero-valued results or parity passes. |
| npm/Python cross-binding output | 52 adapted operations are exact on their applicable rows. QED/Chi1v last-bit differences use the recorded `1e-12` tolerance; 1,575 formula rows have equal tokenized composition but different spelling; 100 ETKDG coordinates differ only within 0.000051 from four-decimal serialization. These are explicitly scoped numeric/representation classifications, not byte identity. MOL writer (2,000 rows), SVG (20/300), and 2D layout (500 rows) remain **non-equivalent representations**; do not credit a speed or output-parity win for them. |
| Legacy 57 reactions | The frozen 57/57 product-set gate remains intact. The separate 83-case extension still has confident published-artifact differences and belongs to P1, not a completed P0.1 reaction-parity claim. |

The v1.0.29/v1.0.30 rowwise differential was executed independently for
each published binding and does not inherit v1.0.29 measurements. The
200 SMARTS cells cannot be relabelled as unsupported after the fact: the
published APIs returned confident match sets. Their chemistry policy and future
implementation work remain open under P1. The four npm gaps likewise require
a future published package before they can be measured. Exposed fixtures do
not become sealed evaluation data.

**Decision:** close P0.1 as an artifact-rebaseline and accounting task with
explicit failures/gaps. Keep strict SMARTS, full CIP, all-API cross-binding,
and broad reaction compatibility **open**. P0.2 speed and P1 reaction exits
are independent and cannot be inferred from this audit.
