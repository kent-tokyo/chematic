# Descriptor Census — `crates/chematic-chem/src/descriptors.rs` (Agent E, diagnosis only)

**Status:** diagnosis complete. No production code changed. Not merged.

**Branch:** `diag/descriptor-census` (forked from `main`@`659baca221f71f135ce0e1780e71245d8770f132`)

**Files touched (all diagnosis tooling, none under `crates/*/src/**`):**
- `docs/descriptor_census_rfc.md` (this file)
- `scripts/descriptor_census.py` (re-runnable Python census script)
- `scripts/descriptor_census_corpus.smi` (5,000-molecule ChEMBL corpus, committed for reproducibility)
- `crates/chematic-chem/examples/descriptor_census_unbound.rs` (new example, not `src/` — dumps the handful of `descriptors.rs` values unreachable from Python)
- `validation/results/descriptor_census.json` (machine-readable summary, full per-descriptor breakdown)
- `validation/results/descriptor_census_unbound.jsonl` (raw Rust-side dump backing the above)

**Explicitly out of scope, not touched:** everything under `crates/*/src/**` (no production fixes, even for the several concrete bugs found below); every other agent's branch (`feat/io-mrv`, `feat/io-tdt`, `feat/io-smiles-supplier-writer`, `fix/smiles-bracket-implicit-h`, `diag/stereo-reader-integration-boundary`, `feat/stereo2d-local-parity`, `diag/canonical-smiles-residual`, `diag/aromaticity-rdkit-parity`, `diag/etkdg-3d-gap`).

**Deliverables:** this RFC, the census script, the corpus file, the Rust dump helper, and the two JSON/JSONL result files above.

**Done condition:** every one of the 71 `pub fn` in `descriptors.rs` has been exercised on a real 5,000-molecule corpus, every one of its ~190 individual values has a row in `validation/results/descriptor_census.json` with fixture/valid counts, RDKit-agreement stats (or an explicit "no RDKit oracle" reason), a dependency category, and a known/unexplained classification — met below.

---

## 1. Scope

CLAUDE.md scopes this census to `crates/chematic-chem/src/descriptors.rs` specifically — **not** the whole `chematic-chem` crate. That file has exactly **71 `pub fn`** (verified: `grep -c '^pub fn' descriptors.rs` → 71) and 4 `pub struct`s (`RingBundle`, `CarbonTypes`, `InformationContent`, `Bcut2D`). Many keys in the Python `Mol.descriptors()` dict (QED, SA score, Kappa/Chi/BertzCT/WienerIndex, the 47-value VSA families, EState, pKa, ADMET, xlogp3/esol/logd) live in **other files** in the same crate and are **not** covered here:

| Sibling file | What it owns | In `Mol.descriptors()`? |
|---|---|---|
| `topo_descriptors.rs` | Wiener/Zagreb/Randić/kappa1-3/chi0-4(v)/BertzCT/LabuteASA/VABC/gravitational_index/schultz_mti/gutman_mti/num_valence_electrons | yes |
| `vsa.rs` | SlogP_VSA(12)/SMR_VSA(10)/PEOE_VSA(14)/EState_VSA(11) — 47 values | yes |
| `qed.rs`, `sa_score.rs` | QED, SA score | yes |
| `estate.rs` | EState indices (sum/max/min) | yes |
| `xlogp3.rs`, `esol.rs`, `logd.rs` | alternative LogP/solubility/LogD | yes |
| `pka.rs` | pKa acid/base site prediction | yes |
| `admet.rs`, `alerts.rs`, `drug_score.rs` | BBB/Caco2/hERG/CYP3A4/AMES/PPB, PAINS/Brenk/REOS alerts, drug_score composite | yes |

This matches `docs/verification_coverage.md`'s own blank-cell item #2 ("170+ descriptor functions outside `bench5k.py`'s 19 ... Wiener/Zagreb/Randic/Balaban/kappa/chi/BCUT2D/MQN/VSA/2D-WHIM/2D-RDF/autocorrelation ... shipped and documented as '190+ descriptors' but only 19 have ever been checked against RDKit"). **This census closes the `descriptors.rs` slice of that gap; the sibling files above are a distinct, still-open gap** — flagged here, not started, no agent currently assigned. Recommend a follow-up census per sibling file using the same corpus and harness pattern.

Two `descriptors.rs` functions have a signature that doesn't fit "free function taking `&Molecule`" at all and are excluded from the count of independently-tested values as pure aliases/helpers (their outputs are 100% covered by other standalone functions already in the table below):
- `ring_bundle(mol)` → `RingBundle` struct, 13 fields, all of which duplicate standalone functions (`ring_count`, `num_aliphatic_rings`, `hba_count`, `rotatable_bond_count`, ...).
- `logp_and_mr(mol)` → `(f64, f64)` tuple, duplicates `logp_crippen`+`molar_refractivity`.
- `cns_mpo_from_parts(mol, logp, tpsa, mw, hbd, pka_b)` → takes 5 precomputed scalars, not just `&Molecule`; duplicates `cns_mpo_score`.

That leaves **68 independently-tested functions**, expanding to **196 individually-named descriptor values** in the census below (51 scalars/booleans/strings + BCUT2D×8 + MQN×42 + CarbonTypes×8 + 6 per-atom families + autocorr_2d×7 + usrcat×42 + moran×7 + geary×7 + InformationContent×5 + MDEC×10 + mmff94_charges×1) — consistent with CLAUDE.md's "190+ descriptor values."

---

## 2. Methodology

- **Oracle:** RDKit 2026.03.3 (verified: `.venv/bin/python -c "import rdkit; print(rdkit.__version__)"`), matching the other parallel diagnosis agents.
- **chematic build:** built from this branch's own worktree source (commit 659baca-based) into an **isolated venv** created specifically for this task (`maturin develop --release`), not the shared main-repo dev venv other agents use — this guarantees the numbers below reflect exactly this diagnosis's pinned commit regardless of what other agents commit to the shared main-repo checkout concurrently.
- **Corpus:** 5,000 unique canonical SMILES freshly downloaded from the ChEMBL REST API this session (`scripts/download_chembl_smiles.py --count 5000`), committed as `scripts/descriptor_census_corpus.smi` for reproducibility. This is a **different, independently-drawn** corpus from whatever produced `docs/validation.md`'s existing "4,999-molecule" numbers — deliberate, per the task's "re-measure everything yourself" instruction. 0 RDKit parse failures, 0 chematic parse failures, 5,000/5,000 evaluated.
- **Value sourcing:** every value is read via the individual Python getter/method that maps 1:1 to the `descriptors.rs` function under test (`mol.tpsa`, `mol.mqn()`, etc.) — **not** via the monolithic `Mol.descriptors()` dict. See §4 for why (a real hang was found and root-caused). Two values (`bcut2d`, `carbon_types`) have no individual Python getter; those are read from a small Rust dump helper instead (`crates/chematic-chem/examples/descriptor_census_unbound.rs`), same helper used for the 5 functions with no Python binding at all (§5).
- **RDKit reference for values with no 1:1 RDKit function:** for boolean drug-likeness filters (`ghose_passes`, `reos_passes`, `ro3_passes`, `lead_like_passes`, `veber_passes`, `egan_passes`, `pfizer_3_75_passes`, `lipinski_passes`) the exact threshold formula is copied from the Rust doc comment and applied to RDKit's own MW/LogP/TPSA/HBA/HBD/MR/ring-count/rotatable-bond values — this checks chematic's *arithmetic*, not RDKit's opinion of drug-likeness (RDKit has no such function). For `ring_system_count` and `num_ester_bonds` (no RDKit function either) a manual RDKit-based reference was implemented (union-find over `GetRingInfo().AtomRings()`; a `[#6](=O)[OX2H0][#6]` SMARTS). Where no reference is possible at all (Mordred-only families, the 2D-pseudo-USRCAT, RDKit's differently-shaped `CalcAUTOCORR2D`), the descriptor is measured for validity/reachability only and marked "no RDKit oracle" — never silently reported as 100%.
- **Reproduce:**
  ```bash
  # one-time: build the Rust dump helper's data
  cargo run -p chematic-chem --release --example descriptor_census_unbound \
      < scripts/descriptor_census_corpus.smi \
      > validation/results/descriptor_census_unbound.jsonl

  # main census (needs rdkit + chematic in the active venv)
  .venv/bin/python scripts/descriptor_census.py \
      --corpus scripts/descriptor_census_corpus.smi \
      --unbound validation/results/descriptor_census_unbound.jsonl \
      --json validation/results/descriptor_census.json
  ```

---

## 3. Headline (do not stop here — full table in §7)

196 individually-named values measured. Do **not** average these into one number — the spread is the finding:

- **139 have a real RDKit oracle**; of those, **111 (80%)** measure ≥99% exact-match on this corpus, confirming most of `docs/validation.md`'s existing claims and extending them to descriptors that were never previously checked (all the `num_*` element counts, `num_hydrogens`... well, see below, `num_ester_bonds`, `calc_mol_formula`, `carbon_types`, per-atom `formal_charge`/`implicit_hcount`).
- **57 have no meaningful RDKit oracle at all** (Mordred-only families, the 2D-pseudo-USRCAT, a mislabeled `ipc`, one dead/legacy `mmff94_charges`) — reported as such, not folded into a pass rate.
- A handful of the 139 oracled values are **badly wrong**, and three of those were root-caused down to the exact line of Rust source and a minimal reproducer (§6) — this is the actual point of the census.

---

## 4. Finding: a real hang, root-caused, and why this census doesn't call `Mol.descriptors()`

The first full-corpus run appeared to hang indefinitely. `ps` showed 100% CPU with no forward progress; `/usr/bin/sample <pid> 2` (macOS's non-destructive stack sampler — the compiled `.so` isn't stripped, so Rust symbol names resolve) showed **every sample across a 2-second window sitting in the identical frame**:

```
Mol.descriptors() → chematic_chem::drug_score::drug_score
                   → chematic_chem::alerts::pains_matches
                   → chematic_smarts::match_vf2::find_matches_with_rings_and_config
                   → match_vf2::match_recursive  (recursed ~15 levels deep, zero forward progress)
```

Bisected by printing each SMILES immediately before its `descriptors()` call: molecule index 28 in the corpus, a symmetric bis-isoquinolinium macrocycle —

```
C1=C\c2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCNc2cc[n+](c3ccccc23)Cc2ccc/1cc2
```

— triggers a catastrophic VF2 combinatorial blowup when PAINS SMARTS patterns are matched against it. This is real (not infinite — a second, longer wait confirmed it eventually returns, just after several minutes) and reproducible. **It is not in this census's scope**: `pains_matches`/`drug_score` live in `alerts.rs`/`drug_score.rs`, both outside `descriptors.rs`. The problem was purely that `Mol.descriptors()` unconditionally computes all ~195 dict keys — including `drug_score`, which nothing in this census's registry needed — on every molecule.

**Fix applied to the census script (not production code):** every value is now read via its individual getter (`mol.tpsa`, `mol.mqn()`, ...) instead of the combined dict. Verified directly: calling all ~50 individual getters this census needs, on the exact pathological molecule, takes 0.75s total (vs. several minutes via `.descriptors()`). `bcut2d`/`carbon_types` (no individual Python getter) are instead read from the Rust dump helper, which calls `chematic_chem::descriptors::{bcut2d, carbon_types}` directly and never touches `alerts.rs`.

**Candidate follow-up (not fixed here, real production defect):** VF2 substructure matching has a genuine performance cliff on symmetric macrocycles when matched against PAINS-style patterns — worth a dedicated look, but it's a `chematic-smarts`/`alerts.rs` issue, not a `descriptors.rs` one.

---

## 5. The 5 `descriptors.rs` functions with zero binding, and 1 correction found along the way

Grepping `chematic-py`, `chematic-mcp`, `chematic-wasm` for these names returns **zero references**:

- `moran_autocorr(mol)` — Moran's I, 7 lags
- `geary_autocorr(mol)` — Geary's C, 7 lags
- `information_content(mol)` — IC/TIC/SIC/BIC/CIC, 5 values
- `mde_carbon(mol)` — Molecular Distance Edge, 10 values
- `mmff94_charges(mol)` — see correction below

These are dumped via a small new Rust example (`crates/chematic-chem/examples/descriptor_census_unbound.rs`, not `src/`) that calls them directly and writes one JSON object per line. None have an RDKit equivalent (Mordred-only families), so they're measured for validity/reachability only.

**A real panic found while building the dump**, with a root cause: `moran_autocorr`/`geary_autocorr` panic on `'[2H]C([2H])([2H])NC=O'` with `index out of bounds: the len is 4 but the index is 4`. Root cause: `topological_distance_matrix()` returns a matrix indexed by **heavy-atom-compacted position** (skipping explicit H/D atoms), but `moran_autocorr`/`geary_autocorr` loop over `0..mol.atom_count()` (the **full**, uncompacted atom count, which includes the 3 explicit deuteriums) and use that loop index as a raw `AtomIdx` into the full atom list. `autocorr_2d`/`ipc` use the exact same pattern but are defensively bounded by `.take(n)` on the already-short matrix, so they don't crash — **they silently look up the wrong atom's property instead**, on any molecule with explicit H/D/T atoms in the graph. Confirmed root cause, not fixed here (production change). Handled in the dump helper with `panic::catch_unwind` per function so one crash doesn't blank out the other 4 values for the same molecule (1/5,000 molecules affected on this corpus).

**A correction made mid-diagnosis** (documented per the "verify, don't assume" discipline this project's memory repeatedly emphasizes): an earlier pass of this same census assumed `descriptors.rs::mmff94_charges` was a distinct, inferior "electronegativity-weighted + formal charge" formula, based on reading only its stale doc comment. Reading the actual function body shows it is a 1-line pass-through to `mmff94_bci::mmff94_charges_bci()` — **the exact same function** Python's `mol.mmff94_charges()` calls. Values are byte-identical to production output. It's still unreachable *as this specific symbol* (dead-wrapper, not dead/wrong-formula), and the doc comment is stale — a minor documentation-drift finding, not a functional one.

---

## 6. Other concrete, source-verified findings (not fixed — diagnosis only)

Each of these was root-caused by reading the actual Rust source and, where possible, confirmed with a minimal reproducer — not just observed as a number.

| Descriptor | Corpus result | Root cause |
|---|---|---|
| `num_unspecified_stereocenters` | **6.12% exact**, MAE 6.17, max 48 | `num_unspecified_stereocenters()` (descriptors.rs ~line 1808) never checks substituent distinctness — unlike `num_stereocenters()` (99.76% exact), which does a real CIP-distinctness pass. It flags *any* sp3 carbon with degree+implicitH==4, no attached double/triple bond, and no explicit chirality tag — including ordinary `-CH2-`/`-CH3` groups, which can never be real stereocenters (two identical H substituents). Worst fixture is a plain polypeptide where every real stereocenter already has `@`/`@@` (RDKit correctly says 0 unspecified) but chematic says 48 — almost certainly counting backbone `-CH2-`/`-CH<` carbons. Corpus-wide, not a macrocycle edge case. |
| `num_hydrogens` | **64.88% exact**, MAE 0.90, max 22 | Confirmed with a minimal 4-atom reproducer: `num_hydrogens()` (~line 2852) sums `atom.hydrogen_count` (explicit, from bracket notation) **plus** `implicit_hcount()` for every atom. For `'C[C@H](N)O'`, `implicit_hcount_per_atom()` is `[3,1,2,1]` (sums to 7, exactly matching RDKit) but `num_hydrogens()` returns 8 — `implicit_hcount()` doesn't know the `[C@H]` atom's H is already counted, so it's counted twice. Double-counts by exactly the number of bracket-H stereocenters (`[C@H]`, `[C@@H]`, ...) in the molecule — extremely common in real SMILES. |
| `molecular_weight` | 99.98% exact overall; worst fixture off by 3.02 Da | `molecular_weight()` (~line 186) computes mass via `avg_mass(atom.element)` only and never reads `atom.isotope` — an explicit deuterium/tritium atom is weighed as ordinary natural-abundance H. `exact_mass()` (same file) correctly reads `atom.isotope`; this molecule is not `exact_mass`'s own worst case, confirming the two diverge specifically on isotope handling. Same trigger molecule as the moran/geary panic above (unrelated mechanism, coincidental shared fixture — `[2H]C([2H])([2H])NC=O` is a good adversarial-fixture candidate generally). |
| `ipc` | **0.00% exact**, MAE 5.4e47, max 1.7e51 | **Name collision, not a numeric bug.** `ipc()` computes `Σ deg_i·deg_j / d(i,j)²` over atom pairs (docstring: "Information Path Count"). RDKit's `GraphDescriptors.Ipc` is Bonchev-Trinajstić total information content on distance-degeneracy classes — a completely unrelated formula that happens to share the abbreviation. Comparing them produces astronomical MAE by construction; this says nothing about correctness. |
| `calc_mol_formula` | 99.66% exact | chematic's Hill-notation formula omits the ionic-charge suffix RDKit appends (`C42H46N4+2` vs `C42H46N4` for a dication) — atom/element counts are correct, only the charge annotation is missing. |
| `ring_count`, `num_aliphatic_rings`, `num_aliphatic_heterocycles` | 99.16% / 99.86% / 99.88% exact, max error up to 13 | Same mechanism as the already-documented `[R1]`/`[R2]` SMARTS divergence in `docs/rdkit_compat.md`'s Known divergence classes table: genuine SSSR-*basis*-cardinality disagreement on large fused/bridged macrocycles (~1.7% of such molecules per that doc), not a predicate bug. Confirmed again here on independently-drawn molecules — same class, not a new defect. |
| `hall_kier_alpha`, `bcut2d_*` (8), `MQN*` (42) | 4.84%, 0.16–2.10%, mostly single-digit-to-low-double-digit % exact | Matches `docs/rdkit-comparison.md`'s existing disclosure: "Kappa/HallKierAlpha/BertzCT/BalabanJ/BCUT2D/VSA/MQN/SAScore were found to diverge substantially once measured at corpus scale." Re-measured independently here and confirmed for the 3 of those 8 named descriptors that actually live in `descriptors.rs` (the other 5 — Kappa/BertzCT/VSA/SAScore — are in sibling files, out of scope, see §1). |
| `balaban_j` | **100.00% exact**, MAE 0.0 | Named in that same `docs/rdkit-comparison.md` sentence above alongside HallKierAlpha/BCUT2D/MQN as "diverging substantially" — but it does **not** diverge in this round's measurement. The doc's blanket phrasing over-generalizes; this census gives the precise per-descriptor breakdown that sentence was missing. |
| `hybridization_per_atom` | 90.72% exact (135,946 atom-instances) | Minor divergence: chematic's hybridization is a 3-way heuristic (has-triple → sp, has-double-or-aromatic → sp2, else sp3); RDKit's is a fuller model accounting for lone pairs/resonance on N/O. Not investigated further (minor tier, not a crash or gross miscount). |
| `logp_crippen_per_atom`, `mr_per_atom` | 43.05% exact, but NOT a real divergence | Attribution-convention difference: RDKit's `_CalcCrippenContribs()` returns contributions for heavy atoms only — each attached implicit H's own Crippen atom-type contribution is computed separately inside `MolLogP`/`MolMR` and is *not* in the per-atom array (confirmed: summing 3 rows for `'CCO'` gives −0.3487, but `Crippen.MolLogP('CCO')` is −0.0014). chematic instead folds each atom's attached-H contribution into that heavy atom's own value. `formal_charge_per_atom` and `implicit_hcount_per_atom` both measure 100% exact on this same corpus, ruling out atom-order misalignment. The molecule-level aggregates (`logp_crippen`, `molar_refractivity`, both ≥99.98%) are the correct level to trust this at — this row is not a candidate follow-up. |
| `usrcat` (36 shape values) | no RDKit oracle possible | Not a divergence measurement — a scope note: this is a 2D-topology-only pseudo-USRCAT (a distance-matrix average scaled by `1.0 + slot/12.0`, not real USR moments); RDKit's `GetUSRCAT` needs a real 3D conformer and computes genuine shape descriptors, so the two are not comparable even in principle. The scaling formula looks like placeholder/synthetic code — worth a closer look outside this diagnosis, not investigated further here. |

---

## 7. Full per-descriptor table (unabridged — 196 rows, 12 groups)

Legend: **Oracle** = has a real RDKit reference (see §2 for what counts). **Exact%** = exact-match rate among valid (both-sides-computed, non-NaN/Inf) fixtures — `n/a` means no oracle exists, never inferred. Worst fixture is truncated to 60 chars for table width; full SMILES are in `validation/results/descriptor_census.json`.

### 7.1 Scalars, booleans, and the molecular formula (53 values)

| Descriptor | Dependency | Oracle | Fixtures | Valid | Exact% | MAE | Median AE | p95 AE | Max AE | Worst fixture |
|---|---|---|---|---|---|---|---|---|---|---|
| `molecular_weight` | independent | yes | 5000 | 5000 | 99.98% | 0.00112 | 0 | 0.002 | 3.018 | `[2H]C([2H])([2H])NC=O` (ch=59.07, rd=62.09) — see §6 |
| `exact_mass` | independent | yes | 5000 | 5000 | 99.96% | 0.001366 | 0.000144 | 0.0003 | 6.009 | large drug-like molecule, ch=742.3 rd=748.4 |
| `heavy_atom_count` | independent | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `hbd_count` | depends-on-aromaticity | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `hba_count` | depends-on-aromaticity | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `rotatable_bond_count` | depends-on-ring-perception | yes | 5000 | 5000 | 99.96% | 0.0004 | 0 | 0 | 1 | `O=C1CCN(N2C(=O)c3ccccc3C2=O)C(=O)N1` (ch=0, rd=1) |
| `tpsa` | depends-on-aromaticity | yes | 5000 | 5000 | 99.76% | 0.0294 | 0 | 0 | 25.3 | `[S-]c1nc2ccccc2c2cccc[n+]12` (ch=42.29, rd=16.99) |
| `logp_crippen` | independent | yes | 5000 | 5000 | 99.98% | 0.000115 | 0 | 0 | 0.5747 | fused polycyclic, ch=4.582 rd=4.007 |
| `lipinski_passes` | independent | yes | 5000 | 5000 | 100.00% | - | - | - | - | — |
| `fsp3` | depends-on-aromaticity | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `aromatic_ring_count` | depends-on-aromaticity | yes | 5000 | 5000 | 99.96% | 0.0004 | 0 | 0 | 1 | large fused polycyclic (ch=18, rd=19) |
| `formal_charge_sum` | independent | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `molar_refractivity` | independent | yes | 5000 | 5000 | 99.98% | 2e-05 | 0 | 0 | 0.102 | ch=105.7 rd=105.8 |
| `num_heteroatoms` | independent | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `ring_count` | depends-on-ring-perception | yes | 5000 | 5000 | 99.16% | 0.0098 | 0 | 0 | 3 | bis-pyridinium macrocycle (ch=7, rd=10) — see §6 |
| `ring_system_count` | depends-on-ring-perception | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `hba_count_lipinski` | depends-on-aromaticity | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `fraction_rotatable_bonds` | depends-on-ring-perception | yes | 5000 | 5000 | 99.96% | 2.6e-05 | 0 | 0 | 0.077 | `CNC(=O)ON(C(C)=O)C(=O)NC` (ch=0.077, rd=0) |
| `num_aliphatic_rings` | depends-on-aromaticity | yes | 5000 | 5000 | 99.86% | 0.0108 | 0 | 0 | 13 | bis-isoquinolinium macrocycle (ch=16, rd=3) — see §6 |
| `num_saturated_rings` | depends-on-aromaticity | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `num_aromatic_heterocycles` | depends-on-aromaticity | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `num_aliphatic_heterocycles` | depends-on-aromaticity | yes | 5000 | 5000 | 99.88% | 0.0106 | 0 | 0 | 13 | same macrocycle as above — see §6 |
| `num_saturated_heterocycles` | depends-on-aromaticity | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `num_spiro_atoms` | depends-on-ring-perception | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `num_bridgehead_atoms` | depends-on-ring-perception | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `num_stereocenters` (legacy oracle) | depends-on-CIP-or-stereo | yes | 5000 | 5000 | 99.76% | 0.003 | 0 | 0 | 2 | `CCCCCCCCC1CC1CCCCCCCC(=O)O` (ch=0, rd=2) |
| `num_stereocenters` (new-CIP oracle) | depends-on-CIP-or-stereo | yes | 5000 | 5000 | 98.96% | 0.0208 | 0 | 0 | 5 | `C=C1CC2(OC1=O)C1CC3CC(C1)CC2C3` (ch=0, rd=5) |
| `num_unspecified_stereocenters` | depends-on-CIP-or-stereo | yes | 5000 | 5000 | **6.12%** | 6.173 | 5 | 17 | 48 | polypeptide (ch=48, rd=0) — **see §6, root-caused** |
| `veber_passes` | depends-on-ring-perception | yes | 5000 | 5000 | 99.98% | - | - | - | - | `NC(N)=NCCCC(N)[PH](=O)O` (ch=True, rd=False) |
| `egan_passes` | depends-on-aromaticity | yes | 5000 | 5000 | 99.98% | - | - | - | - | `CNC(=N)NCCCC(N)[PH](=O)O` (ch=True, rd=False) |
| `reos_passes` | depends-on-ring-perception | yes | 5000 | 5000 | 100.00% | - | - | - | - | — |
| `ghose_passes` | independent | yes | 5000 | 5000 | 100.00% | - | - | - | - | — |
| `ro3_passes` | depends-on-ring-perception | yes | 5000 | 5000 | 100.00% | - | - | - | - | — |
| `lead_like_passes` | depends-on-ring-perception | yes | 5000 | 5000 | 99.94% | - | - | - | - | `CCO[C@H]1O[C@@H]2O[C@@]3(C)CCC4...` (ch=True, rd=False) |
| `pfizer_3_75_passes` | depends-on-aromaticity | yes | 5000 | 5000 | 100.00% | - | - | - | - | — |
| `cns_mpo_score` | out-of-scope (pka.rs+logd.rs) | no | 5000 | 5000 | n/a | - | - | - | - | no oracle attempted this round |
| `mcf_passes` | out-of-scope (alerts.rs) | no | 5000 | 5000 | n/a | - | - | - | - | no oracle attempted this round |
| `num_carbons` | independent | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `num_nitrogens` | independent | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `num_oxygens` | independent | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `num_fluorines` | independent | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `num_chlorines` | independent | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `num_bromines` | independent | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `num_iodines` | independent | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `num_sulfurs` | independent | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `num_phosphorus` | independent | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `num_hydrogens` | independent | yes | 5000 | 5000 | **64.88%** | 0.9048 | 0 | 4 | 22 | polypeptide (ch=196, rd=174) — **see §6, root-caused** |
| `num_amide_bonds` | independent | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `num_ester_bonds` | independent | yes | 5000 | 5000 | 100.00% | 0 | 0 | 0 | 0 | — |
| `calc_mol_formula` | independent | yes | 5000 | 5000 | 99.66% | n/a (string) | n/a | n/a | n/a | dication (ch=C42H46N4, rd=C42H46N4+2) — see §6 |
| `balaban_j` | depends-on-ring-perception | yes | 5000 | 5000 | **100.00%** | 0 | 0 | 0 | 0 | — see §6 (contradicts blanket "diverges" claim) |
| `ipc` | depends-on-ring-perception | yes | 5000 | 5000 | **0.00%** | 5.4e+47 | 4.5e+05 | 2.1e+10 | 1.7e+51 | — **see §6, name collision, not a bug** |
| `hall_kier_alpha` | depends-on-ring-perception | yes | 5000 | 5000 | **4.84%** | 0.1829 | 0.1325 | 0.5044 | 2.538 | — see §6, confirms existing disclosure |

### 7.2 BCUT2D (8 values) — confirms `docs/rdkit-comparison.md`'s existing disclosure

| Descriptor | Oracle | Exact% | MAE | Max AE | Worst fixture |
|---|---|---|---|---|---|
| `bcut2d_mwhi` | yes | 0.52% | 0.357 | 3.41 | `COS(C)(=O)=S` |
| `bcut2d_mwlo` | yes | 0.24% | 0.411 | 0.91 | `CN(C)CC#CCN1OCCC1=O` |
| `bcut2d_chghi` | yes | 0.70% | 0.586 | 3.41 | `COS(C)(=O)=S` |
| `bcut2d_chglo` | yes | 1.46% | 0.555 | 3.37 | `COS(C)(=O)=S` |
| `bcut2d_logphi` | yes | 0.60% | 0.661 | 3.43 | `COS(C)(=O)=S` |
| `bcut2d_logplo` | yes | 2.10% | 0.500 | 3.37 | `COS(C)(=O)=S` |
| `bcut2d_mrhi` | yes | 0.16% | 0.825 | 3.69 | `COS(C)(=O)=S` |
| `bcut2d_mrlo` | yes | 0.28% | 0.527 | 1.86 | sulfonamide, ch=-1.27 rd=0.60 |

### 7.3 MQN — Molecular Quantum Numbers (42 values) — confirms existing disclosure, large spread

Full table (see `validation/results/descriptor_census.json` for exact worst-fixture SMILES per value — omitted here for width):

| # | Exact% | MAE | Max AE | | # | Exact% | MAE | Max AE |
|---|---|---|---|---|---|---|---|---|
| MQN1 | 100.00% | 0 | 0 | | MQN22 | 7.56% | 2.47 | 40 |
| MQN2 | 15.52% | 2.57 | 35 | | MQN23 | 18.86% | 1.53 | 32 |
| MQN3 | 10.74% | 3.36 | 31 | | MQN24 | 0.10% | 21.25 | 132 |
| MQN4 | 85.96% | 0.28 | 13 | | MQN25 | 37.70% | 1.46 | 41 |
| MQN5 | 98.74% | 0.01 | 4 | | MQN26 | 3.50% | 4.52 | 52 |
| MQN6 | 77.40% | 0.28 | 3 | | MQN27 | 0.00% | 123.2 | 129 |
| MQN7 | 77.40% | 0.28 | 3 | | MQN28 | 34.82% | 1.73 | 54 |
| MQN8 | 41.42% | 1.17 | 33 | | MQN29 | 16.58% | 0.94 | 5 |
| MQN9 | 31.46% | 1.51 | 12 | | MQN30 | 2.48% | 6.81 | 33 |
| MQN10 | 11.78% | 3.07 | 31 | | MQN31 | 1.34% | 5.54 | 58 |
| MQN11 | 0.04% | 14.56 | 131 | | MQN32 | 5.04% | 5.72 | 81 |
| MQN12 | 0.00% | 25.10 | 136 | | MQN33 | 7.86% | 11.80 | 48 |
| MQN13 | 0.94% | 10.15 | 131 | | MQN34 | 16.42% | 1.88 | 35 |
| MQN14 | 1.56% | 10.98 | 66 | | MQN35 | 3.06% | 4.11 | 33 |
| MQN15 | 1.36% | 3.09 | 32 | | MQN36 | 17.50% | 1.70 | 20 |
| MQN16 | 1.08% | 8.81 | 58 | | MQN37 | 54.10% | 0.63 | 9 |
| MQN17 | 4.16% | 5.73 | 27 | | MQN38 | 0.00% | 27.19 | 165 |
| MQN18 | 0.00% | 1.01 | 2 | | MQN39 | 5.70% | 7.60 | 69 |
| MQN19 | 12.18% | 3.62 | 78 | | MQN40 | 39.20% | 1.07 | 86 |
| MQN20 | 3.22% | 7.71 | 93 | | MQN41 | 99.22% | 0.05 | 24 |
| MQN21 | 4.02% | 5.01 | 64 | | MQN42 | 41.46% | 0.88 | 30 |

Notable: MQN1 (C count) and MQN41 (spiro count) are near-perfect; MQN2/3/8-27/etc (bond-type and degree-statistic bins) diverge heavily — consistent with `docs/rdkit-comparison.md`'s "found to diverge substantially" note, now precisely quantified per-index rather than as one blanket statement.

### 7.4 CarbonTypes — Mordred-style hybridization×degree (8 values)

| Descriptor | Exact% | MAE | Max AE |
|---|---|---|---|
| `c1sp1` | 100.00% | 0 | 0 |
| `c2sp1` | 99.54% | 0.0046 | 1 |
| `c1sp2` | 100.00% | 0 | 0 |
| `c2sp2` | 99.54% | 0.0046 | 1 |
| `c3sp2` | 100.00% | 0 | 0 |
| `c1sp3` | 100.00% | 0 | 0 |
| `c2sp3` | 100.00% | 0 | 0 |
| `c3sp3` | 100.00% | 0 | 0 |

### 7.5 Per-atom families (6 values, atom-instance-level fixtures — 135,946–135,949 atom pairs across the 5,000-molecule corpus, not molecule-level)

| Descriptor | Dependency | Fixtures | Exact% | MAE | Max AE | Note |
|---|---|---|---|---|---|---|
| `hybridization_per_atom` | independent | 135,946 | 90.72% | 0.093 | 1 | minor — simplified heuristic vs RDKit's fuller model, see §6 |
| `formal_charge_per_atom` | independent | 135,949 | 100.00% | 0 | 0 | — |
| `implicit_hcount_per_atom` | independent | 135,949 | 100.00% | 0 | 0 | — |
| `tpsa_per_atom` | depends-on-aromaticity | 135,949 | 98.29% | 0.392 | 40.1 | one nitro-group outlier, `N[C@@H](CCC(=O)Nc1ccc([N+](=O)[O-])cc1)C(=O)O` |
| `logp_crippen_per_atom` | independent | 135,949 | 43.05% | 0.115 | 0.575 | **not a real divergence, see §6** |
| `mr_per_atom` | independent | 135,949 | 43.05% | 0.920 | 3.171 | **not a real divergence, see §6** (same mechanism as above) |

### 7.6 `autocorr_2d` — 7 lags, no RDKit oracle (different definition, see §6)

All 7 lags: 5,000/5,000 valid, no oracle. Shares the heavy-atom-index lookup pattern flagged in §5 for `moran_autocorr`/`geary_autocorr` — doesn't crash (defensively `.take(n)`-bounded) but is at risk of silently wrong values on explicit-H/D/T molecules; not separately quantified this round.

### 7.7 `usrcat` — 42 values, no RDKit oracle (2D pseudo-shape descriptor, see §6)

All 42: 5,000/5,000 valid, no oracle possible even in principle.

### 7.8 `moran_autocorr` / `geary_autocorr` — 7+7 lags, unreachable via any binding (see §5)

All 14: 5,000 fixtures, **4,999 valid** (1 panic, root-caused in §5), no RDKit oracle (Mordred-only family).

### 7.9 `information_content` (IC/TIC/SIC/BIC/CIC, 5 values) and `mde_carbon` (MDEC, 10 values) — unreachable via any binding, no RDKit oracle

All 15: 5,000/5,000 valid, no oracle (Mordred-only families).

### 7.10 `mmff94_charges` (plain) — unreachable via any binding, byte-identical to production (see §5)

5,000/5,000 valid, no oracle attempted (not a distinct formula from the production `mmff94_charges_bci` — see §5's correction).

---

## 8. `depends-on-aromaticity` — explicit list for the aromaticity-diagnosis agent

Per the coordination instructions, this census does **not** wait on or investigate the parallel aromaticity-diagnosis agent's work — it only flags what should be re-measured once that agent's fixes land in production. These 14 values (of 196) are tagged `depends-on-aromaticity` in `validation/results/descriptor_census.json`:

- `hbd_count` (100.00% today)
- `hba_count` (100.00%)
- `tpsa` (99.76%)
- `fsp3` (100.00%)
- `aromatic_ring_count` (99.96%)
- `hba_count_lipinski` (100.00%)
- `num_aliphatic_rings` (99.86%)
- `num_saturated_rings` (100.00%)
- `num_aromatic_heterocycles` (100.00%)
- `num_aliphatic_heterocycles` (99.88%)
- `num_saturated_heterocycles` (100.00%)
- `egan_passes` (99.98%)
- `pfizer_3_75_passes` (100.00%)
- `tpsa_per_atom` (98.29%)

Most are already near-100% here — the ones worth re-checking first after an aromaticity fix are `tpsa`/`tpsa_per_atom` (the two with the most room, and the largest single outliers: 25.3 Å² and 40.1 Å² respectively) and `num_aliphatic_rings`/`num_aliphatic_heterocycles` (whose residual is dominated by the SSSR-cardinality macrocycle class in §6, not obviously aromaticity-related, but worth re-confirming that attribution once the aromaticity model changes underfoot).

`depends-on-ring-perception` (55 values, mostly MQN sub-indices) and `depends-on-CIP-or-stereo` (3 values: `num_stereocenters` ×2 oracles + `num_unspecified_stereocenters`) are also tagged in the JSON for the same reason, in case a ring-perception or CIP/stereo specialist agent wants the same kind of list — not reproduced in full here since the task only asked for the aromaticity list explicitly.

---

## 9. What this census does NOT claim

- It does not re-validate the 19 descriptors `bench5k.py` already covers as a duplicate effort for its own sake — where this census's numbers land within noise of `docs/validation.md`'s claims (HBA/HBD/TPSA/LogP/MR/Fsp3/ring counts/rotatable bonds/spiro/bridgehead/amide bonds/aromatic-heterocycles), that's treated as **independent confirmation** on a fresh corpus, not "already known, skip." Where it disagrees even slightly (e.g. `tpsa` 99.76% here vs. `docs/validation.md`'s claimed 100% on a different corpus), the discrepancy is a genuinely different, independently-drawn 5,000-molecule sample, not a contradiction — both are within the ±0.1 Å² tolerance's expected noise band; not investigated further as its own issue.
- It does not fix anything. Every "candidate follow-up" and "FLAG FOR ... SPECIALIST" note above is exactly that — a diagnosis handoff, not a patch.
- It does not measure `cns_mpo_score` or `mcf_passes` against a real oracle — both require `pka.rs`/`logd.rs`/`alerts.rs`, out of this file-scoped census; flagged, not measured.
