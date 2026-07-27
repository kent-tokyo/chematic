# chematic-mcp tool inventory

**This is a manually maintained snapshot, not a generated artifact.** It is
not the single source of truth for `chematic-mcp`'s tool catalog — that
remains `crates/chematic-mcp/src/tools.rs` (`list_tools()` / `call_tool()`)
until a PR 2 registry refactor can generate this document (or replace it)
from a single typed definition. If this file and `tools.rs` ever disagree,
`tools.rs` is correct.

**Protocol note (2026-07-28):** the 20 tools and their behavior below are
unchanged by the 2026-07-28 protocol work — every tool now additionally
carries an `outputSchema` (JSON Schema 2020-12) alongside its existing
`inputSchema`, and its result is exposed as `structuredContent` in the
modern protocol era, but the underlying chemistry, determinism,
completeness, and network-access classifications on this page are exactly
as they were. See `docs/mcp/2026-07-28-implementation-rfc.md` for the
protocol-layer design and `crates/chematic-mcp/README.md` for the
legacy/modern wire examples.

## Classification axes

Adapted from the remote-readiness design discussion (not yet implemented as
Rust types — this document uses them as prose labels only):

- **Determinism**: `DeterministicByDesign` (no intentional random sampling is
  used) / `Stochastic` (uses randomization) / `ExternalDependent` (result
  depends on a third-party service's state at call time). `DeterministicByDesign`
  is **not**, by itself, a guarantee of byte-identical results across
  processes, platforms, hash seeds, or timeout-limited executions — it only
  means no `rand`/RNG call is in the path. `HashMap`/`HashSet` iteration
  order, floating-point results across platforms, and timeout-limited search
  cutoffs can all still vary without any randomization being involved.
- **Result method**: `Algorithmic` (exact computation over the input
  structure) / `RuleBasedEstimate` (empirical/heuristic model — screening
  estimate, not a measurement) / `HeuristicSearch` (explores a solution
  space without a completeness guarantee) / `ExternalLookup` (retrieves a
  value from a third-party database rather than computing one).
- **Completeness**: `Complete` (always runs to a definitive result) /
  `TimeoutLimited` (may return a truncated/non-optimal result under a time
  budget) / `BestEffort` (may fail outright — network dependency, or a
  heuristic that doesn't guarantee it found anything).

No tool in this inventory uses randomization — verified by grepping
`chematic-3d`, `chematic-chem`, and `chematic-smarts` for `rand::`/
`thread_rng`/`StdRng`/`SmallRng` (no hits). That rules out *stochastic*
sampling as a source of variation for `generate_3d` and `find_mcs`; it does
not by itself establish byte-identical output across runs, machines, or
timeout-affected executions — see each tool's own entry below for what is
and isn't established.

## Tools

### `parse_smiles`
- Description: Parse a SMILES string → atom count, bond count, molecular weight.
- Local/network: local
- Determinism: DeterministicByDesign — Result method: Algorithmic — Completeness: Complete
- High cost: No
- Explicit limits: none
- Caveats: no SMILES length limit; unbounded input size.

### `calc_properties`
- Description: MW, exact mass, Crippen LogP, TPSA, HBD, HBA, rotatable bond count, heavy atom count, QED.
- Local/network: local
- Determinism: DeterministicByDesign — Result method: **mixed** (MW/exact mass/TPSA/HBD/HBA/rotatable bonds/heavy atom count are exact structural computations; LogP (Crippen atom-contribution) and QED (weighted composite) are empirical `RuleBasedEstimate`) — Completeness: Complete
- High cost: No
- Explicit limits: none
- Caveats: bundles exact counts and empirical estimates in one JSON object without a field-level distinction — a caller cannot currently tell which of the 9 fields are exact vs. modeled without reading this document.

### `ecfp4`
- Description: ECFP4 (Morgan radius-2) circular fingerprint, 2048-bit hex + popcount.
- Local/network: local
- Determinism: DeterministicByDesign (within chematic's own implementation) — Result method: Algorithmic — Completeness: Complete
- High cost: No
- Explicit limits: none
- Caveats: deterministic and reproducible within chematic, but **not** guaranteed RDKit bit-exact by default — chematic's ECFP4 invariant includes aromaticity, RDKit's default does not (see `ecfp4_agreement_methodology` in project history).

### `tanimoto`
- Description: Tanimoto (Jaccard) similarity between two molecules' ECFP4 fingerprints.
- Local/network: local
- Determinism: DeterministicByDesign — Result method: Algorithmic — Completeness: Complete
- High cost: No
- Explicit limits: none
- Caveats: inherits `ecfp4`'s non-RDKit-bit-exact caveat.

### `smarts_match`
- Description: SMARTS substructure search — match count + atom index maps.
- Local/network: local
- Determinism: DeterministicByDesign — Result method: Algorithmic (VF2-based exact isomorphism) — Completeness: Complete (exhaustive match enumeration)
- High cost: **Potentially** — subgraph isomorphism is worst-case exponential; no timeout guard exists for this tool (unlike `find_mcs`).
- Explicit limits: none
- Caveats: an adversarial SMARTS/molecule pair could run long with no internal cutoff.

### `canonical_smiles`
- Description: Canonical SMILES representation.
- Local/network: local
- Determinism: DeterministicByDesign — Result method: Algorithmic — Completeness: Complete
- High cost: No
- Explicit limits: none
- Caveats: deterministic and reproducible within chematic, but the canonical **string** is chematic's own canonical form — not necessarily identical to RDKit's canonical SMILES for the same molecule (different canonicalization algorithms, both valid).

### `find_mcs`
- Description: Maximum common substructure across 2–20 molecules.
- Local/network: local
- Determinism: Deterministic search procedure, **EnvironmentSensitive outcome under timeout** — the search algorithm itself has no randomization, but the 5000ms wall-clock budget means the same input can hit the cutoff at a different point in the search on a slower machine or a loaded runner, changing the returned result. Not `DeterministicByDesign` without qualification. — Result method: Algorithmic — Completeness: **TimeoutLimited** (internal 5000ms budget; a run that hits the budget returns the best substructure found so far, not necessarily the true maximum)
- High cost: Yes — MCS is NP-hard in general; this is one of two tools with an internal timeout for exactly that reason.
- Explicit limits: 2–20 molecules, ≤200 atoms per molecule, 5000ms internal search timeout.
- Caveats: do not describe this tool's result as "the" MCS without qualification — under the timeout it is a best-effort lower bound, and the exact result can vary run-to-run on borderline-slow inputs depending on machine speed/load.

### `generate_3d`
- Description: 3D coordinates via rule-based placement + DREIDING force-field minimization (XYZ).
- Local/network: local
- Determinism: DeterministicByDesign (verified: no RNG anywhere in `chematic-3d`) — Result method: HeuristicSearch (DREIDING minimization is a local energy optimization from a heuristic starting geometry, not a guaranteed global minimum) — Completeness: **BestEffort**
- High cost: **Yes, and unguarded** — no atom-count limit exists for this tool at all, unlike `find_mcs`/`retrosynthesis`.
- Explicit limits: **none** (known gap — see the audit's resource-risk findings)
- Caveats: single deterministic conformer, not an ensemble or global-minimum search; DREIDING is an approximate universal force field, not a high-accuracy method. `BestEffort`, not `Complete`, because: coordinate generation or minimization can fail to produce a usable geometry for some inputs; local minimization is not guaranteed to converge to (or even near) the global minimum; convergence has not been proven or tested across the full space of possible inputs; and there is no atom-count or time budget guarding the minimization step, so behavior on large/pathological inputs is unverified.

### `pains_check`
- Description: PAINS (Pan-Assay Interference Compounds) structural alerts.
- Local/network: local
- Determinism: DeterministicByDesign (SMARTS pattern matching against a fixed alert list is deterministic execution) — Result method: RuleBasedEstimate (the alert list itself is a curated empirical correlation, Baell & Holloway 2010, not a physical law) — Completeness: Complete
- High cost: No
- Explicit limits: none
- Caveats: a match indicates historical association with assay interference in specific HTS contexts, not a guarantee of interference in any given assay; absence of alerts is not a guarantee of clean behavior.

### `brenk_check`
- Description: Brenk structural alerts (toxicity / metabolic instability / reactivity).
- Local/network: local
- Determinism: DeterministicByDesign — Result method: RuleBasedEstimate (Brenk et al. 2008 empirical alert list) — Completeness: Complete
- High cost: No
- Explicit limits: none
- Caveats: same as `pains_check` — empirical alert list, not a predictive toxicology model.

### `sa_score`
- Description: Synthetic accessibility score (Ertl & Schuffenhauer 2009), 1 (easy) to 10 (hard).
- Local/network: local
- Determinism: DeterministicByDesign — Result method: RuleBasedEstimate (fragment-contribution empirical scoring) — Completeness: Complete
- High cost: No
- Explicit limits: none
- Caveats: an empirical estimate correlated with fragment frequency in known compounds, not a retrosynthetic feasibility proof.

### `admet_profile`
- Description: BBB, Caco-2, hERG, CYP3A4, AMES, PPB, hepatic clearance.
- Local/network: local
- Determinism: DeterministicByDesign — Result method: RuleBasedEstimate — Completeness: Complete
- High cost: No
- Explicit limits: none
- Caveats: **every field is a rule-based/QSAR-style screening estimate, not an experimental measurement or clinical/regulatory conclusion.** This is the clearest case in the catalog where the current tool description ("Full ADMET profile") risks reading as more authoritative than the underlying method supports.

### `boiled_egg`
- Description: BOILED-Egg (Daina & Zoete 2016) — GI absorption + BBB zone prediction from LogP/TPSA thresholds.
- Local/network: local
- Determinism: DeterministicByDesign — Result method: RuleBasedEstimate (threshold-based empirical classification) — Completeness: Complete
- High cost: No
- Explicit limits: none
- Caveats: threshold classification, not a mechanistic permeability simulation.

### `lipinski_check`
- Description: Lipinski's Rule of Five with per-rule breakdown.
- Local/network: local
- Determinism: DeterministicByDesign — Result method: RuleBasedEstimate — Completeness: Complete
- High cost: No
- Explicit limits: none
- Caveats: a heuristic screening rule with well-known exceptions (many approved oral drugs violate Ro5); pass/fail is not a bioavailability guarantee in either direction.

### `name_to_smiles`
- Description: Chemical name → isomeric SMILES via the PubChem REST API.
- Local/network: **network** — the only tool in the catalog that makes an external call.
- Determinism: **ExternalDependent** (result depends on PubChem's database content and availability at call time, not solely on chematic's own code) — Result method: **ExternalLookup** (retrieves a value from PubChem's database; this is not chemical-algorithm-driven structure generation — chematic computes nothing here beyond encoding the request and parsing the response) — Completeness: **BestEffort** (fails outright on network/PubChem issues; no local fallback)
- High cost: No (bounded by a 10-second HTTP timeout)
- Explicit limits: name ≤500 characters, 10s HTTP timeout.
- Caveats: input name is sent, URL-encoded, to `pubchem.ncbi.nlm.nih.gov`; see `crates/chematic-mcp/README.md`'s Network & privacy section. Result depends on PubChem's database being current and correct.

### `retrosynthesis`
- Description: One-step BRICS disconnection, all breakable bonds cut individually, ranked by max fragment SA Score.
- Local/network: local
- Determinism: DeterministicByDesign — Result method: **HeuristicSearch** — Completeness: **BestEffort** (a single-step rule-based disconnection proposal, not a complete or exhaustive retrosynthetic search)
- High cost: Potentially, near the atom-count ceiling (BRICS enumeration + per-bond BFS + 2× SA-score per breakable bond).
- Explicit limits: ≤500 atoms, single connected component required.
- Caveats: do not describe this as "the" retrosynthesis or a validated route — it is a rule-based (BRICS) one-step disconnection proposal that does not account for reagent availability, reaction feasibility, or stereochemistry, and does not search beyond one disconnection step.

### `smiles_to_moljson`
- Description: SMILES → MolJSON.
- Local/network: local
- Determinism: DeterministicByDesign — Result method: Algorithmic — Completeness: Complete
- High cost: No
- Explicit limits: none

### `moljson_to_smiles`
- Description: MolJSON → canonical SMILES.
- Local/network: local
- Determinism: DeterministicByDesign — Result method: Algorithmic — Completeness: Complete
- High cost: No
- Explicit limits: **none** — no length cap on the input `json` string.

### `representation_router`
- Description: Route SMILES to the molecular text representation (MolJSON/CML/InChI/canonical SMILES) best suited to a given LLM task.
- Local/network: local
- Determinism: DeterministicByDesign — Result method: **mixed** (the underlying format conversions are Algorithmic/exact; the task→representation mapping itself is a heuristic recommendation from one specific paper (arXiv 2026), not a universal standard) — Completeness: Complete
- High cost: No
- Explicit limits: none
- Caveats: the "best" representation choice reflects one paper's recommendation, not a guarantee that it is optimal for every downstream LLM/task.

### `molecule_context_pack`
- Description: Identifiers, physicochemical properties, drug-likeness flags, ADMET profile, structural alerts, and MolJSON in one object (json/markdown/prompt output).
- Local/network: local
- Determinism: DeterministicByDesign — Result method: **mixed** (bundles exact identifiers/MolJSON with empirical LogP/QED/SA/ADMET/PAINS/Brenk estimates) — Completeness: Complete
- High cost: No (aggregates several already-covered lighter calls; no single expensive step)
- Explicit limits: none
- Caveats: the most heavily mixed tool in the catalog — inherits every caveat listed above for the individual estimates it bundles (Crippen LogP, QED, `sa_score`, `admet_profile`, PAINS/Brenk), none of which are distinguished from the exact fields in the output.

## Known cross-cutting gaps (see the PR 1 audit for full detail)

- No tool declares `readOnlyHint`/`destructiveHint` (all 20 are in fact read-only/non-destructive, but this isn't machine-declared).
- No tool declares the network-access / determinism / completeness axes above as structured metadata — this document is the only place they exist today.
- Only `find_mcs`, `retrosynthesis`, and `name_to_smiles` have any explicit input/size limit; the other 17 are unbounded.
- `generate_3d` is both high-cost and has no limit at all — the largest concrete resource-risk gap in the catalog.
