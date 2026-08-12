# 3D Breakthrough Program — master plan

Status: **Wave 0 complete (this document + frozen baseline). Wave 1 not yet
launched.** This document is the single source of truth for wave order, file
ownership, and the architecture decisions that unblock Wave 1's parallel
agents. It supersedes ad hoc restatement of the program brief in any PR
description — PRs should link here, not re-derive it.

Coordinator has **not** modified any file under `crates/*/src/**` in Wave 0.
This document, `validation/manifests/3d_breakthrough_baseline.json`, and
`validation/manifests/dataset_provenance.json` are the only new artifacts.

## 0. Wave 0 findings (read before touching anything)

- `docs/rfcs/etkdg_3d_gap_rfc.md`'s diagnosis was re-verified, not trusted on
  faith: `scripts/etkdg_vs_rdkit_gap.py` was rerun as **5 independent fresh
  Python processes** against current HEAD (`9f9f459`, v0.7.0 — the RFC's own
  numbers were measured against an older 0.6.0-era build). All 5 runs
  reproduced the RFC's numbers byte-identically:
  `geometrically_valid_rate = 25.9%` (15/58 ok, 39/58 bond-length-blowup,
  4/58 gross-clash). **The diagnosed defect is confirmed live and unchanged
  at v0.7.0**, not a stale claim. Full detail, denominators, and the
  ensemble-metric spread (5-run observation, see below) are in
  `validation/manifests/3d_breakthrough_baseline.json`. The original
  RFC-committed JSON/JSONL were preserved at
  `validation/results/etkdg_baseline_history/etkdg_vs_rdkit_{summary,rows}.RFC-original.{json,jsonl}`
  before being overwritten by the rerun.
- All 11 files/directories the task asked to confirm exist and were verified
  by direct file check (not assumed): `docs/rfcs/etkdg_3d_gap_rfc.md`,
  `scripts/etkdg_vs_rdkit_gap.py`, `crates/chematic-3d/src/{dg_fft.rs,
  etkdg.rs, etkdg_knowledge.rs, conformer.rs, stereo3d.rs}`,
  `crates/chematic-ff/src/`, `crates/chematic-inchi/src/dedup.rs`,
  `crates/chematic-mol/src/{mol2000.rs,mol3000.rs}`.
- **`stereo3d.rs` is not dead code** — `assign_stereo_from_3d` is exported
  and used by the gap script itself to *perceive* stereo from a geometry.
  This is the opposite direction from Agent D's job (*enforcing* declared
  stereo as an embedding constraint). Agent D must not assume this file
  already solves its problem.
- **Ensemble-metric reproducibility is observed, not contracted.** 5
  independent process runs produced bit-identical `chematic_n_kept` and
  diversity numbers for every `ENSEMBLE_SUBSET` molecule, despite the code's
  own comment stating the process-global-atomic-counter PRNG is "not
  bit-reproducible run-to-run." Most likely explanation: the counter starts
  from the same fixed value every fresh process and today's call pattern is
  single-threaded/call-order-deterministic — this is almost certainly
  coincidental to current usage, not a contract. **Agent C must still add an
  explicit, public, call-local seed API** (`EmbedParameters.random_seed`) —
  observed stability today does not substitute for one.
- Open PRs at freeze time: #126/#127/#128 (I/O streaming formats, unrelated),
  #158 (3D scalar descriptors, `chematic-3d` — **potential file-overlap risk
  with Agent G**, see §6), #159 (SMARTS, unrelated). Open issues: #70, #90,
  #91, #92, #107, #139, #149, #161 (none block this program directly).
- `git status` was clean at freeze except 3 unrelated README edits from a
  prior task in this session, which were committed and branched separately
  (`docs/npm-status-sync`, not part of this program).
- `main` is protected (`GH006` on direct push, confirmed earlier this
  session) — Wave 0's own artifacts ship via `chore/3d-breakthrough-wave0-baseline`
  → draft PR → CI green → merge, the same shape as every prior release
  branch. No wave in this program pushes directly to `main`.

## 1. Two architecture decisions that unblock Wave 1 (resolved here, not left to agents)

### 1a. Where does the 3D coordinate type live? (unblocks Agent A)

**Decision: a new, minimal `Coords3D`-shaped type lives in `chematic-core`,
not in `chematic-mol` importing `chematic-3d`.**

`chematic-3d`'s existing `Coords3D` (`crates/chematic-3d/src/coords.rs`) is
unsuitable to import directly into `chematic-mol`: per `CLAUDE.md`'s layering,
`chematic-mol → core, rxn, perception, smiles` has no edge to `chematic-3d`,
and `chematic-3d` itself pulls `chematic-ff`/`chematic-chem`/`chematic-fp`/
`chematic-smarts` — adding that whole subtree as a `chematic-mol` dependency
just to store an `(x, y, z)` array per atom would be a large, unwanted
footprint increase for every consumer of `chematic-mol` (Python, WASM, CLI
tools that only read/write MOL files and never touch 3D generation).

Every crate already depends on `chematic-core`, so putting the plain type
there requires **zero new dependency edges** and needs no
`scripts/check_publish_graph.py` change. Concretely:

- New `crates/chematic-core/src/coords3d.rs`: a plain, algorithm-free
  `pub struct Coords3D(pub Vec<[f64; 3]>)` (name and exact shape are Agent
  A's to finalize with Coordinator sign-off, but it must stay a dumb data
  holder — no embedding logic, no force-field code, nothing that would pull
  `chematic-3d`'s dependencies backward into `chematic-core`).
  Agent A implements this file as part of its own PR (it is additive to
  `chematic-core` and does not touch any Coordinator-only file other than
  one new `pub mod coords3d;` line in `chematic-core/src/lib.rs`, which Agent
  A may add directly since `chematic-core/src/lib.rs` is not on the
  Coordinator-only list).
- `chematic-3d`'s own `Coords3D` (`crates/chematic-3d/src/coords.rs`) is
  **Agent C's file** (see §3) and should grow a `From`/`Into` (or thin
  wrapper) relationship with `chematic_core::Coords3D` rather than being
  replaced outright in Wave 1 — Coordinator will reconcile the two into one
  canonical type at Wave 1→Wave 2 integration time, once both Agent A and
  Agent C's PRs exist and can be diffed together.
- `MolReadReport.conformer: Option<chematic_core::Coords3D>` — additive
  field, existing `MolReadReport` consumers are unaffected (already the
  additive-API precedent from `StereoDiagnostic`/`EzDirectionDiagnostic` in
  v0.7.0).

### 1b. How does Agent C measure its own acceptance gate, given `etkdg.rs` is Coordinator-only? (unblocks Agent C)

**Decision: Agent C's PR is measured through its own new module and its own
new example/test harness — never through `etkdg.rs` or `conformer_ensemble`.**
Wiring the new distance-geometry embedder into the live `etkdg.rs` default
path is explicitly a **Wave 2 integration step performed by Coordinator**,
consistent with the task's own ordering rule ("Cがgreenになる前にEのtorsionを
live defaultへ接続してはいけません" — i.e. C's output isn't live-defaulted
during Wave 1 in the first place).

Concretely:

- Agent C owns a new file, `crates/chematic-3d/src/distance_geometry_v2.rs`,
  with its own public entry point (shape sketched in the original brief,
  Agent C finalizes): `pub fn embed_distance_geometry_v2(mol: &Molecule,
  params: &EmbedParameters) -> Result<Coords3D, EmbedFailureCause>`.
- Coordinator adds **one line**, `pub mod distance_geometry_v2;`, to
  `chematic-3d/src/lib.rs` as a narrow, additive, pre-Wave-1 scaffolding
  commit (tracked as its own tiny PR, opened immediately before Agent C
  starts — not part of Wave 0's docs-only PR, and not "production code" in
  the sense Wave 0 forbids, since it adds no behavior, only a module
  declaration for code that doesn't exist yet).
- Agent C's own PR includes a new, self-contained validation entry point it
  owns outright (e.g. `crates/chematic-3d/examples/distance_geometry_v2_gap_check.rs`,
  mirroring the existing `stereo2d_fixture_dump.rs`-style example pattern
  used elsewhere in this codebase) that calls
  `embed_distance_geometry_v2` directly against the frozen 58-molecule
  corpus and reports `geometrically_valid_rate` the same way
  `scripts/etkdg_vs_rdkit_gap.py` does — but through the new function, not
  through `Mol.conformer_ensemble()`.
- **Agent C's acceptance gate is scoped narrower than the original brief's
  blanket "100% on frozen 58":** see §4 (gate scoping, resolves the RFC's
  unattributed-residual risk the advisor flagged).

## 2. "RDKitを超えた" — unchanged from the original brief

Quality win and Pareto win definitions, and the blinded-holdout requirement,
are as specified in the program brief. Restated in one line since it governs
every wave: **no README/CHANGELOG claim of RDKit-superiority before Wave 3's
blinded benchmark completes** — this document does not relax that gate.

## 3. Hot-file ownership (updated — 3 files the original table left unassigned)

```text
crates/chematic-3d/src/etkdg.rs           Coordinator only
crates/chematic-3d/src/lib.rs             Coordinator only (Agent C gets exactly
                                           one pre-authorized additive line, see §1b)
workspace Cargo.toml                      Coordinator only
CHANGELOG / README                        Coordinator only
crates/chematic-core/src/coords3d.rs      Agent A (new file; also adds one
                                           `pub mod` line to chematic-core/src/lib.rs)
mol2000.rs / mol3000.rs / sdf.rs          Agent A
dedup.rs / native InChI                   Agent B
dg_fft.rs / distance_geometry_v2.rs       Agent C (new module, see §1b)
crates/chematic-3d/src/prng.rs            Agent C (EmbedParameters.random_seed
                                           requires replacing the process-global RNG)
crates/chematic-3d/src/coords.rs          Agent C (existing chematic-3d Coords3D;
                                           reconciled with chematic-core's type at
                                           Wave 1→2 integration, see §1a)
stereo constraint modules                 Agent D
etkdg_knowledge.rs                        Agent E
crates/chematic-3d/src/minimize.rs        Agent F (the crippled bond+angle+VdW-only
                                           energy function — force-field bridge target)
force-field bridge (new modules)          Agent F
conformer.rs / pruning                    Agent G
validation scripts/results                Agent H
Python/WASM bindings                      Agent I
```

No agent edits another agent's owned files. Cross-cutting exports/integration
are Coordinator's job, performed at wave boundaries, not mid-wave.

**Known overlap risk**: open draft PR #158 ("RDKit-compatible G1 3D scalar
descriptors") touches `chematic-3d`. Coordinator will check its diff against
this ownership table before Agent G starts, and will not merge #158 into
`main` while Agent G's PR is open without first diffing for file collisions.

## 4. Agent C's acceptance gate — scoped, not the blanket 100% from the original brief

The RFC pins the 67% blow-up to **three** mechanisms, only one of which is
Agent C's:

1. Under-iterated constraint repair (`etkdg.rs` / `constraints.rs`) — not
   Agent C's file.
2. Ring-unaware torsion rotation (`etkdg_knowledge.rs` / `mol_transforms.rs`)
   — Agent E's file.
3. Silent missing-MMFF94-parameter zero-gradient (`minimize.rs`) — Agent F's
   file.

Plus an **explicitly unattributed** acyclic residual (decane, hexadecane,
but-2-ene, 2-chlorobutane, l-serine/l-threonine blow up under *both* force
fields) that the RFC itself declines to pin to a single mechanism.

Asking Agent C alone to hit "100% geometrically valid on frozen 58" bakes in
failure classes that are structurally not its file's responsibility. Revised
gate for Agent C's Wave 1 merge:

- **First deliverable (cheap, before writing the embedder): a short
  root-cause note attributing the acyclic residual** — either to the same
  local/non-rigid constraint-projection weakness as mechanism 1, or to
  something new. This de-risks the 100% claim before it's made.
- **Merge gate measured pre-minimization, force-field-agnostic**: raw output
  of `embed_distance_geometry_v2` (bounds construction → smoothing → Gram/
  eigendecomposition → bounds-force refinement, no MMFF94/DREIDING
  minimization pass) must be `geometrically_valid_rate = 100%` on frozen 58
  under the same >50%-covalent-radius-blowup check, evaluated on the raw
  embedding directly — this isolates Agent C's actual deliverable (real
  distance geometry replacing the DFS placer) from Agent F's (force-field
  minimization) and Agent E's (torsion knowledge, deferred to Wave 2 per the
  wave order).
- The full, original "100% after the whole live pipeline including MMFF94
  minimization" gate is a **Wave 2 integration gate**, owned by Coordinator,
  measured only after Agent C + Agent F (+ Agent E, rings/macrocycles from
  Agent G) are all merged.

## 5. Everything else from the original program brief is unchanged

Sections not restated here (Wave 1 agent briefs B/F, Wave 2 D/E/G, Wave 3 H,
Wave 4 I/J, no-silent-fallback rules, shared-environment prohibition, PR
rules, common tests, final release gate) apply exactly as specified in the
program brief already delivered to Coordinator. This document only resolves
the two decisions that were genuine blockers and updates the ownership table
and Agent C's gate accordingly. Agents B and F need no further architectural
resolution and are launchable as originally briefed.

## 6. Wave 1 launch sequence (next action after this document merges)

1. Coordinator opens the one-line `distance_geometry_v2` module-stub PR
   (§1b) against `main` — small enough to merge same-day.
2. Coordinator dispatches Agent A, Agent B, Agent C, Agent F in parallel
   (≤4 concurrent), each in its own worktree + branch + draft PR, per the
   original brief's Wave 1 section, using this document's §1/§3/§4
   resolutions as their brief instead of the ambiguous originals.
3. No Wave 1 PR merges to `main` until Coordinator has independently
   reviewed it against this document's ownership table and acceptance gates.
