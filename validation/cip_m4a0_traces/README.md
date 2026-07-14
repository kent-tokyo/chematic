# M4A-0 comparison traces

One representative `ComparisonTrace` per bucket (Rules-1a/2 decision path for the tied
stereocenter's root-children ranking), generated via
`cargo run -p chematic-cip --release --example trace_report -- '<smiles>' <atom_idx>`.
Each is the *shortest* corpus residual row in its bucket, used as the minimal
reproducible case rather than a constructed synthetic molecule.

| bucket | minimal SMILES | atom_idx | trace file |
|---|---|---|---|
| phosphorus | `CN[P@@]1(N)=NP(N2CC2)(N2CC2)=N[P@](N)(NC)=N1` | 2 | `phosphorus.txt` |
| rule5_pseudoasymmetry | `O=c1cc(COC2CCOCC2)occ1OC(=O)[C@]12C[C@H]3C[C@H](C[C@H](C3)C1)C2` | 20 | `rule5_pseudoasymmetry.txt` |
| rule4_candidate | `O=C(O[C@H]1[C@H](O)C[C@](O)(C(=O)O)C[C@H]1O)c1cc(O)c(O)c(O)c1` | 3 | `rule4_candidate.txt` |

Notable, directly observable in the trace text itself (not just inferred from bucket
labels):

- **phosphorus**: the trace fully resolves the digraph to 5 distinct rank groups (no
  tie) -- the residual isn't a rule-insufficiency tie, it's a *definite but wrong*
  ranking (a Rule 1a/2 correctness bug in how the comparator handles this P=N
  duplicate-node structure), distinct from the other two buckets' genuine ties.
- **rule4_candidate**: the top-level ring-branch pair resolves to `Equal (rule leaf)`
  after just 1 level of comparison -- a shallow, clean tie.
- **rule5_pseudoasymmetry**: the same shape of tie, but the comparator explores 424
  decision steps (vs. the quinic case's handful) before exhausting Rules 1a/2 -- the
  cage's fused tricyclic ring system requires far deeper recursion before the tie
  surfaces, consistent with it being a genuinely deep constitutional tie (see
  `crates/chematic-cip/src/tests.rs::diagnose_m4a0_quinic_residual_constitutional_identity`
  for the quinic case's own structural-identity confirmation via `branch_signature`).
