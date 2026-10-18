# RDKit Cleanup: isolated alkali with one hydroxyl site

Date: 2026-10-18. Candidate source: `5fec6a07`.
The comparator is RDKit `2026.09.1` on the pinned 4,072-row stress corpus.
The [row-level result](2026-10-18-rdkit-cleanup-alkali-oxygen-followup.json)
records source hashes, the corpus hash, and every remaining difference.

Both parsers agree on all 4,072 inputs. RDKit Cleanup accepts 4,069 of them.
Exact Cleanup output agreement increases from **3,922 to 3,923/4,069
(96.41%)** against the previous source result. Row 2796 improves, zero rows
regress, and 146 output differences remain.

The previous isolated-metal rule charged a bare K atom in the same mixture as
a phenol but left the phenol neutral. The output gained a net positive formal
charge. RDKit pairs the new K cation with a phenoxide anion in that row. This
follow-up applies that charge pairing only to one originally isolated bare
Li/Na/K atom and one unambiguous neutral carbon-bound OH site in an otherwise
uncharged mixture. Multiple metal atoms, competing OH sites, explicit OH
hydrogen counts, and already charged mixtures remain outside this rule.
Focused Rust tests check alcohol, phenol, carboxylic acid, water, an existing
Na cation, competing OH sites, and Cleanup idempotence.

Reproduce with a Python environment containing RDKit `2026.09.1`:

```sh
cargo build -p chematic-chem --example cleanup_rdkit_dump --offline
target/debug/examples/cleanup_rdkit_dump \
  < validation/cosmolkit_comparison/corpora/stress-v1.smi \
  > /tmp/chematic-cleanup-alkali-oxygen.jsonl
python scripts/rdkit_cleanup_source_gate.py \
  --dump /tmp/chematic-cleanup-alkali-oxygen.jsonl \
  --previous-source benchmarks/2026-10-18-rdkit-cleanup-isolated-metals-followup.json \
  --output /tmp/chematic-cleanup-alkali-oxygen-gate.json
```

Local checks: `chematic-chem` unit tests 985 passed and 1 ignored; Clippy,
rustfmt, and refactor budget passed. This is output agreement on a selected
corpus, not independent chemical truth or broad accuracy superiority over
RDKit. Reaction and R-group source behavior is unchanged by this follow-up.

## Charge conservation in two remaining disagreements

The frozen row-level result also exposes two simple metal-disconnection cases
where exact RDKit agreement would lose a useful invariant. For row 410,
`O=[Mn](=O)(=O)[O-].[K+]` has total formal charge 0. RDKit Cleanup produces
`[K+].[Mn+8].[O-2].[O-2].[O-2].[O-2]` (charge +1); CheMatic produces the same
separated ions with `[Mn+7]` (charge 0). For row 455, `C(#N)[Cu]` starts at
charge 0. RDKit produces `[C-]#N.[Cu]` (charge -1); CheMatic produces
`[C-]#N.[Cu+]` (charge 0). A source unit test now pins both outputs and checks
charge conservation. These are two narrow, independently checkable advantages
on the recorded examples; they do not establish general Cleanup superiority.

The [formal-charge audit](2026-10-18-rdkit-cleanup-charge-invariant-audit.json)
reconstructs all 4,069 comparable RDKit `2026.09.1` Cleanup outputs from the
frozen comparison and verifies the source dump hash. RDKit `2025.09.3` was used
only to parse the recorded SMILES and sum formal charges. Of 4,066 rows that
parser accepts in all three forms, both outputs preserve input charge in 4,061,
only CheMatic does in 3 (indices 410, 455, and 3514), only RDKit does in 0,
and neither does in 2. Three rows cannot be scored by that parser. The two
rows where neither preserves charge are bare-metal cases where both tools add
charge. Charge preservation is one necessary check for the two disconnection
examples above, not a complete correctness metric.

After generating the source dump with the command above, reproduce the audit:

```sh
python scripts/rdkit_cleanup_charge_invariant_audit.py \
  --dump /tmp/chematic-cleanup-alkali-oxygen.jsonl \
  --comparison benchmarks/2026-10-18-rdkit-cleanup-alkali-oxygen-followup.json \
  --output /tmp/chematic-cleanup-charge-audit.json
```
