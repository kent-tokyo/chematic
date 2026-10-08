# Probes for xsmarts-autoconf (#754)

`smirks-probes-v1.json` holds minimized SMIRKS behaviour probes found by
chematic's property-based fuzz lane (`scripts/smirks_property_fuzz.py`) and
the #754 comparison. Each probe is a SMIRKS, a reactant SMILES, a hydrogen
mode and RDKit 2026.03.6's sanitized product sets, so it runs on any engine
without chematic. They are offered to
[xsmarts-autoconf](https://github.com/swamidasslab/xsmarts-autoconf) as
candidate behaviour flags; nothing here is vendored from that project.

| Probe | What it pins |
|---|---|
| `smarts.lowercase_pair_is_aromatic_atom_plus_primitive` | `[ca]` is aromatic carbon AND aromatic, not calcium |
| `smarts.lowercase_ring_size` | `[cr6]` is a ring-size primitive, not chromium |
| `smirks.element_change_resets_charge` | an element change takes the template's charge |
| `smirks.element_change_recomputes_h` | and re-derives hydrogens for the new element |
| `smirks.element_change_keeps_aromatic_flag` | the rewritten ring atom stays aromatic and kekulizes hypervalent |
| `smirks.explicit_h_element_change` | `AddHs` reactants keep the H atom on the new element |
| `sanitize.kekulize_depends_on_atom_order` | RDKit's Kekulé structure of a ring cation depends on atom order |

Regenerate the RDKit products with the snippet in
`benchmarks/2026-10-08-754-734-followups-batch24.md`.
