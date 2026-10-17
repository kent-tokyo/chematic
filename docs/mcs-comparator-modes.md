# Typed MCS comparator modes

`find_mcs` and `find_mcs_checked` accept the Python string enums
`McsAtomCompare` and `McsBondCompare`. These typed choices address the enum side
of the comparator-input work tracked by
[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601). Existing valid string
values continue to work. RDKit comparator subclasses and Python callbacks are
outside this API's supported modes.

```python
import chematic as ch

molecules = [ch.from_smiles("CCC"), ch.from_smiles("NNN")]
result, timed_out = ch.find_mcs_checked(
    molecules,
    atom_compare=ch.McsAtomCompare.ANY_HEAVY_ATOM,
    bond_compare=ch.McsBondCompare.ORDER_OR_AROMATIC,
)
assert not timed_out
assert result is not None and result.atom_count == 3
assert ch.find_mcs(molecules, atom_compare=ch.McsAtomCompare.ELEMENTS) is None
```

| Atom mode | Legacy string | Matching rule |
| --- | --- | --- |
| `McsAtomCompare.ELEMENTS` | `"elements"` | Equal atomic number; the default |
| `McsAtomCompare.ANY_HEAVY_ATOM` | `"any_heavy_atom"` | Ignore element identity, exclude hydrogen |
| `McsAtomCompare.ANY` | `"any"` | Ignore element identity, include hydrogen |

| Bond mode | Legacy string | Matching rule |
| --- | --- | --- |
| `McsBondCompare.ORDER_OR_AROMATIC` | `"order_or_aromatic"` | Equal normalized order; directional single bonds normalize to single |
| `McsBondCompare.ANY` | `"any"` | Ignore bond order, retain required connectivity |

The modes are native chematic choices. They do not supply RDKit's isotope-only
atom comparator, arbitrary comparison callbacks, or RDKit's distinction between
all bond-order comparison modes. `match_bonds=False` continues to bypass bond
order comparison. Existing charge, isotope, ring and chirality options retain
their native behavior.

`ANY_HEAVY_ATOM` retains chematic's existing heavy-only meaning. RDKit's
`CompareAnyHeavyAtom` instead permits hydrogen-to-hydrogen matches as well as
arbitrary heavy-to-heavy matches. Explicit-H inputs can therefore have different
MCS sizes even with similarly named modes. Use `ANY` when all atom identities
including hydrogen may be ignored.

Use the public aliases `McsAtomComparison` and `McsBondComparison` to annotate
variables accepting either an enum or a valid literal string. Invalid literal
strings and the wrong enum family are rejected by type checkers. Runtime remains
compatible with strings: unknown values raise `ValueError`, and unsupported
non-string objects raise `TypeError`. As string subclasses, enum values retain
string equality; runtime does not distinguish a wrong-family enum whose value is
also a valid string in the target family. Convert dynamic configuration using
`McsAtomCompare(value)` or `McsBondCompare(value)` to validate its vocabulary.

## Comparator-aware search bounds

The native MCS search previously bounded every branch using element-specific
atom counts. That underestimates possible matches when element identity is
ignored: `CCC` and `NNN` returned one atom, while the common three-atom topology
is valid. `CCCO` and `NCCN` similarly returned three instead of four atoms.

Both the seed bound and the remaining-frontier bound now use the selected atom
comparison. `ELEMENTS` retains element buckets. `ANY` uses one shared bucket;
`ANY_HEAVY_ATOM` uses one bucket for non-hydrogen atoms. Additional charge/isotope/
ring constraints are ignored by the bound, making it conservative; the existing
compatibility checks still enforce them during matching. Search-node and time
limits remain in place. Searches allowing more atom pairings may explore more
branches; this change makes no performance claim.

The Python return remains a concrete `Mol` projected from the first input side.
It is not a generalized SMARTS query encoding every permitted atom/bond choice
and must not be assumed to match every input after reparsing its SMILES. This
change validates MCS **size/counts**, not generalized query output, atom mappings,
canonical winners, or full RDKit FMCS equivalence. The existing checked boolean
reports timeouts; it does not report the native Rust resource-limit outcome.

The Rust comparison modes are `chematic_smarts::AtomCompare` and `BondCompare`.
`find_mcs_with_config_checked` retains the separate `McsOutcome` variants for
exhaustive, timed-out and resource-limited searches.
