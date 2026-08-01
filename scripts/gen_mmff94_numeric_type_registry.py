#!/usr/bin/env python3
"""Generate crates/chematic-ff/src/mmff94_numeric_type_registry.rs from the
frozen RDKit source tables in scripts/mmff94_provenance/.

Parses `rdkit_defaultMMFFDef.txt` and `rdkit_defaultMMFFProp.txt` using
*exactly* the same rules as RDKit's own C++ parsers
(`MMFFDefCollection`/`MMFFPropCollection` in `Code/ForceField/MMFF/Params.cpp`
at commit e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f, tag Release_2026_03_3 --
see scripts/mmff94_provenance/PROVENANCE.md), not a re-derivation: lines
starting with '*' are secondary/alias symbols and are skipped, matching
`inLine[0] != '*'` in RDKit's parser.

This is a generated-file generator per the program's "generated results are
regenerated from the generator, never hand-edited" rule -- do not edit
mmff94_numeric_type_registry.rs directly.

Run: python3 scripts/gen_mmff94_numeric_type_registry.py
"""

import argparse
import re
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
PROV = ROOT / "scripts" / "mmff94_provenance"
OUT = ROOT / "crates" / "chematic-ff" / "src" / "mmff94_numeric_type_registry.rs"

# Atomic-number -> chematic_core::Element variant name, only for elements
# that actually appear in MMFF94 (defaultMMFFProp `aspec` column values).
ATOMIC_NUMBER_TO_ELEMENT = {
    1: "H", 3: "LI", 6: "C", 7: "N", 8: "O", 9: "F", 11: "NA", 12: "MG",
    14: "SI", 15: "P", 16: "S", 17: "CL", 19: "K", 20: "CA", 26: "FE",
    29: "CU", 30: "ZN", 35: "BR", 53: "I",
}


def extract_cstring_array(content, varname):
    """Extract and concatenate a C++ `const std::string NAME = "a" "b" ...;`
    literal's parts, decoding C escapes. Matches only consecutive quoted
    string tokens (not a naive "up to the first semicolon" scan) -- some
    RDKit MMFF description text contains a literal `;` character (e.g. type
    55's "N IN +N=C-N:;\tQ=1/2"), which would otherwise truncate the match.
    """
    m = re.search(r"const std::string " + varname + r"\s*=\s*", content)
    if not m:
        raise ValueError("not found: " + varname)
    pos = m.end()
    token_re = re.compile(r'\s*"((?:[^"\\]|\\.)*)"')
    parts = []
    while True:
        tm = token_re.match(content, pos)
        if not tm:
            break
        parts.append(tm.group(1))
        pos = tm.end()
    raw = "".join(parts)
    return raw.encode().decode("unicode_escape")


def parse_def(text):
    """Returns {type_id: (symbol, eq_levels[4], description)}."""
    out = {}
    for line in text.split("\n"):
        if not line or line[0] == "*":
            continue
        cols = line.split("\t")
        symbol = cols[0]
        type_id = int(cols[1])
        eq_levels = tuple(int(x) for x in cols[2:6])
        description = " ".join(c for c in cols[6:] if c)
        if type_id not in out:
            out[type_id] = (symbol, eq_levels, description)
    return out


def parse_prop(text):
    """Returns {type_id: (atno, crd, val, pilp, mltb, arom, lin, sbmb)}."""
    out = {}
    for line in text.split("\n"):
        if not line or line[0] == "*":
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        type_id = int(cols[0])
        out[type_id] = tuple(int(x) for x in cols[1:9])
    return out


def refresh_frozen_tables(params_cpp_path):
    """Re-extract the frozen rdkit_defaultMMFF{Def,Prop}.txt tables from a
    local `Code/ForceField/MMFF/Params.cpp` copy at the pinned commit -- see
    scripts/mmff94_provenance/PROVENANCE.md for how to fetch that copy."""
    content = pathlib.Path(params_cpp_path).read_text()
    defs_raw = extract_cstring_array(content, "defaultMMFFDef")
    props_raw = extract_cstring_array(content, "defaultMMFFProp")
    (PROV / "rdkit_defaultMMFFDef.txt").write_text(defs_raw)
    (PROV / "rdkit_defaultMMFFProp.txt").write_text(props_raw)
    print(f"refreshed {PROV / 'rdkit_defaultMMFFDef.txt'} and rdkit_defaultMMFFProp.txt")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--extract",
        metavar="PARAMS_CPP",
        help="re-extract the frozen provenance .txt tables from a local "
        "Params.cpp copy before generating (see PROVENANCE.md)",
    )
    args = ap.parse_args()
    if args.extract:
        refresh_frozen_tables(args.extract)

    defs = parse_def((PROV / "rdkit_defaultMMFFDef.txt").read_text())
    props = parse_prop((PROV / "rdkit_defaultMMFFProp.txt").read_text())

    type_ids = sorted(set(defs) | set(props))
    rows = []
    for tid in type_ids:
        symbol, eq_levels, description = defs.get(tid, ("?", (0, 0, 0, 0), "unknown"))
        prop = props.get(tid)
        if prop is None:
            # A handful of MMFFDEF.PAR entries (metal cations, ions) have no
            # MMFFPROP.PAR row in RDKit's own tables either -- skip rather
            # than fabricate property data chematic never assigns anyway.
            continue
        atno, crd, val, pilp, mltb, arom, lin, sbmb = prop
        element = ATOMIC_NUMBER_TO_ELEMENT.get(atno)
        if element is None:
            continue
        rows.append((tid, symbol, element, atno, crd, val, bool(pilp), mltb,
                      bool(arom), bool(lin), bool(sbmb), eq_levels, description))

    lines = []
    lines.append("//! GENERATED FILE -- DO NOT EDIT.")
    lines.append("//!")
    lines.append(
        "//! Regenerate with `python3 scripts/gen_mmff94_numeric_type_registry.py`."
    )
    lines.append(
        "//! Source: RDKit `Code/ForceField/MMFF/Params.cpp` "
        "(`defaultMMFFDef` + `defaultMMFFProp`), tag `Release_2026_03_3`,"
    )
    lines.append(
        "//! commit `e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f`. Full provenance: "
        "`scripts/mmff94_provenance/PROVENANCE.md`."
    )
    lines.append("//!")
    lines.append(
        "//! Authoritative metadata for every MMFF94 numeric atom type chematic can"
    )
    lines.append(
        "//! assign or look up. `element` is the semantic-compatibility gate's ground"
    )
    lines.append(
        "//! truth: a parameter row matched under a type whose registered `element`"
    )
    lines.append(
        "//! disagrees with the atom's real element is a collision, not a hit (issue"
    )
    lines.append("//! #227's `furan` finding).")
    lines.append("")
    lines.append("use chematic_core::Element;")
    lines.append("")
    lines.append("/// One MMFF94 numeric atom type's authoritative metadata.")
    lines.append("#[derive(Debug, Clone, Copy)]")
    lines.append("pub struct Mmff94NumericTypeInfo {")
    lines.append("    /// MMFF94 numeric type ID (1-99).")
    lines.append("    pub id: u8,")
    lines.append("    /// RDKit/Halgren short symbol, e.g. \"CB\", \"C5A\", \"NPYD\".")
    lines.append("    pub symbol: &'static str,")
    lines.append("    /// The only element this numeric type may legitimately represent.")
    lines.append("    pub element: Element,")
    lines.append("    /// Atomic number (redundant with `element`, kept for direct")
    lines.append("    /// comparison against RDKit's `aspec` column without a round-trip).")
    lines.append("    pub atomic_number: u8,")
    lines.append("    /// Coordination number (MMFFPROP `crd`).")
    lines.append("    pub coordination: u8,")
    lines.append("    /// Valence (MMFFPROP `val`).")
    lines.append("    pub valence: u8,")
    lines.append("    /// Has a lone pair capable of conjugation (MMFFPROP `pilp`).")
    lines.append("    pub has_pi_lone_pair: bool,")
    lines.append("    /// Number of attached multiple bonds (MMFFPROP `mltb`, 0-3).")
    lines.append("    pub multiple_bond_count: u8,")
    lines.append("    /// Aromatic atom type (MMFFPROP `arom`).")
    lines.append("    pub aromatic: bool,")
    lines.append("    /// Linear geometry, e.g. nitrile/acetylenic (MMFFPROP `lin`).")
    lines.append("    pub linear: bool,")
    lines.append("    /// Single bond adjacent to a multiple bond (MMFFPROP `sbmb`).")
    lines.append("    pub single_bond_multiple_bond: bool,")
    lines.append("    /// MMFF94 atom-type equivalence levels 2-5 (`MMFFDEF.PAR`), used")
    lines.append("    /// for the equivalence-class parameter fallback (not yet wired")
    lines.append("    /// into any chematic resolver -- see the Phase 1B-0 PR body's")
    lines.append("    /// bond-fallback classification for why).")
    lines.append("    pub equivalence_levels: [u8; 4],")
    lines.append("    /// Human-readable Halgren description, e.g. \"AROMATIC C\".")
    lines.append("    pub description: &'static str,")
    lines.append("}")
    lines.append("")
    lines.append(f"pub static MMFF94_NUMERIC_TYPE_REGISTRY: &[Mmff94NumericTypeInfo] = &[")
    for (tid, symbol, element, atno, crd, val, pilp, mltb, arom, lin, sbmb,
         eq_levels, description) in rows:
        desc_escaped = description.replace("\\", "\\\\").replace("\"", "\\\"")
        lines.append(
            "    Mmff94NumericTypeInfo { "
            f"id: {tid}, symbol: \"{symbol}\", element: Element::{element}, "
            f"atomic_number: {atno}, coordination: {crd}, valence: {val}, "
            f"has_pi_lone_pair: {str(pilp).lower()}, multiple_bond_count: {mltb}, "
            f"aromatic: {str(arom).lower()}, linear: {str(lin).lower()}, "
            f"single_bond_multiple_bond: {str(sbmb).lower()}, "
            f"equivalence_levels: [{eq_levels[0]}, {eq_levels[1]}, {eq_levels[2]}, {eq_levels[3]}], "
            f"description: \"{desc_escaped}\" }},"
        )
    lines.append("];")
    lines.append("")
    lines.append("/// Look up a numeric MMFF94 type's authoritative metadata.")
    lines.append("///")
    lines.append("/// `MMFF94_NUMERIC_TYPE_REGISTRY` is sorted by ascending `id` (types")
    lines.append("/// are emitted in the source table's natural order, which is already")
    lines.append("/// ascending), so this binary-searches rather than scanning linearly.")
    lines.append("pub fn mmff94_numeric_type_info(id: u8) -> Option<&'static Mmff94NumericTypeInfo> {")
    lines.append("    MMFF94_NUMERIC_TYPE_REGISTRY")
    lines.append("        .binary_search_by_key(&id, |info| info.id)")
    lines.append("        .ok()")
    lines.append("        .map(|idx| &MMFF94_NUMERIC_TYPE_REGISTRY[idx])")
    lines.append("}")
    lines.append("")
    lines.append("#[cfg(test)]")
    lines.append("mod tests {")
    lines.append("    use super::*;")
    lines.append("")
    lines.append("    #[test]")
    lines.append("    fn registry_is_sorted_by_id() {")
    lines.append("        let ids: Vec<u8> = MMFF94_NUMERIC_TYPE_REGISTRY.iter().map(|i| i.id).collect();")
    lines.append("        let mut sorted = ids.clone();")
    lines.append("        sorted.sort_unstable();")
    lines.append("        assert_eq!(ids, sorted, \"generator must emit ascending id order\");")
    lines.append("    }")
    lines.append("")
    lines.append("    #[test]")
    lines.append("    fn benzene_type_37_is_aromatic_carbon() {")
    lines.append("        let info = mmff94_numeric_type_info(37).expect(\"type 37\");")
    lines.append("        assert_eq!(info.symbol, \"CB\");")
    lines.append("        assert_eq!(info.element, Element::C);")
    lines.append("        assert!(info.aromatic);")
    lines.append("    }")
    lines.append("")
    lines.append("    #[test]")
    lines.append("    fn furan_alpha_carbon_type_63_is_carbon_not_the_old_wrong_symbol() {")
    lines.append("        let info = mmff94_numeric_type_info(63).expect(\"type 63\");")
    lines.append("        assert_eq!(info.symbol, \"C5A\");")
    lines.append("        assert_eq!(info.element, Element::C);")
    lines.append("    }")
    lines.append("")
    lines.append("    #[test]")
    lines.append("    fn pyridine_n_type_38_is_nitrogen() {")
    lines.append("        let info = mmff94_numeric_type_info(38).expect(\"type 38\");")
    lines.append("        assert_eq!(info.symbol, \"NPYD\");")
    lines.append("        assert_eq!(info.element, Element::N);")
    lines.append("    }")
    lines.append("}")
    lines.append("")

    OUT.write_text("\n".join(lines))
    print(f"wrote {OUT} ({len(rows)} numeric types)")


if __name__ == "__main__":
    main()
