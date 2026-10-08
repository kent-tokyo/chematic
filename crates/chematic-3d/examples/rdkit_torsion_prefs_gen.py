import re, sys
base = "/tmp/claude-0/-home-claude/c9395eac-d598-5a86-b7de-b9824d14929d/scratchpad/rdsrc/rdkit/Code/GraphMol/ForceFieldHelpers/CrystalFF/"
def literals(path):
    src = open(base + path).read()
    # strip /* */ and // comments outside strings
    out, i, n = [], 0, len(src)
    lits = []
    while i < n:
        c = src[i]
        if src.startswith("/*", i):
            i = src.index("*/", i) + 2; continue
        if src.startswith("//", i):
            i = src.index("\n", i); continue
        if c == '"':
            j = i + 1; s = []
            while src[j] != '"':
                if src[j] == "\\":
                    e = src[j+1]
                    s.append({"n": "\n", "t": "\t", "\\": "\\", '"': '"'}[e]); j += 2
                else:
                    s.append(src[j]); j += 1
            lits.append("".join(s)); i = j + 1; continue
        i += 1
    return "".join(lits)
def rust_str(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n") + '"'
out = ["//! ETKDG experimental torsion preferences, verbatim from RDKit 2026.03",
       "//! (`Code/GraphMol/ForceFieldHelpers/CrystalFF/torsionPreferences_*.in`,",
       "//! BSD licence; J. Chem. Inf. Model. 56, 1 (2016) and 10.1021/acs.jcim.0c00025).",
       "//! Generated; do not edit by hand.", ""]
for name, f in (("V1", "torsionPreferences_v1.in"), ("V2", "torsionPreferences_v2.in"),
                ("SMALL_RINGS", "torsionPreferences_smallrings.in"), ("MACROCYCLES", "torsionPreferences_macrocycles.in")):
    s = literals(f)
    out.append(f"pub(crate) const TORSION_PREFERENCES_{name}: &str = concat!(")
    for line in s.splitlines(keepends=True):
        out.append("    " + rust_str(line) + ",")
    out.append(");\n")
open(sys.argv[1], "w").write("\n".join(out))
for f in ("torsionPreferences_v2.in", "torsionPreferences_macrocycles.in"):
    s = literals(f); print(f, len(s.splitlines()), sum(1 for l in s.splitlines() if l.startswith("#")))
