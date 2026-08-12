#!/usr/bin/env python3
"""Classify mmff94_strict's Angle/Bond "table_gap" instances (issue #227,
i.e. `mmff94_term_coverage_audit`'s `Angle`/`Bond` rows whose
`present_at_different_classification` is `null`) by which RDKit resolution
path -- source-verified, not guessed -- actually produces a value for the
same (angleType/bondType, atom-type-tuple).

Ported directly from the pinned RDKit source (Release_2026_03_3 /
e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f, see PROVENANCE.md):
  - Code/ForceField/MMFF/Params.h   (MMFFAngleCollection/MMFFBondCollection
    operator() -- the REAL table lookup, including Angle's 4-stage eqLevel
    canonical-type-substitution ladder, which chematic has no equivalent of
    at all; Bond has NO ladder, confirmed by direct read)
  - Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp (empirical rules,
    eq 18/19/20, MMFF.V page 625/627/628; isAngleInRingOfSize3or4)

Classification per unique tuple:
  - direct_table:     found at eqLevel stage 0 (types unchanged)
  - equivalence_table: found only after eqLevel substitution (stage 1-3)
  - empirical_rule:    not found in table at any stage -> RDKit's real
                       algorithm falls to getMMFFAngleBendEmpiricalRuleParams
                       / getMMFFBondStretchEmpiricalRuleParams (eq 18-20)
  - type_mismatch:     RDKit assigns a DIFFERENT numeric MMFF type to one of
                       i/j/k than chematic did for the same atom index --
                       not a genuine table gap, a chematic-side typing bug
  - unsupported:       neither table nor empirical rule can produce a value
                       (e.g. missing CovRadPauEle for an exotic element)

Live-oracle cross-check: for every non-type_mismatch tuple, also calls
GetMMFFAngleBendParams/GetMMFFBondStretchParams and, for empirical_rule
tuples, independently recomputes eq 18-20 in this script and diffs against
the live oracle value -- if they match, that's positive proof the empirical
path is what produced the number, not just an absence-of-table inference.

Usage:
    cargo run --release -p chematic-3d --example mmff94_term_coverage_audit \\
      > validation/results/mmff94_coverage_227_term_audit.jsonl \\
      2> /dev/null
    python3 scripts/mmff94_angle_bond_gap_classify.py \\
      validation/results/mmff94_coverage_227_term_audit.jsonl
"""

import collections
import json
import math
import re
import subprocess
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

RDKIT_SHA = "e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f"
PARAMS_CPP_URL = (
    f"https://raw.githubusercontent.com/rdkit/rdkit/{RDKIT_SHA}/"
    "Code/ForceField/MMFF/Params.cpp"
)

# ---------------------------------------------------------------------------
# Generic C++ adjacent-string-literal block extractor
# ---------------------------------------------------------------------------


def extract_cpp_string_block(text, var_name, is_array=False):
    """Extract and concatenate a `const std::string NAME = "..." "..." ;`
    (or `const std::string NAME[] = { "...", "...", "EOS" };` for angle
    data, which RDKit splits across 5 array entries to dodge an old MSVC
    string-literal-length limit) into one decoded Python string."""
    if is_array:
        pat = re.compile(
            r"const std::string " + re.escape(var_name) + r"\[\]\s*=\s*\{(.*?)\};",
            re.S,
        )
    else:
        pat = re.compile(
            r"const std::string " + re.escape(var_name) + r"\s*=\s*(.*?);", re.S
        )
    m = pat.search(text)
    if not m:
        raise ValueError(f"block {var_name} not found")
    body = m.group(1)
    lit_pat = re.compile(r'"((?:[^"\\]|\\.)*)"')
    raw = "".join(lit_pat.findall(body))
    return raw.replace("\\t", "\t").replace("\\n", "\n").replace('\\"', '"')


def parse_table(raw_text, ncols):
    rows = []
    for line in raw_text.split("\n"):
        if not line or line.startswith("*"):
            continue
        cols = [c for c in line.split("\t") if c != ""]
        if len(cols) >= ncols:
            rows.append(cols)
    return rows


def load_rdkit_tables():
    # `curl`, not `urllib`, matching PROVENANCE.md's own documented
    # regeneration commands -- more reliable than Python's urllib against
    # environments with an incomplete local CA bundle.
    cpp = subprocess.run(
        ["curl", "-sL", PARAMS_CPP_URL], check=True, capture_output=True, text=True
    ).stdout

    eqlevel = {}
    for cols in parse_table(extract_cpp_string_block(cpp, "defaultMMFFDef"), 6):
        try:
            atom_type = int(cols[1])
            levels = [int(cols[2]), int(cols[3]), int(cols[4]), int(cols[5])]
        except (ValueError, IndexError):
            continue
        eqlevel.setdefault(atom_type, levels)

    prop = {}
    for cols in parse_table(extract_cpp_string_block(cpp, "defaultMMFFProp"), 9):
        try:
            prop[int(cols[0])] = dict(
                atno=int(cols[1]),
                crd=int(cols[2]),
                val=int(cols[3]),
                pilp=int(cols[4]),
                mltb=int(cols[5]),
                arom=int(cols[6]),
                linh=int(cols[7]),
                sbmb=int(cols[8]),
            )
        except (ValueError, IndexError):
            continue

    bond = {}
    for cols in parse_table(extract_cpp_string_block(cpp, "defaultMMFFBond"), 5):
        try:
            bt, i, j = int(cols[0]), int(cols[1]), int(cols[2])
            kb, r0 = float(cols[3]), float(cols[4])
        except (ValueError, IndexError):
            continue
        lo, hi = (i, j) if i <= j else (j, i)
        bond[(bt, lo, hi)] = (kb, r0)

    bndk = {}
    for cols in parse_table(extract_cpp_string_block(cpp, "defaultMMFFBndk"), 4):
        try:
            i, j = int(cols[0]), int(cols[1])
            r0, kb = float(cols[2]), float(cols[3])
        except (ValueError, IndexError):
            continue
        lo, hi = (i, j) if i <= j else (j, i)
        bndk[(lo, hi)] = (r0, kb)

    hl = {}
    for cols in parse_table(
        extract_cpp_string_block(cpp, "defaultMMFFHerschbachLaurie"), 5
    ):
        try:
            i, j = int(cols[0]), int(cols[1])
            a_ij, d_ij, dp_ij = float(cols[2]), float(cols[3]), float(cols[4])
        except (ValueError, IndexError):
            continue
        lo, hi = (i, j) if i <= j else (j, i)
        hl[(lo, hi)] = (a_ij, d_ij, dp_ij)

    covrad = {}
    for cols in parse_table(
        extract_cpp_string_block(cpp, "defaultMMFFCovRadPauEle"), 3
    ):
        try:
            covrad[int(cols[0])] = (float(cols[1]), float(cols[2]))
        except (ValueError, IndexError):
            continue

    angle = {}
    for cols in parse_table(
        extract_cpp_string_block(cpp, "defaultMMFFAngleData", is_array=True), 6
    ):
        if cols[0] == "EOS":
            continue
        try:
            at, i, j, k = int(cols[0]), int(cols[1]), int(cols[2]), int(cols[3])
            ka, theta0 = float(cols[4]), float(cols[5])
        except (ValueError, IndexError):
            continue
        lo, hi = (i, k) if i <= k else (k, i)
        angle[(at, lo, j, hi)] = (ka, theta0)

    return eqlevel, prop, bond, bndk, hl, covrad, angle


EQLEVEL, PROP, BOND, BNDK, HL, COVRAD, ANGLE = load_rdkit_tables()
print(
    f"parsed: eqLevel={len(EQLEVEL)} prop={len(PROP)} bond={len(BOND)} "
    f"bndk={len(BNDK)} hl={len(HL)} covrad={len(COVRAD)} angle={len(ANGLE)}",
    file=sys.stderr,
)

# ---------------------------------------------------------------------------
# RDKit's REAL lookup / empirical-rule algorithms, ported
# ---------------------------------------------------------------------------


def angle_table_lookup(angle_type, ti, tj, tk):
    """Mirrors MMFFAngleCollection::operator() exactly: angleType and tj
    fixed, ti/tk canonical-substituted via eqLevel stages 0..3 (RDKit's
    Level 2/3/4/5; Level 1 is skipped, identical to Level 2)."""
    for stage in range(4):
        ci = EQLEVEL.get(ti, [ti] * 4)[stage]
        ck = EQLEVEL.get(tk, [tk] * 4)[stage]
        lo, hi = (ci, ck) if ci <= ck else (ck, ci)
        hit = ANGLE.get((angle_type, lo, tj, hi))
        if hit is not None:
            return stage, hit
    return None, None


def bond_table_lookup(bond_type, ti, tj):
    lo, hi = (ti, tj) if ti <= tj else (tj, ti)
    return BOND.get((bond_type, lo, hi))


def periodic_row_hl(atomic_num):
    row = 0
    if atomic_num == 2:
        row = 1
    elif 3 <= atomic_num <= 10:
        row = 2
    elif 11 <= atomic_num <= 18:
        row = 3
    elif 19 <= atomic_num <= 36:
        row = 4
    elif 37 <= atomic_num <= 54:
        row = 5
    if (21 <= atomic_num <= 30) or (39 <= atomic_num <= 48):
        row *= 10
    return row


def bond_empirical(atomic_num_i, atomic_num_j):
    """eq 18/19, MMFF.V page 625. Returns (kb, r0) or None if CovRadPauEle
    is missing for either element (RDKit would throw a PRECONDITION;
    treated here as 'unsupported')."""
    if atomic_num_i not in COVRAD or atomic_num_j not in COVRAD:
        return None
    r0_i, chi_i = COVRAD[atomic_num_i]
    r0_j, chi_j = COVRAD[atomic_num_j]
    c = 0.050 if (atomic_num_i == 1 or atomic_num_j == 1) else 0.085
    n = 1.4
    r0 = r0_i + r0_j - c * (abs(chi_i - chi_j) ** n)
    lo, hi = (
        (atomic_num_i, atomic_num_j)
        if atomic_num_i <= atomic_num_j
        else (atomic_num_j, atomic_num_i)
    )
    if (lo, hi) in BNDK:
        bndk_r0, bndk_kb = BNDK[(lo, hi)]
        kb = bndk_kb * ((bndk_r0 / r0) ** 6)
    else:
        rowlo, rowhi = periodic_row_hl(atomic_num_i), periodic_row_hl(atomic_num_j)
        rlo, rhi = (rowlo, rowhi) if rowlo <= rowhi else (rowhi, rowlo)
        if (rlo, rhi) not in HL:
            return None
        a_ij, d_ij, _ = HL[(rlo, rhi)]
        kb = 10.0 ** (-(r0 - a_ij) / d_ij)
    return kb, r0


Z_TABLE = {
    1: 1.395, 6: 2.494, 7: 2.711, 8: 3.045, 9: 2.847, 14: 2.350, 15: 2.350,
    16: 2.980, 17: 2.909, 35: 3.017, 53: 3.086,
}
C_TABLE = {6: 1.016, 7: 1.113, 8: 1.337, 14: 0.811, 15: 1.068, 16: 1.249, 17: 1.078}


def angle_empirical(
    atomic_num_i, atomic_num_j_central, atomic_num_k,
    prop_central, r0_ij, r0_jk, ring_size, old_theta0=None,
):
    """eq 20, MMFF.V page 628, plus the theta0 empirical rule immediately
    above it (same RDKit function, `getMMFFAngleBendEmpiricalRuleParams`).
    `old_theta0`: RDKit reuses a found-but-ka==0 table row's OWN theta0
    verbatim (`mmffAngleParams->theta0 = oldMMFFAngleParams->theta0`)
    instead of running the theta0 sub-rule at all -- only a genuine "no row
    found at any eqLevel stage" case runs the sub-rule below."""
    if old_theta0 is not None:
        return _angle_empirical_ka(
            atomic_num_i, atomic_num_j_central, atomic_num_k, r0_ij, r0_jk,
            ring_size, old_theta0,
        )
    theta0 = 120.0
    crd = prop_central["crd"]
    if crd == 4:
        theta0 = 109.45
    elif crd == 2:
        if atomic_num_j_central == 8:
            theta0 = 105.0
        elif prop_central["linh"] == 1:
            theta0 = 180.0
    elif crd == 3:
        if prop_central["val"] == 3 and prop_central["mltb"] == 0:
            theta0 = 107.0 if atomic_num_j_central == 7 else 92.0
    if ring_size == 3:
        theta0 = 60.0
    elif ring_size == 4:
        theta0 = 90.0
    return _angle_empirical_ka(
        atomic_num_i, atomic_num_j_central, atomic_num_k, r0_ij, r0_jk,
        ring_size, theta0,
    )


def _angle_empirical_ka(
    atomic_num_i, atomic_num_j_central, atomic_num_k, r0_ij, r0_jk, ring_size, theta0
):
    beta = 1.75
    if ring_size == 4:
        beta *= 0.85
    elif ring_size == 3:
        beta *= 0.05
    zi, zk = Z_TABLE.get(atomic_num_i, 0.0), Z_TABLE.get(atomic_num_k, 0.0)
    cj = C_TABLE.get(atomic_num_j_central, 0.0)
    d = (r0_ij - r0_jk) ** 2 / (r0_ij + r0_jk) ** 2
    theta0_rad = math.radians(theta0)
    ka = beta * zi * cj * zk / ((r0_ij + r0_jk) * theta0_rad * theta0_rad * math.exp(2.0 * d))
    return ka, theta0


def angle_ring_size(mol, i1, i2, i3):
    """isAngleInRingOfSize3or4 (`AtomTyper.cpp`): local bond adjacency, NOT
    SSSR-based."""
    if mol.GetBondBetweenAtoms(i1, i2) is None or mol.GetBondBetweenAtoms(i2, i3) is None:
        return 0
    if mol.GetBondBetweenAtoms(i3, i1) is not None:
        return 3
    s1 = {n.GetIdx() for n in mol.GetAtomWithIdx(i1).GetNeighbors() if n.GetIdx() != i2}
    s2 = {n.GetIdx() for n in mol.GetAtomWithIdx(i3).GetNeighbors() if n.GetIdx() != i2}
    return 4 if (s1 & s2) else 0


# ---------------------------------------------------------------------------
# Load chematic's own table_gap rows and dedupe to unique tuples
# ---------------------------------------------------------------------------


def load_unique_table_gap_rows(audit_jsonl_path, term_kind):
    seen = {}
    with open(audit_jsonl_path) as f:
        for line in f:
            row = json.loads(line)
            if row.get("term_kind") != term_kind:
                continue
            if row.get("present_at_different_classification") is not None:
                continue  # routing-bug candidate, not a genuine table_gap
            key = tuple(row["lookup_key_before_normalization"])
            seen.setdefault(key, row)
    return list(seen.values())


def classify_angle(row):
    at, cti, ctj, ctk = row["lookup_key_before_normalization"]
    smiles = row["smiles"]
    i1, i2, i3 = [a["index"] for a in row["atoms"]]

    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return dict(**row, verdict="parse_failed_in_rdkit")
    props = AllChem.MMFFGetMoleculeProperties(mol)
    if props is None:
        return dict(**row, verdict="rdkit_mmff_setup_failed")

    rti, rtj, rtk = (
        props.GetMMFFAtomType(i1),
        props.GetMMFFAtomType(i2),
        props.GetMMFFAtomType(i3),
    )
    if (rti, rtj, rtk) != (cti, ctj, ctk):
        return dict(
            **row,
            verdict="type_mismatch",
            rdkit_types=[rti, rtj, rtk],
            chematic_types=[cti, ctj, ctk],
            any_aromatic_atom_involved=any(
                mol.GetAtomWithIdx(i).GetIsAromatic() for i in (i1, i2, i3)
            ),
        )

    stage, table_hit = angle_table_lookup(at, cti, ctj, ctk)
    live = props.GetMMFFAngleBendParams(mol, i1, i2, i3)

    if table_hit is not None and table_hit[0] != 0.0:
        return dict(
            **row,
            verdict="direct_table" if stage == 0 else "equivalence_table",
            eqlevel_stage=stage,
            table_ka_theta0=list(table_hit),
            live_oracle=list(live) if live else None,
        )

    # Table miss (or ka==0.0, RDKit's own isDoubleZero empirical trigger) ->
    # empirical rule. Need both flanking bonds' (kb, r0) via the SAME public
    # API RDKit's own getMMFFAngleBendParams calls internally.
    an_i, an_j, an_k = (
        mol.GetAtomWithIdx(i1).GetAtomicNum(),
        mol.GetAtomWithIdx(i2).GetAtomicNum(),
        mol.GetAtomWithIdx(i3).GetAtomicNum(),
    )
    r0_bonds = []
    for a, b in ((i1, i2), (i2, i3)):
        bres = props.GetMMFFBondStretchParams(mol, a, b)
        if bres is None:
            return dict(**row, verdict="unsupported_bond_empirical_failed")
        r0_bonds.append(bres[2])

    prop_central = PROP.get(ctj)
    if prop_central is None:
        return dict(**row, verdict="unsupported_no_prop_row_for_central_type")

    ring_size = angle_ring_size(mol, i1, i2, i3)
    old_theta0 = table_hit[1] if table_hit is not None else None
    ka_pred, theta0_pred = angle_empirical(
        an_i, an_j, an_k, prop_central, r0_bonds[0], r0_bonds[1], ring_size,
        old_theta0=old_theta0,
    )

    matches_live = None
    if live is not None:
        _, ka_live, theta0_live = live
        matches_live = abs(ka_live - ka_pred) < 1e-2 and abs(theta0_live - theta0_pred) < 1e-1

    return dict(
        **row,
        verdict="empirical_rule",
        empirical_trigger="zero_ka_table_row" if table_hit is not None else "no_table_row",
        ring_size=ring_size,
        predicted_ka_theta0=[round(ka_pred, 4), round(theta0_pred, 4)],
        live_oracle=list(live) if live else None,
        predicted_matches_live_oracle=matches_live,
    )


def classify_bond(row):
    bt, cti, ctj = row["lookup_key_before_normalization"]
    smiles = row["smiles"]
    i1, i2 = [a["index"] for a in row["atoms"]]

    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return dict(**row, verdict="parse_failed_in_rdkit")
    props = AllChem.MMFFGetMoleculeProperties(mol)
    if props is None:
        return dict(**row, verdict="rdkit_mmff_setup_failed")

    rti, rtj = props.GetMMFFAtomType(i1), props.GetMMFFAtomType(i2)
    if (rti, rtj) != (cti, ctj):
        return dict(
            **row,
            verdict="type_mismatch",
            rdkit_types=[rti, rtj],
            chematic_types=[cti, ctj],
            any_aromatic_atom_involved=any(
                mol.GetAtomWithIdx(i).GetIsAromatic() for i in (i1, i2)
            ),
        )

    table_hit = bond_table_lookup(bt, cti, ctj)
    live = props.GetMMFFBondStretchParams(mol, i1, i2)
    if table_hit is not None:
        return dict(
            **row, verdict="direct_table", table_kb_r0=list(table_hit),
            live_oracle=list(live) if live else None,
        )

    an_i, an_j = mol.GetAtomWithIdx(i1).GetAtomicNum(), mol.GetAtomWithIdx(i2).GetAtomicNum()
    emp = bond_empirical(an_i, an_j)
    if emp is None:
        return dict(**row, verdict="unsupported_empirical_failed")
    kb_pred, r0_pred = emp
    matches_live = None
    if live is not None:
        _, kb_live, r0_live = live
        matches_live = abs(kb_live - kb_pred) < 1e-2 and abs(r0_live - r0_pred) < 1e-3
    return dict(
        **row,
        verdict="empirical_rule",
        predicted_kb_r0=[round(kb_pred, 4), round(r0_pred, 4)],
        live_oracle=list(live) if live else None,
        predicted_matches_live_oracle=matches_live,
    )


def main():
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} <mmff94_term_coverage_audit.jsonl>", file=sys.stderr)
        sys.exit(1)
    audit_path = sys.argv[1]

    angle_rows = load_unique_table_gap_rows(audit_path, "Angle")
    bond_rows = load_unique_table_gap_rows(audit_path, "Bond")

    angle_results = [classify_angle(r) for r in angle_rows]
    bond_results = [classify_bond(r) for r in bond_rows]

    print(json.dumps({"angle": angle_results, "bond": bond_results}))

    print(f"\n=== ANGLE ({len(angle_rows)} unique tuples) ===", file=sys.stderr)
    c = collections.Counter(r["verdict"] for r in angle_results)
    for k, v in c.most_common():
        print(f"  {k}: {v}", file=sys.stderr)
    n_emp = sum(1 for r in angle_results if r["verdict"] == "empirical_rule")
    n_emp_matches = sum(
        1 for r in angle_results
        if r["verdict"] == "empirical_rule" and r.get("predicted_matches_live_oracle")
    )
    print(
        f"  -> empirical_rule predicted-value matches live oracle: {n_emp_matches}/{n_emp}",
        file=sys.stderr,
    )

    print(f"\n=== BOND ({len(bond_rows)} unique tuples) ===", file=sys.stderr)
    c2 = collections.Counter(r["verdict"] for r in bond_results)
    for k, v in c2.most_common():
        print(f"  {k}: {v}", file=sys.stderr)


if __name__ == "__main__":
    main()
