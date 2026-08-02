#!/usr/bin/env python3
"""Issue #227 Priority 1A-2, Phase 1.1: exhaustive per-atom evidence for the
79 atoms remaining in the Priority-1A residual after this PR's two
production fixes (NC=C, 129 atoms; the C=C/C=O double-bond-partner
umbrella, 39 atoms -- both closed to 0, verified by corpus re-measurement
and Rust regression tests, and therefore not re-documented here on a
per-atom basis). Classified into exclusive buckets, 0 unclassified.

Ground truth: `validation/results/mmff94_chematic_numeric_types.jsonl`
(post-fix chematic dump) joined against
`validation/results/mmff94_rdkit_type_oracle.jsonl` by molecule + atom
index, plus per-atom structural features from the chematic Python bindings.

Run: .venv/bin/python scripts/mmff94_atom_type_bucket_evidence_227.py \
  > validation/results/mmff94_atom_type_bucket_evidence_227.jsonl
"""

import json
import re
import sys

import chematic

REGISTRY_PATH = "crates/chematic-ff/src/mmff94_numeric_type_registry.rs"
RDKIT_ORACLE_PATH = "validation/results/mmff94_rdkit_type_oracle.jsonl"
CHEMATIC_DUMP_PATH = "validation/results/mmff94_chematic_numeric_types.jsonl"

# Bucket classification per (chematic_symbol, rdkit_symbol) group for the
# 79-atom residual remaining after this PR's two fixes. See this PR's body
# for the full evidence trail (SSSR comparisons, minimal-fragment RDKit
# probes, AtomTyper.cpp line citations) behind each entry.
BUCKET_INFO = {
    ("O=C", "O2CM"): dict(
        bucket="terminal_oxygen_o2cm_umbrella_gap",
        responsible_rdkit_rule="AtomTyper.cpp setMMFFHeavyAtomType case 8 (aliphatic O), "
        "1-neighbor branch, lines ~1600-1720 at the pinned commit "
        "(isCarboxylateO / isNitroO / isNOxideO / isThioSulfinateO / isSulfateO / "
        "isPhosphateOrPerchlorateO union -> type 32)",
        proposed_fix="not implemented -- requires porting the full sulfate/sulfone/"
        "phosphate/perchlorate/carboxylate/nitro terminal-O disambiguation "
        "(a materially larger, separate port than this PR's two fixes); "
        "spot-checked examples are sulfonyl (S(=O)(=O)) terminal oxygens, "
        "not literal carboxylates, despite the registry's descriptive name",
        blast_radius=35,
    ),
    ("CB", "C=C"): dict(
        bucket="charged_macrocycle_sssr_divergence",
        responsible_rdkit_rule="not a typing rule -- upstream ring-perception divergence "
        "feeding compute_mmff94_aromatic_view (Aromaticity.cpp setMMFFAromaticity)",
        proposed_fix="not implemented -- confirmed via direct SSSR comparison "
        "(chematic_perception::find_sssr vs RDKit RingInfo) on the responsible "
        "molecule: both engines agree on the two 6-membered pyridinium/benzo "
        "rings themselves, but chematic's SSSR finds 1 macrocyclic ring (28 atoms) "
        "for the whole bis-pyridinium-macrocycle topology where RDKit's SSSR finds "
        "2 -- a tie-break difference in minimum-cycle-basis selection for a "
        "degenerate/symmetric macrocycle, which then perturbs "
        "compute_mmff94_aromatic_view's cross-ring convergence loop for the "
        "otherwise-identically-perceived small ring. NOT an exocyclic-multiple-bond "
        "override (the originally-suspected mechanism from PR #238's body) -- "
        "the relevant ring atoms carry no exocyclic multiple bond; this is a "
        "correction to that earlier characterization. Root-causing the SSSR "
        "tie-break itself is out of this PR's scope (a chematic-perception "
        "core-algorithm change, not an mmff94_numeric.rs typing fix). Two isolated "
        "minimal-fragment RDKit probes (2-methylisoquinolinium, with and without an "
        "aniline substituent) both came back fully aromatic, confirming the "
        "divergence needs this molecule's full macrocyclic topology to reproduce, "
        "not just the local charged-ring pattern.",
        blast_radius=16,
    ),
    ("NPD+", "N+=C"): dict(
        bucket="charged_macrocycle_sssr_divergence",
        responsible_rdkit_rule="same as CB-vs-C=C above -- same 6 molecules "
        "(chembl_tier_b_0009/0023/0028/0029/0030/0034), same SSSR-divergence "
        "mechanism, confirmed via direct comparison",
        proposed_fix="same as CB-vs-C=C above",
        blast_radius=8,
    ),
    ("CB", "C=O"): dict(
        bucket="charged_macrocycle_sssr_divergence",
        responsible_rdkit_rule="same as CB-vs-C=C above",
        proposed_fix="same as CB-vs-C=C above",
        blast_radius=8,
    ),
    ("N=C", "NSP"): dict(
        bucket="nitrile_n_approximation",
        responsible_rdkit_rule="AtomTyper.cpp case 7, 1-neighbor branch (isNSP, "
        "triple-bonded N), type 42",
        proposed_fix="not implemented -- already explicitly flagged as a known "
        "approximation in assign_n_type's own pre-existing comment "
        "(\"close approximation for nitrile\"), not a new finding",
        blast_radius=2,
    ),
    ("NR", "NSO2"): dict(
        bucket="nso2_sulfonamide_cyano_n_unimplemented",
        responsible_rdkit_rule="AtomTyper.cpp case 7, isNSO2orNSO3orNCN "
        "(ipso N bonded to P/S with >=2 terminal O, or to cyano-C), type 43",
        proposed_fix="not implemented -- deliberately excluded from this PR's "
        "NC=C port scope (see classify_n_c3_carbon_context's doc: "
        "is_cyano_like gate falls through rather than mis-typing as 40)",
        blast_radius=2,
    ),
    ("C=O", "CSP"): dict(
        bucket="heterocumulene_degree2_gate_missing",
        responsible_rdkit_rule="AtomTyper.cpp case 6, 2-neighbor branch "
        "(CSP/allenic carbon, type 4) -- takes priority over the 3-neighbor "
        "double-bond-partner check this PR's fix extended",
        proposed_fix="not implemented -- chematic's assign_c_type checks "
        "double_bonds>0 before any degree gate, so a cumulated heterocumulene "
        "carbon (e.g. N=C=S isothiocyanate) with degree 2 falls into the "
        "3-neighbor-shaped double-bond-partner branch instead of a dedicated "
        "degree==2 check; pre-existing, present before and after this PR's "
        "fixes (unchanged count), not introduced or worsened by them",
        blast_radius=2,
    ),
    ("OM", "O2CM"): dict(
        bucket="terminal_oxygen_o2cm_umbrella_gap",
        responsible_rdkit_rule="same umbrella as O=C-vs-O2CM above, reached via a "
        "different chematic assign_o_type branch (anionic-oxygen path)",
        proposed_fix="same as O=C-vs-O2CM above",
        blast_radius=2,
    ),
    ("NR+", "=N="): dict(
        bucket="azide_diazo_typing_not_ported",
        responsible_rdkit_rule="AtomTyper.cpp case 7, 2-neighbor branch "
        "(central N of C=N=N/N=N=N, type 53)",
        proposed_fix="not implemented -- chematic's charge>0 short-circuit in "
        "assign_n_type fires before any azide-specific structural check "
        "could run",
        blast_radius=1,
    ),
    ("N=C", "NAZT"): dict(
        bucket="azide_diazo_typing_not_ported",
        responsible_rdkit_rule="AtomTyper.cpp case 7, 1-neighbor branch "
        "(terminal N of azido/diazo group, type 47)",
        proposed_fix="same as NR+-vs-=N= above -- same molecule, same azide group",
        blast_radius=1,
    ),
    ("NR+", "NO2"): dict(
        bucket="charge_shortcircuit_masks_nitro_n",
        responsible_rdkit_rule="AtomTyper.cpp case 7, 3-neighbor branch, "
        "nTermObondedToN>=2 -> type 45 (nitro N)",
        proposed_fix="not implemented -- assign_n_type's `atom.charge > 0` "
        "check returns type 34 unconditionally before the existing nitro-N "
        "check (double_o>=2) ever runs; a formally-neutral-by-resonance "
        "nitro N written with a +1 formal charge on N (paired with -1 on one "
        "O) is caught by this early short-circuit",
        blast_radius=1,
    ),
    ("S", "S=O"): dict(
        bucket="charged_sulfoxide_s_unhandled",
        responsible_rdkit_rule="AtomTyper.cpp case 16 (sulfur), not read in full "
        "this session -- structural signature is a charged sulfoxide "
        "([S+]([O-])) in a beta-lactam/cephalosporin-like ring",
        proposed_fix="not implemented -- assign_s_type does not yet special-case "
        "charged sulfoxide sulfur",
        blast_radius=1,
    ),
}


def load_registry_symbols():
    text = open(REGISTRY_PATH).read()
    symbols = {}
    for m in re.finditer(r"id:\s*(\d+),\s*symbol:\s*\"([^\"]+)\"", text):
        symbols[int(m.group(1))] = m.group(2)
    return symbols


def main():
    symbols = load_registry_symbols()

    rdkit_rows = {}
    with open(RDKIT_ORACLE_PATH) as f:
        for line in f:
            row = json.loads(line)
            if row.get("status") == "ok":
                rdkit_rows[row["name"]] = row

    chematic_rows = {}
    with open(CHEMATIC_DUMP_PATH) as f:
        for line in f:
            row = json.loads(line)
            if row.get("status") == "ok":
                chematic_rows[row["name"]] = row

    rows_out = []
    counts = {}
    mol_cache = {}

    for name, c_row in chematic_rows.items():
        r = rdkit_rows.get(name)
        if r is None:
            continue
        r_by_idx = {a["index"]: a for a in r["atom_types"]}

        for a in c_row["atoms"]:
            idx = a["index"]
            ra = r_by_idx.get(idx)
            if ra is None:
                continue
            chematic_type_id = a["chematic_numeric_type"]
            rdkit_type_id = ra["rdkit_mmff_type"]
            if chematic_type_id == rdkit_type_id:
                continue  # exact match, not a residual atom

            chematic_symbol = symbols.get(chematic_type_id)
            rdkit_symbol = symbols.get(rdkit_type_id)
            key = (chematic_symbol, rdkit_symbol)
            info = BUCKET_INFO.get(key)
            if info is None:
                print(
                    f"UNCLASSIFIED: {name}#{idx} chematic={chematic_symbol}({chematic_type_id}) "
                    f"rdkit={rdkit_symbol}({rdkit_type_id})",
                    file=sys.stderr,
                )
                continue

            if name not in mol_cache:
                mol_cache[name] = chematic.from_smiles(c_row["smiles"])
            mol = mol_cache[name]

            atbl = mol.atom_table[idx]
            symbol_, atomic_num, charge, aromatic, implicit_h, degree, in_ring = atbl
            ring_sizes = mol.ring_sizes_for_atom(idx)

            neighbor_sig = []
            multiple_bonds = []
            for b in mol.bond_table:
                a1, a2, order, barom = b
                if a1 == idx or a2 == idx:
                    other = a2 if a1 == idx else a1
                    other_sym = mol.atom_table[other][0]
                    entry = f"{other_sym}({order}{'[arom]' if barom else ''})"
                    neighbor_sig.append(entry)
                    if order in ("DOUBLE", "TRIPLE"):
                        multiple_bonds.append(entry)

            counts[key] = counts.get(key, 0) + 1
            rows_out.append(
                {
                    "molecule": name,
                    "smiles": c_row["smiles"],
                    "atom_index": idx,
                    "element": symbol_,
                    "charge": charge,
                    "aromatic": aromatic,
                    "degree": degree,
                    "implicit_h": implicit_h,
                    "ring_sizes": sorted(ring_sizes),
                    "neighbor_signature": sorted(neighbor_sig),
                    "multiple_bonds": sorted(multiple_bonds),
                    "chematic_type": f"{chematic_symbol}({chematic_type_id})",
                    "rdkit_type": f"{rdkit_symbol}({rdkit_type_id})",
                    "bucket": info["bucket"],
                    "responsible_rdkit_rule": info["responsible_rdkit_rule"],
                    "proposed_fix": info["proposed_fix"],
                    "blast_radius": info["blast_radius"],
                }
            )

    for row in rows_out:
        print(json.dumps(row))

    total = len(rows_out)
    unclassified = sum(
        1
        for name, c_row in chematic_rows.items()
        for a in c_row["atoms"]
        if (r := rdkit_rows.get(name))
        and (ra := {x["index"]: x for x in r["atom_types"]}.get(a["index"]))
        and a["chematic_numeric_type"] != ra["rdkit_mmff_type"]
        and (
            symbols.get(a["chematic_numeric_type"]),
            symbols.get(ra["rdkit_mmff_type"]),
        )
        not in BUCKET_INFO
    )
    print(
        json.dumps(
            {
                "_summary": True,
                "total_rows": total,
                "unclassified": unclassified,
                "counts": {f"{k[0]}->{k[1]}": v for k, v in counts.items()},
            }
        ),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
