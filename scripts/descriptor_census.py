#!/usr/bin/env python3
"""
Full census of every descriptor value exposed by
crates/chematic-chem/src/descriptors.rs (71 functions, ~190 values), measured
against RDKit where an oracle exists.

Scope note: this script covers ONLY descriptors.rs. Many `mol.descriptors()`
dict keys (QED, SA score, Kappa/Chi/BertzCT/WienerIndex, VSA families,
EState, pKa, ADMET, xlogp3/esol/logd) live in *other* files in the
chematic-chem crate (qed.rs, sa_score.rs, topo_descriptors.rs, vsa.rs,
estate.rs, pka.rs, admet.rs, xlogp3.rs, esol.rs, logd.rs) and are
out of scope here — see docs/rfcs/descriptor_census_rfc.md.

Five descriptors.rs functions have NO Python/WASM/MCP binding at all
(moran_autocorr, geary_autocorr, information_content, mde_carbon, and the
plain mmff94_charges, shadowed by mmff94_charges_bci) — those are read from
validation/results/descriptor_census_unbound.jsonl, produced by:

    cargo run -p chematic-chem --release --example descriptor_census_unbound \
        < scripts/descriptor_census_corpus.smi \
        > validation/results/descriptor_census_unbound.jsonl

Usage:
    .venv/bin/python scripts/descriptor_census.py \
        --corpus scripts/descriptor_census_corpus.smi \
        --unbound validation/results/descriptor_census_unbound.jsonl \
        --json validation/results/descriptor_census.json
"""
from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
from collections import defaultdict

try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors, Crippen
except ImportError:
    sys.exit("rdkit not installed")

try:
    import chematic
except ImportError:
    sys.exit("chematic not installed (see docs/rfcs/descriptor_census_rfc.md for the isolated-venv build steps)")


# ---------------------------------------------------------------------------
# Generic accumulator for one scalar descriptor across the corpus.
# ---------------------------------------------------------------------------
class Acc:
    def __init__(self, name, kind, dependency, oracle, note=""):
        self.name = name
        self.kind = kind  # "float" | "int" | "bool" | "str"
        self.dependency = dependency
        self.oracle = oracle  # human description, or None if no RDKit oracle
        # `note` is a root-cause explanation, NOT a final classification --
        # the actual classification tier (confirms-doc / minor / major /
        # unexplained / no-oracle) is derived in finalize() from the measured
        # exact-match rate, with `note` overlaid when a root cause is known.
        # A static pre-assigned classification would silently keep saying
        # "confirms doc" on rows that measure as heavily divergent.
        self.note = note
        self.fixture_count = 0
        self.valid_count = 0  # both ch and rd produced a usable value
        self.ch_error = 0
        self.rd_error = 0
        self.nan_inf_both = 0
        self.nan_inf_mismatch = 0
        self.exact_match = 0
        self.abs_errors = []
        self.worst = None  # (abs_err, smiles, ch_val, rd_val)

    def record_no_oracle(self, smi, ch_val, ch_err):
        self.fixture_count += 1
        if ch_err:
            self.ch_error += 1
        else:
            self.valid_count += 1
            if isinstance(ch_val, float) and (math.isnan(ch_val) or math.isinf(ch_val)):
                self.nan_inf_both += 1  # just tracked as "non-finite count" here

    def record(self, smi, ch_val, rd_val, ch_err=None, rd_err=None, tol=0.0):
        self.fixture_count += 1
        if ch_err is not None:
            self.ch_error += 1
            return
        if rd_err is not None:
            self.rd_error += 1
            return

        ch_nan = isinstance(ch_val, float) and (math.isnan(ch_val) or math.isinf(ch_val))
        rd_nan = isinstance(rd_val, float) and (math.isnan(rd_val) or math.isinf(rd_val))
        if ch_nan or rd_nan:
            if ch_nan and rd_nan:
                self.nan_inf_both += 1
                self.valid_count += 1
            else:
                self.nan_inf_mismatch += 1
                if self.worst is None or 1e18 > self.worst[0]:
                    pass
                # NaN/Inf disagreement is always "worst" in its own right
                self.worst = self.worst or (float("inf"), smi, ch_val, rd_val)
            return

        self.valid_count += 1

        if self.kind in ("bool", "str"):
            if ch_val == rd_val:
                self.exact_match += 1
            else:
                self.abs_errors.append(1.0)
                if self.worst is None:
                    self.worst = (1.0, smi, ch_val, rd_val)
            return

        # numeric
        err = abs(ch_val - rd_val)
        self.abs_errors.append(err)
        if err <= tol:
            self.exact_match += 1
        if self.worst is None or err > self.worst[0]:
            self.worst = (err, smi, ch_val, rd_val)

    def finalize(self, tol_desc="exact"):
        n = len(self.abs_errors)
        # exact_match_pct's denominator must be every valid (both-sides-computed)
        # fixture, not len(abs_errors) -- for bool/str kinds we only push an
        # entry onto abs_errors on a *mismatch* (abs_errors is only used for
        # MAE/percentile reporting, which is meaningless for bool/str anyway),
        # so len(abs_errors) under-counts the denominator for those kinds.
        match_denom = self.valid_count - self.nan_inf_both
        exact_pct = round(100 * self.exact_match / match_denom, 3) if (self.oracle and match_denom) else None

        # --- classification tier, derived from measured agreement, not pre-assigned ---
        if self.oracle is None:
            tier = "no RDKit oracle"
        elif self.nan_inf_mismatch:
            tier = "major divergence (NaN/Inf disagreement present)"
        elif exact_pct is None:
            tier = "insufficient data"
        elif exact_pct >= 99.0:
            tier = "confirms doc / near-exact"
        elif exact_pct >= 90.0:
            tier = "minor divergence"
        else:
            tier = "major divergence"

        if self.note:
            classification = f"known: {self.note}" if tier not in ("no RDKit oracle",) else f"{tier} -- {self.note}"
        elif tier in ("major divergence", "major divergence (NaN/Inf disagreement present)"):
            classification = f"unexplained -- {tier} (candidate follow-up)"
        else:
            classification = tier

        d = {
            "dependency": self.dependency,
            "kind": self.kind,
            "rdkit_oracle": self.oracle,
            "tolerance": tol_desc,
            "fixture_count": self.fixture_count,
            "valid_count": self.valid_count,
            "chematic_errors": self.ch_error,
            "rdkit_errors": self.rd_error,
            "nan_inf_both": self.nan_inf_both,
            "nan_inf_mismatch": self.nan_inf_mismatch,
            "exact_match_count": self.exact_match if self.oracle else None,
            "exact_match_pct": exact_pct,
            "mae": round(statistics.mean(self.abs_errors), 6) if (self.oracle and n and self.kind not in ("bool", "str")) else None,
            "median_ae": round(statistics.median(self.abs_errors), 6) if (self.oracle and n and self.kind not in ("bool", "str")) else None,
            "p95_ae": round(sorted(self.abs_errors)[int(0.95 * (n - 1))], 6) if (self.oracle and n and self.kind not in ("bool", "str")) else None,
            "max_ae": round(max(self.abs_errors), 6) if (self.oracle and n and self.kind not in ("bool", "str")) else None,
            "worst_fixture": None,
            "classification": classification,
        }
        if self.worst is not None:
            err, smi, chv, rdv = self.worst
            d["worst_fixture"] = {
                "smiles": smi,
                "chematic": chv if not isinstance(chv, float) else round(chv, 6),
                "rdkit": rdv if not isinstance(rdv, float) else round(rdv, 6),
                "abs_error": None if math.isinf(err) else round(err, 6),
            }
        return d


# ---------------------------------------------------------------------------
# RDKit-side helpers for descriptors.rs functions RDKit has no named
# equivalent for, but whose formula is fully documented in the Rust source
# (so we can independently re-derive the "should be" value from RDKit's own
# MW/LogP/TPSA/etc and check chematic's arithmetic, not just borrow its
# inputs).
# ---------------------------------------------------------------------------
def rd_ring_system_count(rm):
    """Count fused-ring clusters (union-find over SSSR rings sharing atoms)."""
    ri = rm.GetRingInfo()
    rings = [set(r) for r in ri.AtomRings()]
    n = len(rings)
    if n == 0:
        return 0
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i in range(n):
        for j in range(i + 1, n):
            if rings[i] & rings[j]:
                union(i, j)
    return len({find(i) for i in range(n)})


def rd_num_ester_bonds(rm):
    # matches chematic's ester SMARTS intent: C(=O)-O-C, excluding carbonate/carbamate
    # via the same simple pattern used across the codebase for "ester bond" counts.
    patt = Chem.MolFromSmarts("[#6](=O)[OX2H0][#6]")
    return len(rm.GetSubstructMatches(patt))


def rd_hybridization_per_atom(rm):
    """Map RDKit hybridization -> chematic's 1=sp,2=sp2,3=sp3 encoding."""
    mapping = {
        Chem.HybridizationType.SP: 1,
        Chem.HybridizationType.SP2: 2,
        Chem.HybridizationType.SP3: 3,
    }
    out = []
    for a in rm.GetAtoms():
        h = a.GetHybridization()
        out.append(mapping.get(h, None))  # None = no comparable RDKit hybridization (SP3D, S, UNSPECIFIED, ...)
    return out


LIPINSKI = lambda mw, hbd, hba, logp: mw <= 500.0 and hbd <= 5 and hba <= 10 and logp <= 5.0
VEBER = lambda tpsa, rotb: tpsa <= 140.0 and rotb <= 10
EGAN = lambda tpsa, logp: tpsa <= 131.6 and logp <= 5.88
REOS = lambda mw, logp, hbd, hba, fc, rotb, hac: (
    200.0 <= mw <= 500.0 and -5.0 <= logp <= 5.0 and 0 <= hbd <= 5 and 0 <= hba <= 10
    and -2 <= fc <= 2 and 0 <= rotb <= 8 and 15 <= hac <= 50
)
GHOSE = lambda mw, logp, hac, mr: (160.0 <= mw <= 480.0 and -0.4 <= logp <= 5.6 and 20.0 <= hac <= 70.0 and 40.0 <= mr <= 130.0)
RO3 = lambda mw, logp, hbd, hba, rotb: mw <= 300.0 and logp <= 3.0 and hbd <= 3 and hba <= 3 and rotb <= 3
LEAD_LIKE = lambda mw, logp, rotb, rings: mw <= 450.0 and -3.5 <= logp <= 4.5 and rotb <= 10 and 1 <= rings <= 4
PFIZER = lambda logp, tpsa: not (logp > 3.0 and tpsa < 75.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--corpus", default="scripts/descriptor_census_corpus.smi")
    ap.add_argument("--unbound", default="validation/results/descriptor_census_unbound.jsonl")
    ap.add_argument("--json", default="validation/results/descriptor_census.json")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    with open(args.corpus) as f:
        smiles_list = [l.strip() for l in f if l.strip()]
    if args.limit:
        smiles_list = smiles_list[: args.limit]

    # unbound_rows[i] is assumed aligned with smiles_list[i] -- both read the
    # same corpus file, in the same order, filtering blank lines the same way
    # (see descriptor_census_unbound.rs). The main loop double-checks this
    # with an exact SMILES-string comparison before trusting it (see
    # `unbound_row` below) rather than assuming the alignment silently holds.
    unbound_rows = []
    try:
        with open(args.unbound) as f:
            for line in f:
                line = line.strip()
                if line:
                    unbound_rows.append(json.loads(line))
    except FileNotFoundError:
        print(f"[warn] {args.unbound} not found; unbound-function rows will be skipped", file=sys.stderr)

    accs = {}

    def A(name, kind, dependency, oracle, note=""):
        accs[name] = Acc(name, kind, dependency, oracle, note)
        return accs[name]

    HEAVY_INDEX_BUG = (
        "topo_dist_usize()/topological_distance_matrix() returns a matrix indexed by "
        "HEAVY-ATOM-COMPACTED position, but this function treats that position as a raw "
        "AtomIdx when looking up per-atom properties (atomic_valence/degree). Confirmed via "
        "a direct reproducer: '[2H]C([2H])([2H])NC=O' panics with 'index out of bounds: the "
        "len is 4 but the index is 4' in moran_autocorr/geary_autocorr (raw 0..n loop, not "
        "iterator-bounded). autocorr_2d/ipc use the same pattern but are defensively bounded "
        "by .take(n), so they don't crash -- they silently look up the WRONG atom's property "
        "instead, on any molecule with explicit H/D/T atoms in the graph. Root cause "
        "identified, NOT fixed here (production change, out of scope for diagnosis)."
    )

    # ---- registry: simple scalar descriptors with a 1:1 RDKit oracle ----
    A("molecular_weight", "float", "independent", "Descriptors.MolWt",
      note="99.98% exact overall; the single worst fixture ([2H]C([2H])([2H])NC=O, off by "
           "3.02 Da) has a clean root cause: molecular_weight() (crates/chematic-chem/src/"
           "descriptors.rs, ~line 186) computes mass via avg_mass(atom.element) only and "
           "never reads atom.isotope, so an explicit deuterium/tritium atom is weighed as "
           "ordinary natural-abundance H. exact_mass() (same file) DOES read atom.isotope "
           "correctly -- this molecule is not in exact_mass's own worst-fixture list, "
           "confirming the two functions diverge in isotope-handling, not in general "
           "accuracy. Same trigger molecule as the moran_autocorr/geary_autocorr panic "
           "(unrelated mechanism, coincidental shared fixture).")
    A("exact_mass", "float", "independent", "Descriptors.ExactMolWt")
    A("heavy_atom_count", "int", "independent", "GetNumHeavyAtoms()")
    A("hbd_count", "int", "depends-on-aromaticity", "rdMolDescriptors.CalcNumHBD")
    A("hba_count", "int", "depends-on-aromaticity", "rdMolDescriptors.CalcNumHBA")
    A("rotatable_bond_count", "int", "depends-on-ring-perception", "rdMolDescriptors.CalcNumRotatableBonds")
    A("tpsa", "float", "depends-on-aromaticity", "rdMolDescriptors.CalcTPSA(includeSandP=True)")
    A("logp_crippen", "float", "independent", "Crippen.MolLogP")
    A("lipinski_passes", "bool", "independent", "derived: MW<=500 and HBD<=5 and HBA<=10 and LogP<=5.0 (formula reproduction on RDKit inputs)")
    A("fsp3", "float", "depends-on-aromaticity", "rdMolDescriptors.CalcFractionCSP3")
    A("aromatic_ring_count", "int", "depends-on-aromaticity", "manual: SSSR rings, all atoms GetIsAromatic()")
    A("formal_charge_sum", "int", "independent", "Chem.GetFormalCharge")
    A("molar_refractivity", "float", "independent", "Crippen.MolMR")
    A("num_heteroatoms", "int", "independent", "Lipinski.NumHeteroatoms")
    A("ring_count", "int", "depends-on-ring-perception", "rdMolDescriptors.CalcNumRings",
      note="SSSR-basis-cardinality disagreement on large fused/bridged macrocycles -- same class "
           "as the [R1]/[R2] SMARTS divergence already documented in docs/rdkit_compat.md's Known "
           "divergence classes table (~1.7% of such molecules). Worst fixture this round is a "
           "16-atom-ring bis-pyridinium macrocycle (chematic=7, RDKit=10 rings) -- not a fresh bug.")
    A("ring_system_count", "int", "depends-on-ring-perception", "manual: union-find over SSSR rings sharing atoms")
    A("hba_count_lipinski", "int", "depends-on-aromaticity", "rdMolDescriptors.CalcNumLipinskiHBA")
    A("fraction_rotatable_bonds", "float", "depends-on-ring-perception", "derived: RotB / HeavyAtoms")
    A("num_aliphatic_rings", "int", "depends-on-aromaticity", "rdMolDescriptors.CalcNumAliphaticRings",
      note="same SSSR-basis-cardinality mechanism as ring_count above -- concentrated in large "
           "polycyclic macrocycles, not a fresh bug.")
    A("num_saturated_rings", "int", "depends-on-aromaticity", "rdMolDescriptors.CalcNumSaturatedRings")
    A("num_aromatic_heterocycles", "int", "depends-on-aromaticity", "rdMolDescriptors.CalcNumAromaticHeterocycles")
    A("num_aliphatic_heterocycles", "int", "depends-on-aromaticity", "rdMolDescriptors.CalcNumAliphaticHeterocycles",
      note="same SSSR-basis-cardinality mechanism as ring_count above.")
    A("num_saturated_heterocycles", "int", "depends-on-aromaticity", "rdMolDescriptors.CalcNumSaturatedHeterocycles")
    A("num_spiro_atoms", "int", "depends-on-ring-perception", "rdMolDescriptors.CalcNumSpiroAtoms")
    A("num_bridgehead_atoms", "int", "depends-on-ring-perception", "rdMolDescriptors.CalcNumBridgeheadAtoms")
    A("num_stereocenters", "int", "depends-on-CIP-or-stereo", "rdMolDescriptors.CalcNumAtomStereoCenters (legacy oracle)")
    A("num_stereocenters_new_cip", "int", "depends-on-CIP-or-stereo", "FindPotentialStereo (modern CIP oracle)")
    A("num_unspecified_stereocenters", "int", "depends-on-CIP-or-stereo", "rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters",
      note="ROOT CAUSE FOUND (read the source, crates/chematic-chem/src/descriptors.rs "
           "num_unspecified_stereocenters(), ~line 1808): unlike num_stereocenters() (which "
           "does a real 2-pass provisional-CIP substituent-distinctness check and measures "
           "99.76% here), this function's test for 'is this atom a stereocenter' is only "
           "'sp3 carbon, non-aromatic, no explicit chirality tag, degree+implicitH==4, no "
           "attached double/triple bond' -- it NEVER checks whether the 4 substituents are "
           "actually distinct. Every ordinary -CH2- or -CH3 group (two identical H "
           "substituents, can never be a real stereocenter) passes this test and gets counted "
           "as an 'unspecified stereocenter'. This is corpus-wide, not a macrocycle edge case: "
           "worst fixture is a plain polypeptide chain (lots of NC(=O)C backbone -CH2-/-CH< "
           "carbons) where every REAL stereocenter already has @/@@ (RDKit correctly reports "
           "0 unspecified) but chematic reports 48 -- almost certainly counting ordinary "
           "backbone carbons. FLAG FOR STEREO/CIP SPECIALIST to fix (production change, out "
           "of scope for diagnosis): the fix is to reuse num_stereocenters()'s substituent-"
           "distinctness check and only then test for a missing chirality tag.")
    A("veber_passes", "bool", "depends-on-ring-perception", "derived: TPSA<=140 and RotB<=10")
    A("egan_passes", "bool", "depends-on-aromaticity", "derived: TPSA<=131.6 and LogP<=5.88")
    A("reos_passes", "bool", "depends-on-ring-perception", "derived: 6-part range filter (see rust doc)")
    A("ghose_passes", "bool", "independent", "derived: 4-part range filter (see rust doc)")
    A("ro3_passes", "bool", "depends-on-ring-perception", "derived: 5-part range filter (see rust doc)")
    A("lead_like_passes", "bool", "depends-on-ring-perception", "derived: 4-part range filter (see rust doc)")
    A("pfizer_3_75_passes", "bool", "depends-on-aromaticity", "derived: NOT(LogP>3 and TPSA<75)")
    A("cns_mpo_score", "float", "depends-on-pKa(out-of-scope pka.rs)+LogD(out-of-scope logd.rs)", None,
      note="no independent oracle this round -- formula depends on pka.rs (pKa) and logd.rs "
           "(ionisation class for cLogD), both out of this file-scoped census. Needs a joint "
           "census with those files to check independently.")
    A("mcf_passes", "bool", "depends-on-alerts(out-of-scope alerts.rs PAINS/Brenk)", None,
      note="no independent oracle this round -- composes pains_passes/brenk_passes from "
           "alerts.rs, out of this file-scoped census.")
    A("num_carbons", "int", "independent", "count atomic_number==6")
    A("num_nitrogens", "int", "independent", "count atomic_number==7")
    A("num_oxygens", "int", "independent", "count atomic_number==8")
    A("num_fluorines", "int", "independent", "count atomic_number==9")
    A("num_chlorines", "int", "independent", "count atomic_number==17")
    A("num_bromines", "int", "independent", "count atomic_number==35")
    A("num_iodines", "int", "independent", "count atomic_number==53")
    A("num_sulfurs", "int", "independent", "count atomic_number==16")
    A("num_phosphorus", "int", "independent", "count atomic_number==15")
    A("num_hydrogens", "int", "independent", "sum(GetTotalNumHs()) + explicit H atoms",
      note="ROOT CAUSE CONFIRMED with a minimal 4-atom reproducer (not just a hypothesis): "
           "num_hydrogens() (crates/chematic-chem/src/descriptors.rs ~line 2852) computes "
           "atom.hydrogen_count.unwrap_or(0) [explicit, from bracket notation] PLUS "
           "implicit_hcount(mol, idx) for every atom and sums both. For 'C[C@H](N)O', "
           "implicit_hcount_per_atom() is [3,1,2,1] (sums to 7, exactly matching RDKit's "
           "GetTotalNumHs() sum) -- but num_hydrogens() returns 8. implicit_hcount() does not "
           "know the [C@H] atom's H is already counted via its explicit hydrogen_count, so "
           "that one atom's H is counted TWICE. This double-counts by exactly the number of "
           "bracket atoms with an explicit H count whose implicit_hcount() independently "
           "returns >0 -- i.e. any molecule with @/@@ stereocenters written in bracket-H form, "
           "extremely common in real SMILES. FLAG FOR FOLLOW-UP to fix (production change, "
           "out of scope here): num_hydrogens() should not add implicit_hcount() on top of an "
           "already-explicit hydrogen_count.")
    A("num_amide_bonds", "int", "independent", "rdMolDescriptors.CalcNumAmideBonds")
    A("num_ester_bonds", "int", "independent", "manual SMARTS [#6](=O)[OX2H0][#6]")
    A("calc_mol_formula", "str", "independent", "rdMolDescriptors.CalcMolFormula",
      note="chematic's formula omits the ionic-charge suffix RDKit appends for charged species "
           "(e.g. RDKit 'C42H46N4+2' vs chematic 'C42H46N4' for a bis-pyridinium dication) -- "
           "atom/element counts themselves are correct, only the charge annotation is missing.")
    A("balaban_j", "float", "depends-on-ring-perception", "Descriptors.BalabanJ")
    A("ipc", "float", "depends-on-ring-perception", "GraphDescriptors.Ipc",
      note="NAME COLLISION, not a numeric bug: chematic's ipc() computes "
           "Sum(deg_i*deg_j/d(i,j)^2) over atom pairs (docstring: 'Information Path Count'), "
           "which is NOT RDKit's Ipc (Bonchev-Trinajstic total information content on distance-"
           "degeneracy classes) despite the shared name/abbreviation. The two formulas are "
           "unrelated; comparing them produces the trillion-scale MAE below by construction. "
           "Also shares the heavy-atom-index lookup pattern -- see autocorr_2d/moran/geary note.")
    A("hall_kier_alpha", "float", "depends-on-ring-perception", "rdMolDescriptors.CalcHallKierAlpha",
      note="matches docs/rdkit-comparison.md's existing disclosure that Kappa/HallKierAlpha/"
           "BertzCT/BalabanJ/BCUT2D/VSA/MQN/SAScore 'were found to diverge substantially once "
           "measured at corpus scale' -- re-measured here, confirmed for hall_kier_alpha "
           "specifically (BalabanJ, also in that same doc sentence, does NOT diverge in this "
           "round's numbers -- see balaban_j row; the doc's blanket phrasing over-generalizes).")

    # array-family: bcut2d (8), mqn (42), carbon_types (8)
    BCUT2D_NOTE = ("matches docs/rdkit-comparison.md's existing disclosure of substantial "
                   "corpus-scale divergence for BCUT2D; re-measured and confirmed here.")
    BCUT2D_NAMES = ["bcut2d_mwhi", "bcut2d_mwlo", "bcut2d_chghi", "bcut2d_chglo",
                    "bcut2d_logphi", "bcut2d_logplo", "bcut2d_mrhi", "bcut2d_mrlo"]
    for n in BCUT2D_NAMES:
        A(n, "float", "independent", "rdMolDescriptors.BCUT2D", note=BCUT2D_NOTE)
    MQN_NOTE = ("matches docs/rdkit-comparison.md's existing disclosure of substantial "
                "corpus-scale divergence for MQN; re-measured and confirmed here.")
    MQN_NAMES = [f"MQN{i}" for i in range(1, 43)]
    for n in MQN_NAMES:
        A(n, "int", "depends-on-ring-perception", "rdMolDescriptors.MQNs_", note=MQN_NOTE)
    CT_NAMES = ["c1sp1", "c2sp1", "c1sp2", "c2sp2", "c3sp2", "c1sp3", "c2sp3", "c3sp3"]
    for n in CT_NAMES:
        A(n, "int", "independent", "manual: RDKit hybridization x heavy-degree (Mordred CarbonTypes)")

    # per-atom families (aggregate across all atom-instances in the corpus)
    A("hybridization_per_atom", "int", "independent", "atom.GetHybridization() (SP/SP2/SP3 mapped to 1/2/3)")
    A("formal_charge_per_atom", "int", "independent", "atom.GetFormalCharge()")
    A("implicit_hcount_per_atom", "int", "independent", "atom.GetTotalNumHs()")
    A("tpsa_per_atom", "float", "depends-on-aromaticity", "rdMolDescriptors._CalcTPSAContribs(includeSandP=True)")
    A("logp_crippen_per_atom", "float", "independent", "rdMolDescriptors._CalcCrippenContribs()[i][0] (see note)",
      note="attribution-convention difference, NOT a numeric divergence: RDKit's "
           "_CalcCrippenContribs() returns contributions for HEAVY atoms only -- each attached "
           "implicit H's own Crippen atom-type contribution is computed separately inside "
           "MolLogP and is NOT included in the per-atom array (confirmed: summing the 3 "
           "returned rows for 'CCO' gives -0.3487, but Crippen.MolLogP('CCO') is -0.0014). "
           "chematic instead folds each atom's attached-H contribution into that heavy atom's "
           "own value, so raw contribs are not directly comparable atom-for-atom. "
           "formal_charge_per_atom and implicit_hcount_per_atom BOTH measure 100% exact on this "
           "same corpus, which rules out atom-order misalignment as the cause. The molecule-level "
           "aggregate (logp_crippen, above) already validates at ~100%, which is the "
           "correct level to trust this at.")
    A("mr_per_atom", "float", "independent", "rdMolDescriptors._CalcCrippenContribs()[i][1] (see note)",
      note="same attribution-convention difference as logp_crippen_per_atom above -- "
           "molar_refractivity (molecule-level aggregate) already validates at ~100%.")

    # autocorr_2d (7), usrcat (42) — no RDKit equivalent (different definitions)
    AUTOCORR_NOTE = ("no RDKit oracle: RDKit's CalcAUTOCORR2D is a 192-value, differently-"
                      "defined descriptor family (confirmed by direct call: len==192), not "
                      "comparable to this 7-value single-property (atomic valence) lag "
                      "autocorrelation. Shares the heavy-atom-index lookup pattern flagged "
                      "for moran_autocorr/geary_autocorr/ipc -- see HEAVY_INDEX_BUG.")
    for i in range(1, 8):
        A(f"autocorr_2d_lag{i}", "float", "independent", None, note=AUTOCORR_NOTE)
    USRCAT_NOTE = ("no RDKit oracle: this is a 2D-topology-only pseudo-USRCAT (distance-matrix "
                    "average scaled by an arbitrary per-slot factor, not real USR moments); "
                    "RDKit's GetUSRCAT requires a real 3D conformer and computes genuine shape "
                    "descriptors -- the two are not comparable even in principle. Candidate "
                    "follow-up: the 36-value shape block looks like placeholder/synthetic-formula "
                    "code (scale = 1.0 + slot/12.0), worth a closer look outside this diagnosis.")
    for i in range(1, 43):
        A(f"usrcat_{i}", "float", "depends-on-3D-embedding", None, note=USRCAT_NOTE)

    # unbound (no python binding at all) — no RDKit oracle available for any of these
    UNREACHABLE_NOTE = "no RDKit oracle (Mordred-only family); also UNREACHABLE via any binding (Python/WASM/MCP -- confirmed by grep, zero references)."
    for i in range(1, 8):
        A(f"moran_autocorr_lag{i}", "float", "independent", None,
          note=f"{UNREACHABLE_NOTE} {HEAVY_INDEX_BUG}")
    for i in range(1, 8):
        A(f"geary_autocorr_lag{i}", "float", "independent", None,
          note=f"{UNREACHABLE_NOTE} {HEAVY_INDEX_BUG}")
    for n in ["ic", "tic", "sic", "bic", "cic"]:
        A(f"information_content_{n}", "float", "independent", None, note=UNREACHABLE_NOTE)
    for i in range(1, 11):
        A(f"mde_carbon_{i}", "float", "independent", None, note=UNREACHABLE_NOTE)
    A("mmff94_charges_plain", "float", "independent", None,
      note="CORRECTION to an earlier assumption in this same diagnosis: descriptors.rs's "
           "mmff94_charges() is NOT a distinct/simplified formula -- reading the actual function "
           "body (not just its stale doc comment, which still describes an old "
           "'electronegativity-weighted + formal charge' formula) shows it is a 1-line "
           "pass-through to mmff94_bci::mmff94_charges_bci(), the EXACT SAME function Python's "
           "mol.mmff94_charges() calls directly. Values are byte-identical to production output. "
           "It is still UNREACHABLE as this specific symbol (no binding calls "
           "descriptors::mmff94_charges by name), but that is a dead-wrapper finding, not a "
           "dead/wrong-formula finding -- doc-comment drift is the only real defect here.")

    # ---- main loop ----
    total = 0
    rd_parse_fail = 0
    ch_parse_fail = 0
    extended_skipped = 0

    for i, smi in enumerate(smiles_list):
        rd_mol = Chem.MolFromSmiles(smi)
        if rd_mol is None:
            rd_parse_fail += 1
            continue
        try:
            ch_mol = chematic.from_smiles(smi)
        except Exception:
            ch_parse_fail += 1
            continue
        total += 1

        # NOTE: deliberately NOT calling ch_mol.descriptors() here. That dict
        # unconditionally computes ~130 out-of-scope values too (QED, SA
        # score, PAINS, drug_score, ...), and on one pathological symmetric
        # macrocycle in this corpus, drug_score()'s PAINS/VF2 substructure
        # match takes several minutes (root-caused via `sample`+bisection --
        # see docs/rfcs/descriptor_census_rfc.md's VF2 performance finding).
        # Every value below instead comes from the individual getter/method
        # that maps 1:1 to the descriptors.rs function under test, which
        # confirmed fast (<1s total) on the same pathological molecule.
        unbound_row = unbound_rows[i] if i < len(unbound_rows) else None
        if unbound_row is not None and unbound_row.get("smiles") != smi:
            unbound_row = None  # index misalignment guard -- see note where unbound_rows is loaded

        # RDKit reference values (computed once per molecule)
        rd_MW = __import__("rdkit.Chem.Descriptors", fromlist=["MolWt"]).MolWt(rd_mol)
        rd_ExactMW = rdMolDescriptors.CalcExactMolWt(rd_mol)
        rd_HAC = rd_mol.GetNumHeavyAtoms()
        rd_HBD = rdMolDescriptors.CalcNumHBD(rd_mol)
        rd_HBA = rdMolDescriptors.CalcNumHBA(rd_mol)
        rd_RotB = rdMolDescriptors.CalcNumRotatableBonds(rd_mol)
        rd_TPSA = rdMolDescriptors.CalcTPSA(rd_mol, includeSandP=True)
        rd_LogP = Crippen.MolLogP(rd_mol)
        rd_MR = Crippen.MolMR(rd_mol)
        rd_Fsp3 = rdMolDescriptors.CalcFractionCSP3(rd_mol)
        rd_ARC = sum(1 for ring in rd_mol.GetRingInfo().AtomRings()
                     if all(rd_mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring))
        rd_FC = Chem.GetFormalCharge(rd_mol)
        rd_NHet = __import__("rdkit.Chem.Lipinski", fromlist=["NumHeteroatoms"]).NumHeteroatoms(rd_mol)
        rd_RingCount = rdMolDescriptors.CalcNumRings(rd_mol)
        rd_RingSys = rd_ring_system_count(rd_mol)
        rd_HBALip = rdMolDescriptors.CalcNumLipinskiHBA(rd_mol)
        rd_FracRotB = (rd_RotB / rd_HAC) if rd_HAC else 0.0
        rd_NAlR = rdMolDescriptors.CalcNumAliphaticRings(rd_mol)
        rd_NSatR = rdMolDescriptors.CalcNumSaturatedRings(rd_mol)
        rd_NAHet = rdMolDescriptors.CalcNumAromaticHeterocycles(rd_mol)
        rd_NAlHet = rdMolDescriptors.CalcNumAliphaticHeterocycles(rd_mol)
        rd_NSatHet = rdMolDescriptors.CalcNumSaturatedHeterocycles(rd_mol)
        rd_NSpiro = rdMolDescriptors.CalcNumSpiroAtoms(rd_mol)
        rd_NBridge = rdMolDescriptors.CalcNumBridgeheadAtoms(rd_mol)
        rd_NSC_legacy = rdMolDescriptors.CalcNumAtomStereoCenters(rd_mol)
        _find_stereo = getattr(Chem, "FindPotentialStereo", None)
        if _find_stereo is not None:
            rd_NSC_new = sum(1 for s in _find_stereo(rd_mol) if str(s.type).endswith("Atom_Tetrahedral"))
        else:
            rd_NSC_new = rd_NSC_legacy
        rd_NUSC = rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(rd_mol)
        rd_NAmide = rdMolDescriptors.CalcNumAmideBonds(rd_mol)
        rd_NEster = rd_num_ester_bonds(rd_mol)
        rd_Formula = rdMolDescriptors.CalcMolFormula(rd_mol)
        rd_BalabanJ = __import__("rdkit.Chem.Descriptors", fromlist=["BalabanJ"]).BalabanJ(rd_mol) if rd_HAC >= 2 else 0.0
        rd_HallKier = rdMolDescriptors.CalcHallKierAlpha(rd_mol)
        try:
            from rdkit.Chem import GraphDescriptors
            rd_Ipc = GraphDescriptors.Ipc(rd_mol)
        except Exception:
            rd_Ipc = None
        elem_counts = defaultdict(int)
        for a in rd_mol.GetAtoms():
            elem_counts[a.GetAtomicNum()] += 1
        rd_numH = sum(a.GetTotalNumHs() for a in rd_mol.GetAtoms()) + elem_counts.get(1, 0)

        # --- extended metrics that can fail on rare elements (Gasteiger/BCUT/etc) ---
        rd_extended_ok = True
        try:
            rd_bcut = rdMolDescriptors.BCUT2D(rd_mol)
            rd_mqn = rdMolDescriptors.MQNs_(rd_mol)
        except Exception:
            rd_extended_ok = False
            extended_skipped += 1

        accs["molecular_weight"].record(smi, ch_mol.mw, rd_MW, tol=0.01)
        accs["exact_mass"].record(smi, ch_mol.exact_mass, rd_ExactMW, tol=0.01)
        accs["heavy_atom_count"].record(smi, ch_mol.heavy_atoms, rd_HAC, tol=0)
        accs["hbd_count"].record(smi, ch_mol.hbd, rd_HBD, tol=0)
        accs["hba_count"].record(smi, ch_mol.hba, rd_HBA, tol=0)
        accs["rotatable_bond_count"].record(smi, ch_mol.rotatable_bonds, rd_RotB, tol=0)
        accs["tpsa"].record(smi, ch_mol.tpsa, rd_TPSA, tol=0.1)
        accs["logp_crippen"].record(smi, ch_mol.logp, rd_LogP, tol=0.01)
        accs["fsp3"].record(smi, ch_mol.fsp3, rd_Fsp3, tol=0.001)
        accs["aromatic_ring_count"].record(smi, ch_mol.aromatic_ring_count, rd_ARC, tol=0)
        accs["formal_charge_sum"].record(smi, ch_mol.formal_charge, rd_FC, tol=0)
        accs["molar_refractivity"].record(smi, ch_mol.molar_refractivity, rd_MR, tol=0.01)
        accs["num_heteroatoms"].record(smi, ch_mol.num_heteroatoms, rd_NHet, tol=0)
        accs["ring_count"].record(smi, ch_mol.ring_count, rd_RingCount, tol=0)
        accs["ring_system_count"].record(smi, ch_mol.ring_system_count, rd_RingSys, tol=0)
        accs["hba_count_lipinski"].record(smi, ch_mol.hba_count_lipinski, rd_HBALip, tol=0)
        accs["fraction_rotatable_bonds"].record(smi, ch_mol.fraction_rotatable_bonds, rd_FracRotB, tol=1e-6)
        accs["num_aliphatic_rings"].record(smi, ch_mol.num_aliphatic_rings, rd_NAlR, tol=0)
        accs["num_saturated_rings"].record(smi, ch_mol.num_saturated_rings, rd_NSatR, tol=0)
        accs["num_aromatic_heterocycles"].record(smi, ch_mol.num_aromatic_heterocycles, rd_NAHet, tol=0)
        accs["num_aliphatic_heterocycles"].record(smi, ch_mol.num_aliphatic_heterocycles, rd_NAlHet, tol=0)
        accs["num_saturated_heterocycles"].record(smi, ch_mol.num_saturated_heterocycles, rd_NSatHet, tol=0)
        accs["num_spiro_atoms"].record(smi, ch_mol.num_spiro_atoms, rd_NSpiro, tol=0)
        accs["num_bridgehead_atoms"].record(smi, ch_mol.num_bridgehead_atoms, rd_NBridge, tol=0)
        accs["num_stereocenters"].record(smi, ch_mol.num_stereocenters, rd_NSC_legacy, tol=0)
        accs["num_stereocenters_new_cip"].record(smi, ch_mol.num_stereocenters, rd_NSC_new, tol=0)
        accs["num_unspecified_stereocenters"].record(smi, ch_mol.num_unspecified_stereocenters, rd_NUSC, tol=0)
        accs["veber_passes"].record(smi, ch_mol.veber_passes, VEBER(rd_TPSA, rd_RotB))
        accs["egan_passes"].record(smi, ch_mol.egan_passes, EGAN(rd_TPSA, rd_LogP))
        accs["reos_passes"].record(smi, ch_mol.reos_passes, REOS(rd_MW, rd_LogP, rd_HBD, rd_HBA, rd_FC, rd_RotB, rd_HAC))
        accs["ghose_passes"].record(smi, ch_mol.ghose_passes, GHOSE(rd_MW, rd_LogP, rd_HAC, rd_MR))
        accs["ro3_passes"].record(smi, ch_mol.ro3_passes, RO3(rd_MW, rd_LogP, rd_HBD, rd_HBA, rd_RotB))
        accs["lead_like_passes"].record(smi, ch_mol.lead_like_passes, LEAD_LIKE(rd_MW, rd_LogP, rd_RotB, rd_RingCount))
        accs["pfizer_3_75_passes"].record(smi, ch_mol.pfizer_3_75_passes, PFIZER(rd_LogP, rd_TPSA))
        accs["lipinski_passes"].record(smi, ch_mol.lipinski_passes, LIPINSKI(rd_MW, rd_HBD, rd_HBA, rd_LogP))
        accs["cns_mpo_score"].record_no_oracle(smi, ch_mol.cns_mpo_score, None)
        accs["mcf_passes"].record_no_oracle(smi, ch_mol.mcf_passes, None)
        accs["num_carbons"].record(smi, ch_mol.num_carbons, elem_counts.get(6, 0), tol=0)
        accs["num_nitrogens"].record(smi, ch_mol.num_nitrogens, elem_counts.get(7, 0), tol=0)
        accs["num_oxygens"].record(smi, ch_mol.num_oxygens, elem_counts.get(8, 0), tol=0)
        accs["num_fluorines"].record(smi, ch_mol.num_fluorines, elem_counts.get(9, 0), tol=0)
        accs["num_chlorines"].record(smi, ch_mol.num_chlorines, elem_counts.get(17, 0), tol=0)
        accs["num_bromines"].record(smi, ch_mol.num_bromines, elem_counts.get(35, 0), tol=0)
        accs["num_iodines"].record(smi, ch_mol.num_iodines, elem_counts.get(53, 0), tol=0)
        accs["num_sulfurs"].record(smi, ch_mol.num_sulfurs, elem_counts.get(16, 0), tol=0)
        accs["num_phosphorus"].record(smi, ch_mol.num_phosphorus, elem_counts.get(15, 0), tol=0)
        accs["num_hydrogens"].record(smi, ch_mol.num_hydrogens, rd_numH, tol=0)
        accs["num_amide_bonds"].record(smi, ch_mol.num_amide_bonds, rd_NAmide, tol=0)
        accs["num_ester_bonds"].record(smi, ch_mol.num_ester_bonds, rd_NEster, tol=0)
        accs["calc_mol_formula"].record(smi, ch_mol.formula, rd_Formula)
        accs["balaban_j"].record(smi, ch_mol.balaban_j, rd_BalabanJ, tol=0.01)
        if rd_Ipc is not None:
            accs["ipc"].record(smi, ch_mol.ipc, rd_Ipc, tol=0.01)
        accs["hall_kier_alpha"].record(smi, ch_mol.hall_kier_alpha, rd_HallKier, tol=0.01)

        if rd_extended_ok and unbound_row is not None:
            for name, rv in zip(BCUT2D_NAMES, rd_bcut):
                accs[name].record(smi, unbound_row[name], rv, tol=0.01)
            ch_mqn = ch_mol.mqn()
            for name, cv, rv in zip(MQN_NAMES, ch_mqn, rd_mqn):
                accs[name].record(smi, int(cv), int(rv), tol=0)

        ct = rd_hybridization_per_atom(rd_mol)
        ct_ct = {n: 0 for n in CT_NAMES}
        deg_map = {1: "c1", 2: "c2", 3: "c3"}
        for atom, hyb in zip(rd_mol.GetAtoms(), ct):
            if atom.GetAtomicNum() != 6 or hyb is None:
                continue
            heavy_deg = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() != 1)
            key = f"{deg_map.get(heavy_deg)}sp{hyb}" if heavy_deg in deg_map else None
            if key in ct_ct:
                ct_ct[key] += 1
        if unbound_row is not None:
            for name in CT_NAMES:
                accs[name].record(smi, unbound_row[name], ct_ct[name], tol=0)

        # per-atom families
        ch_hyb = ch_mol.hybridization_per_atom()
        rd_hyb = ct
        ch_fc_pa = ch_mol.formal_charge_per_atom()
        rd_fc_pa = [a.GetFormalCharge() for a in rd_mol.GetAtoms()]
        ch_ih_pa = ch_mol.implicit_hcount_per_atom()
        rd_ih_pa = [a.GetTotalNumHs() for a in rd_mol.GetAtoms()]
        ch_tpsa_pa = ch_mol.tpsa_per_atom()
        rd_tpsa_pa = rdMolDescriptors._CalcTPSAContribs(rd_mol, False, True)
        ch_logp_pa = ch_mol.logp_per_atom()
        ch_mr_pa = ch_mol.mr_per_atom()
        rd_crippen_pa = rdMolDescriptors._CalcCrippenContribs(rd_mol)

        n_atoms = rd_mol.GetNumAtoms()
        for i2 in range(min(n_atoms, len(ch_hyb))):
            if rd_hyb[i2] is not None:
                accs["hybridization_per_atom"].record(smi, ch_hyb[i2], rd_hyb[i2], tol=0)
            accs["formal_charge_per_atom"].record(smi, ch_fc_pa[i2], rd_fc_pa[i2], tol=0)
            accs["implicit_hcount_per_atom"].record(smi, ch_ih_pa[i2], rd_ih_pa[i2], tol=0)
            accs["tpsa_per_atom"].record(smi, ch_tpsa_pa[i2], rd_tpsa_pa[i2], tol=0.1)
            accs["logp_crippen_per_atom"].record(smi, ch_logp_pa[i2], rd_crippen_pa[i2][0], tol=0.01)
            accs["mr_per_atom"].record(smi, ch_mr_pa[i2], rd_crippen_pa[i2][1], tol=0.01)

        # autocorr_2d / usrcat — no oracle, just validity
        ch_ac2d = ch_mol.autocorr_2d()
        for i2 in range(min(7, len(ch_ac2d))):
            accs[f"autocorr_2d_lag{i2+1}"].record_no_oracle(smi, ch_ac2d[i2], None)
        ch_usrcat = ch_mol.usrcat()
        for i2 in range(min(42, len(ch_usrcat))):
            accs[f"usrcat_{i2+1}"].record_no_oracle(smi, ch_usrcat[i2], None)

    # unbound rows (aligned by line index with the corpus file, NOT filtered by parse success above -
    # the Rust harness has its own independent parse step)
    for row in unbound_rows:
        smi = row.get("smiles", "")
        if not row.get("parse_ok"):
            continue
        moran = row.get("moran_autocorr")
        if moran is not None:
            for i2 in range(min(7, len(moran))):
                accs[f"moran_autocorr_lag{i2+1}"].record_no_oracle(smi, moran[i2], None)
        else:
            for i2 in range(7):
                accs[f"moran_autocorr_lag{i2+1}"].fixture_count += 1
                accs[f"moran_autocorr_lag{i2+1}"].ch_error += 1
        geary = row.get("geary_autocorr")
        if geary is not None:
            for i2 in range(min(7, len(geary))):
                accs[f"geary_autocorr_lag{i2+1}"].record_no_oracle(smi, geary[i2], None)
        else:
            for i2 in range(7):
                accs[f"geary_autocorr_lag{i2+1}"].fixture_count += 1
                accs[f"geary_autocorr_lag{i2+1}"].ch_error += 1
        for n in ["ic", "tic", "sic", "bic", "cic"]:
            v = row.get(n)
            if v is not None:
                accs[f"information_content_{n}"].record_no_oracle(smi, v, None)
            else:
                accs[f"information_content_{n}"].fixture_count += 1
                accs[f"information_content_{n}"].ch_error += 1
        mde = row.get("mde_carbon")
        if mde is not None:
            for i2 in range(min(10, len(mde))):
                accs[f"mde_carbon_{i2+1}"].record_no_oracle(smi, mde[i2], None)
        else:
            for i2 in range(10):
                accs[f"mde_carbon_{i2+1}"].fixture_count += 1
                accs[f"mde_carbon_{i2+1}"].ch_error += 1
        mmff = row.get("mmff94_charges_sum")
        if mmff is not None:
            accs["mmff94_charges_plain"].record_no_oracle(smi, mmff, None)
        else:
            accs["mmff94_charges_plain"].fixture_count += 1
            accs["mmff94_charges_plain"].ch_error += 1

    # ---- finalize ----
    out = {}
    for name, acc in accs.items():
        tol_desc = {
            "molecular_weight": "±0.01", "exact_mass": "±0.01", "tpsa": "±0.1",
            "logp_crippen": "±0.01", "molar_refractivity": "±0.01", "fsp3": "±0.001",
            "fraction_rotatable_bonds": "±1e-6", "balaban_j": "±0.01", "ipc": "±0.01",
            "hall_kier_alpha": "±0.01", "tpsa_per_atom": "±0.1",
            "logp_crippen_per_atom": "±0.01", "mr_per_atom": "±0.01",
        }.get(name, "exact" if acc.oracle else "n/a (no oracle)")
        out[name] = acc.finalize(tol_desc)

    summary = {
        "generated_at": __import__("datetime").datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
        "chematic_version": getattr(chematic, "__version__", "unknown"),
        "rdkit_version": __import__("rdkit").__version__,
        "scope": "crates/chematic-chem/src/descriptors.rs only (71 functions) -- see docs/rfcs/descriptor_census_rfc.md",
        "corpus": {
            "source": args.corpus,
            "total_lines": len(smiles_list),
            "rdkit_parse_failures": rd_parse_fail,
            "chematic_parse_failures": ch_parse_fail,
            "evaluated": total,
            "extended_metrics_skipped": extended_skipped,
        },
        "descriptors": out,
    }

    with open(args.json, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Evaluated {total}/{len(smiles_list)} molecules.")
    print(f"JSON written to {args.json}")

    # quick human-readable summary to stdout
    print(f"\n{'name':<30} {'oracle':<8} {'exact%':>8} {'MAE':>10} {'max_AE':>10}  classification")
    for name, dd in out.items():
        oracle_flag = "yes" if dd["rdkit_oracle"] else "no"
        exact_pct = f"{dd['exact_match_pct']:.2f}" if dd["exact_match_pct"] is not None else "-"
        mae = f"{dd['mae']:.4g}" if dd["mae"] is not None else "-"
        maxae = f"{dd['max_ae']:.4g}" if dd["max_ae"] is not None else "-"
        print(f"{name:<30} {oracle_flag:<8} {exact_pct:>8} {mae:>10} {maxae:>10}  {dd['classification']}")


if __name__ == "__main__":
    main()
