"""
AiZynthFinder + chematic Integration Tutorial
=============================================

This script demonstrates how to use chematic as a backend for
AiZynthFinder-based retrosynthetic planning workflows.

Integration points
------------------
1. Pre-processing  — chematic standardises / validates targets before
                     feeding them to AiZynthFinder
2. Fast one-step   — chematic's BRICS disconnection runs instantly
                     alongside AiZynthFinder's ML search
3. Post-scoring    — chematic scores AiZynthFinder's building-block
                     suggestions (SA score, Tanimoto to known BBs,
                     drug-likeness, ADMET)
4. Route ranking   — combine chematic scores with AiZynthFinder scores

AiZynthFinder is optional: the script runs fully with chematic alone
and activates the AiZynthFinder sections when it is installed and
configured.

Install
-------
    pip install chematic                   # always required
    pip install aizynthfinder              # optional (ML multi-step)

    # Download policy/filter model data (if using AiZynthFinder):
    download_public_data aizynthfinder_data/

Usage
-----
    python examples/aizynthfinder_integration.py
    python examples/aizynthfinder_integration.py --smiles "O=C(O)c1ccc(N)cc1"
    python examples/aizynthfinder_integration.py --config aizynthfinder_data/config.yml
"""

import argparse
import sys
from typing import Optional

# ── chematic (required) ───────────────────────────────────────────────────────
try:
    import chematic
except ImportError:
    sys.exit("chematic is not installed.  Run: pip install chematic")

# ── AiZynthFinder (optional) ──────────────────────────────────────────────────
AIZYNTHFINDER_AVAILABLE = False
try:
    from aizynthfinder.aizynthfinder import AiZynthFinder  # type: ignore
    AIZYNTHFINDER_AVAILABLE = True
except ImportError:
    pass  # graceful degradation — chematic-only mode

# ── pandas (optional, for display) ───────────────────────────────────────────
try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False


# =============================================================================
# Section 1 — Target preparation with chematic
# =============================================================================

def prepare_target(smiles: str):
    """Parse, validate, and profile a target molecule using chematic.

    Returns a chematic Mol on success, None if the SMILES is invalid.
    """
    print(f"\n{'='*60}")
    print("SECTION 1 — Target preparation with chematic")
    print(f"{'='*60}")

    try:
        mol = chematic.from_smiles(smiles)
    except Exception as exc:
        print(f"  [ERROR] Could not parse SMILES: {exc}")
        return None

    d = mol.descriptors()
    print(f"  Input SMILES      : {smiles}")
    print(f"  Canonical SMILES  : {mol.smiles}")
    print(f"  Molecular formula : {mol.formula}")
    print(f"  Heavy atoms       : {d['heavy_atoms']}")

    # Drug-likeness profile
    print(f"\n  Descriptor profile:")
    print(f"    MW          = {d['mw']:.1f}  (Lipinski ≤500)")
    print(f"    LogP        = {d['logp']:.2f}  (Lipinski ≤5)")
    print(f"    HBD / HBA   = {d['hbd']} / {d['hba']}")
    print(f"    TPSA        = {d['tpsa']:.1f} Ų")
    print(f"    QED         = {d['qed']:.3f}  (0–1, higher=better)")
    print(f"    SA score    = {d['sa_score']:.2f}  (1=easy, 10=hard to synthesise)")
    print(f"    Rot. bonds  = {d['rotatable_bonds']}")
    print(f"    Arom. rings = {d['aromatic_ring_count']}")
    print(f"    Lipinski    : {'✅ passes' if d['lipinski_passes'] else '❌ fails'}")

    # ADMET quick profile
    admet = mol.admet()
    print(f"\n  ADMET profile:")
    print(f"    BBB penetrant  : {admet.get('bbb_penetrant', '?')}")
    print(f"    GI absorbed    : {d.get('gi_absorbed', '?')}")
    print(f"    hERG risk      : {admet.get('herg_risk', '?')}")
    print(f"    CYP3A4 risk    : {admet.get('cyp3a4_risk', '?')}")

    # SA-score interpretation
    sa = d['sa_score']
    if sa < 3.0:
        print(f"\n  ✅ SA score {sa:.2f} — straightforward synthesis expected")
    elif sa < 6.0:
        print(f"\n  ⚠️  SA score {sa:.2f} — moderate synthetic complexity")
    else:
        print(f"\n  ❌ SA score {sa:.2f} — complex target; retrosynthesis strongly recommended")

    return mol


# =============================================================================
# Section 2 — Fast one-step retrosynthesis with chematic (BRICS)
# =============================================================================

def chematic_retrosynthesis(mol) -> list:
    """Run chematic's BRICS-based one-step disconnection.

    Returns a list of building block info dicts sorted by SA score.
    """
    print(f"\n{'='*60}")
    print("SECTION 2 — Fast one-step retrosynthesis (chematic BRICS)")
    print(f"{'='*60}")

    # brics_fragments() returns a list of Mol objects with dummy-atom attachment points
    try:
        frag_mols = mol.brics_fragments()
    except Exception:
        frag_mols = []

    if not frag_mols:
        print("  No BRICS disconnections found for this molecule.")
        return []

    print(f"  Found {len(frag_mols)} BRICS building blocks:")

    candidates = []
    for frag in frag_mols:
        try:
            d = frag.descriptors()
            candidates.append({
                "mol": frag,
                "smiles": frag.smiles,
                "mw": d["mw"],
                "sa_score": d["sa_score"],
                "logp": d["logp"],
                "lipinski": d["lipinski_passes"],
                "source": "BRICS",
            })
        except Exception:
            pass

    # Sort by SA score (easiest first)
    candidates.sort(key=lambda x: x["sa_score"])

    for i, c in enumerate(candidates[:8]):
        print(f"    [{i+1}] SA={c['sa_score']:.2f}  MW={c['mw']:.1f}  {c['smiles']}")

    if len(candidates) > 8:
        print(f"    ... and {len(candidates)-8} more")

    return candidates


# =============================================================================
# Section 3 — AiZynthFinder multi-step retrosynthesis (optional)
# =============================================================================

def _mock_aizynthfinder_routes(smiles: str) -> list:
    """Return realistic-looking mock routes when AiZynthFinder is absent.

    This demonstrates the data structure you'd work with from the real tool.
    """
    return [
        {
            "route_id": 1,
            "score": 0.87,
            "num_steps": 2,
            "building_blocks": ["c1ccc(N)cc1", "CC(=O)Cl"],
            "reactions": ["amide coupling"],
            "source": "AiZynthFinder (mock)",
        },
        {
            "route_id": 2,
            "score": 0.74,
            "num_steps": 3,
            "building_blocks": ["c1ccc(N)cc1", "OCC(=O)O", "ClCCl"],
            "reactions": ["esterification", "amide coupling"],
            "source": "AiZynthFinder (mock)",
        },
    ]


def aizynthfinder_search(smiles: str, config_file: Optional[str] = None) -> list:
    """Run AiZynthFinder multi-step retrosynthesis.

    Falls back to a mock when AiZynthFinder is not installed.
    """
    print(f"\n{'='*60}")
    print("SECTION 3 — Multi-step retrosynthesis (AiZynthFinder)")
    print(f"{'='*60}")

    if not AIZYNTHFINDER_AVAILABLE:
        print("  AiZynthFinder not installed — showing mock output.")
        print("  Install: pip install aizynthfinder")
        print("  Then download models: download_public_data aizynthfinder_data/")
        routes = _mock_aizynthfinder_routes(smiles)
    elif config_file is None:
        print("  No config file provided — showing mock output.")
        print("  Pass --config path/to/config.yml to use real AiZynthFinder.")
        routes = _mock_aizynthfinder_routes(smiles)
    else:
        print(f"  Running AiZynthFinder with config: {config_file}")
        try:
            finder = AiZynthFinder(configfile=config_file)
            finder.target_smiles = smiles
            finder.tree_search()
            finder.build_routes()

            routes = []
            for i, route in enumerate(finder.routes):
                # Extract leaf molecules (building blocks)
                bb_smiles = []
                try:
                    for node in route.molecule_nodes:
                        if node.is_leaf:
                            bb_smiles.append(node.mol.smiles)
                except AttributeError:
                    pass

                routes.append({
                    "route_id": i + 1,
                    "score": getattr(route, "score", 0.0),
                    "num_steps": getattr(route, "num_steps", 0),
                    "building_blocks": bb_smiles,
                    "reactions": [],
                    "source": "AiZynthFinder",
                })
        except Exception as exc:
            print(f"  [ERROR] AiZynthFinder failed: {exc}")
            routes = _mock_aizynthfinder_routes(smiles)

    print(f"\n  Found {len(routes)} retrosynthetic route(s):")
    for r in routes:
        print(f"\n    Route {r['route_id']}  (score={r['score']:.2f}, {r['num_steps']} steps)")
        print(f"      Building blocks: {r['building_blocks']}")
        if r['reactions']:
            print(f"      Reactions: {' → '.join(r['reactions'])}")

    return routes


# =============================================================================
# Section 4 — Score building blocks with chematic
# =============================================================================

# Representative building-block library (eMolecules / Sigma-Aldrich subset)
KNOWN_BB_LIBRARY = [
    "c1ccc(N)cc1",          # aniline
    "CC(=O)Cl",             # acetyl chloride
    "OCC(=O)O",             # glycolic acid
    "c1ccc(O)cc1",          # phenol
    "NCc1ccccc1",           # benzylamine
    "CC(=O)OC(C)=O",        # acetic anhydride
    "c1ccncc1",             # pyridine
    "NCCc1ccccc1",          # 2-phenylethylamine
    "OC(=O)c1ccccc1",       # benzoic acid
    "CC(N)C(=O)O",          # alanine
]

_BB_FPS_CACHE = None


def _get_bb_fps():
    global _BB_FPS_CACHE
    if _BB_FPS_CACHE is None:
        fps = []
        for s in KNOWN_BB_LIBRARY:
            try:
                fps.append(chematic.from_smiles(s).ecfp4())
            except Exception:
                fps.append(None)
        _BB_FPS_CACHE = fps
    return _BB_FPS_CACHE


def score_building_blocks(routes: list) -> list:
    """Score each building block in every route using chematic."""
    print(f"\n{'='*60}")
    print("SECTION 4 — Building block scoring (chematic)")
    print(f"{'='*60}")

    bb_fps_library = _get_bb_fps()
    scored_routes = []

    for route in routes:
        scored_bbs = []
        for smi in route["building_blocks"]:
            # Strip BRICS dummy atoms (*) if present for clean parsing
            clean_smi = smi.replace("[*]", "[H]")
            try:
                mol = chematic.from_smiles(clean_smi)
                d = mol.descriptors()
                fp = mol.ecfp4()

                # Tanimoto similarity to known BB library
                sims = []
                for lib_fp in bb_fps_library:
                    if lib_fp is not None:
                        try:
                            sims.append(chematic.tanimoto(fp, lib_fp))
                        except Exception:
                            sims.append(0.0)
                    else:
                        sims.append(0.0)

                max_sim = max(sims) if sims else 0.0
                best_idx = sims.index(max_sim) if sims else 0
                best_match = KNOWN_BB_LIBRARY[best_idx] if sims else "?"

                bb_info = {
                    "smiles": smi,
                    "sa_score": d["sa_score"],
                    "mw": d["mw"],
                    "logp": d["logp"],
                    "lipinski_ok": d["lipinski_passes"],
                    "max_sim_to_bb": max_sim,
                    "best_match": best_match,
                    "available": max_sim >= 0.85,
                }
                scored_bbs.append(bb_info)
            except Exception as exc:
                scored_bbs.append({"smiles": smi, "error": str(exc)})

        route_copy = dict(route)
        route_copy["scored_bbs"] = scored_bbs
        scored_routes.append(route_copy)

        print(f"\n  Route {route['route_id']} building blocks:")
        header = f"    {'SMILES':<30} {'SA':>5} {'MW':>7} {'Lipinski':>9} {'MaxSim':>7} {'Available':>10}"
        print(header)
        print("    " + "-" * 62)
        for bb in scored_bbs:
            if "error" in bb:
                print(f"    {bb['smiles']:<30}  (parse error: {bb['error'][:30]})")
                continue
            avail = "✅ yes" if bb["available"] else "❌ no"
            lip = "✅" if bb["lipinski_ok"] else "❌"
            print(
                f"    {bb['smiles']:<30} {bb['sa_score']:>5.2f} {bb['mw']:>7.1f}"
                f" {lip:>9} {bb['max_sim_to_bb']:>7.3f} {avail:>10}"
            )
            if bb["available"]:
                print(f"       → closest library match: {bb['best_match']}")

    return scored_routes


# =============================================================================
# Section 5 — Route ranking and recommendation
# =============================================================================

def rank_routes(scored_routes: list) -> None:
    """Rank retrosynthetic routes by overall feasibility score."""
    print(f"\n{'='*60}")
    print("SECTION 5 — Route ranking")
    print(f"{'='*60}")

    rankings = []
    for route in scored_routes:
        bbs = [bb for bb in route.get("scored_bbs", []) if "error" not in bb]
        if not bbs:
            continue

        avg_sa = sum(b["sa_score"] for b in bbs) / len(bbs)
        frac_available = sum(1 for b in bbs if b.get("available", False)) / len(bbs)

        # Composite score: lower is better
        composite = (
            avg_sa * 0.4
            + (1.0 - frac_available) * 3.0
            + route.get("num_steps", 1) * 0.2
            - route.get("score", 0.0) * 0.5
        )

        rankings.append({
            "route_id": route["route_id"],
            "avg_sa": avg_sa,
            "frac_available": frac_available,
            "num_steps": route.get("num_steps", "?"),
            "aizf_score": route.get("score", 0.0),
            "composite": composite,
        })

    rankings.sort(key=lambda x: x["composite"])

    print(f"\n  Route ranking (lower composite = more feasible):\n")
    print(f"  {'Rank':>4}  {'Route':>5}  {'Steps':>5}  {'Avg SA':>7}  "
          f"{'BBs avail':>9}  {'AiZf score':>10}  {'Composite':>9}")
    print("  " + "-" * 60)
    for rank, r in enumerate(rankings, 1):
        avail_pct = f"{r['frac_available']*100:.0f}%"
        print(
            f"  {rank:>4}  {r['route_id']:>5}  {str(r['num_steps']):>5}  "
            f"{r['avg_sa']:>7.2f}  {avail_pct:>9}  {r['aizf_score']:>10.2f}  "
            f"{r['composite']:>9.3f}"
        )

    if rankings:
        best = rankings[0]
        print(f"\n  ✅ Recommended route: Route {best['route_id']}")
        print(f"     Avg SA score  : {best['avg_sa']:.2f}")
        print(f"     BBs available : {best['frac_available']*100:.0f}%")
        print(f"     Steps         : {best['num_steps']}")


# =============================================================================
# Section 6 — chematic-only workflow summary
# =============================================================================

def chematic_only_summary(mol, brics_candidates: list) -> None:
    """Show a self-contained chematic workflow summary."""
    print(f"\n{'='*60}")
    print("SECTION 6 — chematic standalone summary")
    print(f"{'='*60}")

    if not brics_candidates:
        print("  No BRICS candidates to display.")
        return

    top = brics_candidates[:5]
    print(f"\n  Top {len(top)} building blocks by SA score (via BRICS):\n")

    if PANDAS_AVAILABLE:
        rows = []
        for bb in top:
            try:
                d = bb["mol"].descriptors()
                rows.append({
                    "SMILES": bb["smiles"],
                    "SA score": round(d["sa_score"], 2),
                    "MW": round(d["mw"], 1),
                    "LogP": round(d["logp"], 2),
                    "HBD/HBA": f"{d['hbd']}/{d['hba']}",
                    "Lipinski": "✅" if d["lipinski_passes"] else "❌",
                })
            except Exception:
                pass
        if rows:
            df = pd.DataFrame(rows)
            print(df.to_string(index=False))
    else:
        for bb in top:
            print(f"    SA={bb['sa_score']:.2f}  MW={bb['mw']:.1f}  {bb['smiles']}")

    if not AIZYNTHFINDER_AVAILABLE:
        print(f"\n  To add AiZynthFinder multi-step retrosynthesis:")
        print(f"    pip install aizynthfinder")
        print(f"    download_public_data aizynthfinder_data/")
        print(f"    python examples/aizynthfinder_integration.py \\")
        print(f"        --smiles '{mol.smiles}' \\")
        print(f"        --config aizynthfinder_data/config.yml")


# =============================================================================
# Main
# =============================================================================

# Sulfanilamide derivative: moderate complexity, drug-like
DEFAULT_TARGET = "O=C(Nc1ccc(S(N)(=O)=O)cc1)c1ccc(N)cc1"


def main():
    parser = argparse.ArgumentParser(
        description="AiZynthFinder + chematic retrosynthesis tutorial"
    )
    parser.add_argument(
        "--smiles",
        default=DEFAULT_TARGET,
        help="Target molecule SMILES (default: sulfanilamide derivative)",
    )
    parser.add_argument(
        "--config",
        default=None,
        help="AiZynthFinder config.yml path (omit for mock/chematic-only mode)",
    )
    args = parser.parse_args()

    print("\n" + "=" * 60)
    print(" AiZynthFinder + chematic Integration Tutorial")
    print("=" * 60)
    print(f"  chematic version : {getattr(chematic, '__version__', 'installed')}")
    print(f"  AiZynthFinder    : {'installed ✅' if AIZYNTHFINDER_AVAILABLE else 'not installed (mock mode)'}")
    print(f"  pandas           : {'installed ✅' if PANDAS_AVAILABLE else 'not installed'}")

    # 1. Prepare target
    mol = prepare_target(args.smiles)
    if mol is None:
        sys.exit(1)

    # 2. chematic BRICS one-step (always runs)
    brics_candidates = chematic_retrosynthesis(mol)

    # 3. AiZynthFinder multi-step (real or mock)
    routes = aizynthfinder_search(mol.smiles, args.config)

    # 4. Score building blocks with chematic
    scored_routes = score_building_blocks(routes)

    # 5. Rank routes
    rank_routes(scored_routes)

    # 6. Standalone summary
    chematic_only_summary(mol, brics_candidates)

    print(f"\n{'='*60}")
    print(" Done.")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
