#!/usr/bin/env python3
"""Generate the frozen Rust/Python/WASM cross-binding parity fixture for
Mol.embed_pipeline_v2() / embed_pipeline_v2_json().

Computes each fixture via the Python binding (crates/chematic-py/src/pipeline_v2.rs)
and dumps only the fields the three bindings are required to agree on: success/
failure, atom count, stage, cause kind, stereo counts, force-field requested/
actual/fallback presence, final soundness. Wall-clock timing is deliberately
excluded (never required to match across bindings or even across repeated runs).

Output is read by:
- crates/chematic-wasm/src/pipeline_v2.rs's own #[cfg(test)] suite (WASM vs. this
  reference, and Rust's raw chematic_3d::pipeline_v2::embed_pipeline_v2 vs. this
  reference -- both bindings call the same underlying Rust function, so this
  mainly catches a JSON-conversion mistake in either binding, not an algorithm
  divergence).
- crates/chematic-wasm/tests/pipeline_v2_parity.test.mjs (WASM, via node).

Regenerate with:
    .venv/bin/maturin develop --release -m crates/chematic-py/Cargo.toml
    .venv/bin/python scripts/gen_pipeline_v2_wasm_parity_fixtures.py
"""

import json
import pathlib

import chematic

OUT_PATH = pathlib.Path(__file__).resolve().parent.parent / "validation" / "pipeline_v2_wasm_parity_fixtures.json"

# (name, smiles, config kwargs for PipelineV2Config.safe)
FIXTURES = [
    ("decane", "CCCCCCCCCC", dict(force_field="none", stereo_policy="ignore", ring_torsion_policy="fail_closed")),
    ("naphthalene", "c1ccc2ccccc2c1", dict(force_field="none", stereo_policy="ignore", ring_torsion_policy="fail_closed")),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O", dict(force_field="none", stereo_policy="ignore", ring_torsion_policy="fail_closed")),
    ("branched", "CCC(C)C", dict(force_field="none", stereo_policy="ignore", ring_torsion_policy="fail_closed")),
    ("declared_tetrahedral_stereo", "F[C@@](Cl)(Br)I", dict(force_field="none", stereo_policy="verify_only", ring_torsion_policy="fail_closed")),
    ("ez_double_bond", "C/C=C/C", dict(force_field="none", stereo_policy="verify_only", ring_torsion_policy="fail_closed")),
    ("ring_torsion_fail_closed", "C1CCCCC1CCCCCCCCCCCC", dict(force_field="dreiding", stereo_policy="ignore", ring_torsion_policy="fail_closed", use_small_ring_torsions=True)),
    ("force_field_fallback", "c1ccc2ccccc2c1", dict(force_field="mmff94_with_uff_fallback", stereo_policy="ignore", ring_torsion_policy="fail_closed")),
]

EMBED_SEED = 7


def config_json_camel_case(kwargs: dict) -> dict:
    """The same config, expressed as the camelCase JSON object WASM's
    embed_pipeline_v2_json() and Rust's own PipelineV2ConfigJson expect --
    filled in with .safe()'s own conservative defaults for anything not
    explicit in `kwargs`, so all three bindings run the literal same config."""
    return {
        "embedSeed": EMBED_SEED,
        "maxAttempts": 8,
        "embedTimeoutMs": None,
        "useExpTorsions": kwargs.get("use_exp_torsions", False),
        "useSmallRingTorsions": kwargs.get("use_small_ring_torsions", False),
        "useMacrocycleTorsions": kwargs.get("use_macrocycle_torsions", False),
        "useMacrocycle14Bounds": kwargs.get("use_macrocycle_14_bounds", False),
        "includeLegacyTorsionHeuristic": False,
        "stereoPolicy": kwargs["stereo_policy"],
        "failOnUnevaluableStereo": False,
        "forceFieldPolicy": kwargs["force_field"],
        "forceFieldMaxIterations": 200,
        "gateMmff94TorsionOop": False,
        "gateMmff94StretchBend": kwargs.get("gate_mmff94_stretch_bend", False),
        "ringTorsionPolicy": kwargs["ring_torsion_policy"],
        "totalTimeoutMs": None,
    }


def run_fixture(name: str, smiles: str, kwargs: dict) -> dict:
    mol = chematic.from_smiles(smiles)
    config = chematic.PipelineV2Config.safe(embed_seed=EMBED_SEED, **kwargs)
    entry = {
        "name": name,
        "smiles": smiles,
        "config": config_json_camel_case(kwargs),
        "atomCount": mol.heavy_atoms,
    }
    try:
        r = mol.embed_pipeline_v2(config)
        entry["ok"] = True
        entry["coordsLength"] = len(r["coords"])
        entry["stereoBeforeDeclared"] = r["stereo_before"]["n_declared"]
        entry["stereoAfterDeclared"] = r["final_stereo"]["n_declared"]
        entry["stereoAfterViolations"] = r["final_stereo"]["n_violations"]
        entry["forceFieldRequested"] = r["force_field"]["requested_force_field"]
        entry["forceFieldActual"] = r["force_field"]["actual_force_field_used"]
        entry["hasFallback"] = r["force_field"]["fallback_reason"] is not None
        entry["sound"] = r["final_validation"]["sound"]
    except chematic.PipelineV2Error as e:
        d = e.diagnostics
        entry["ok"] = False
        entry["stage"] = d["stage"]
        entry["causeKind"] = d["cause"]["kind"]
    return entry


def main() -> None:
    fixtures = [run_fixture(name, smiles, kwargs) for name, smiles, kwargs in FIXTURES]
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps({"fixtures": fixtures}, indent=2) + "\n")
    print(f"wrote {len(fixtures)} fixtures to {OUT_PATH}")
    for f in fixtures:
        print(" ", f["name"], "ok" if f["ok"] else f"FAIL({f['causeKind']})")


if __name__ == "__main__":
    main()
