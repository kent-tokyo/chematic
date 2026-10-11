"""Repeat the source-extension quality diagnostic; no published-wheel claim."""

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import chematic
from scripts.public_package_3d_chematic import load_rows, pipeline_config, run_pipeline


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = pipeline_config(chematic, "mmff94_bond_angle_strict", stereo_safe=True)
    rows = load_rows("AB")
    with args.output.open("w") as handle:
        for index, row in enumerate(rows):
            result = run_pipeline(chematic.from_smiles(row["smiles"]), config)
            record = {
                **row,
                "row_index": index,
                "engine": "chematic_source_candidate_python",
                "arm": "chematic_pipeline_v2_mmff94_strict_stereo_safe",
                **result,
            }
            handle.write(json.dumps(record, separators=(",", ":")) + "\n")
            handle.flush()
    print(f"completed {len(rows)} source-candidate quality rows")


if __name__ == "__main__":
    main()
