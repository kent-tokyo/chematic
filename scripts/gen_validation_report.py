#!/usr/bin/env python3
"""Generate docs/validation.md from a bench5k JSON results file.

Usage:
    python3 scripts/gen_validation_report.py validation/results/latest.json
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent


def badge(pct: float) -> str:
    color = "brightgreen" if pct >= 99.0 else "green" if pct >= 95.0 else "yellow"
    return f"![{pct:.1f}%](https://img.shields.io/badge/agreement-{pct:.1f}%25-{color})"


def main() -> None:
    if len(sys.argv) < 2:
        sys.exit(f"Usage: {sys.argv[0]} <results.json>")

    data = json.loads(Path(sys.argv[1]).read_text())
    ts = data.get("generated_at", "unknown")
    ver = data.get("chematic_version", "?")
    corpus = data["corpus"]
    m = data["metrics"]

    rows = [
        ("HBA (hydrogen bond acceptors)", m["hba"]["agreement_pct"], m["hba"]["tolerance"],
         corpus["total"], m["hba"]["match"]),
        ("HBD (hydrogen bond donors)", m["hbd"]["agreement_pct"], m["hbd"]["tolerance"],
         corpus["total"], m["hbd"]["match"]),
        ("ARC (aromatic ring count)", m["arc"]["agreement_pct"], m["arc"]["tolerance"],
         corpus["total"], m["arc"]["match"]),
        ("[nH] SMARTS", m["nh_smarts"]["agreement_pct"], "exact",
         corpus["total"], round(m["nh_smarts"]["agreement_pct"] * corpus["total"] / 100)),
        ("TPSA", m["tpsa"]["agreement_pct"], m["tpsa"]["tolerance"],
         corpus["total"], m["tpsa"]["match"]),
        ("LogP (Crippen)", m["logp"]["agreement_pct"], m["logp"]["tolerance"],
         corpus["total"], m["logp"]["match"]),
    ]

    table_rows = "\n".join(
        f"| {name} | **{pct:.1f}%** ({match}/{n}) | {tol} |"
        for name, pct, tol, n, match in rows
    )

    nh = m["nh_smarts"]

    md = f"""\
# Validation Report

RDKit parity dashboard for chematic v{ver}.

**Environment**: Python 3.12, Apple M-series, RDKit 2026.03.3 (`includeSandP=True`)
**Corpus**: {corpus['total']:,} molecules (ChEMBL-derived subset)
**Generated**: {ts}
**Reproduce**: `python3 scripts/bench5k.py ~/Downloads/SMILES.csv --json validation/results/latest.json`

---

## Descriptor Agreement vs RDKit

| Descriptor | Agreement | Tolerance |
|---|---|---|
{table_rows}

### [nH] SMARTS precision/recall

| Metric | Value |
|---|---|
| Precision (no false positives) | {nh['precision_pct']:.1f}% |
| Recall (no false negatives) | {nh['recall_pct']:.1f}% |
| TP / TN / FP / FN | {nh['tp']} / {nh['tn']} / {nh['fp']} / {nh['fn']} |

---

## Known Limitations

| Area | Status |
|---|---|
| TPSA on complex P/S heterocycles | ~99.4% (phosphazene/sulfonimide edge cases) |
| LogP for squaramide and aminocoumarin | ~99.7% (unusual ring systems) |
| 3D conformer quality vs RDKit ETKDGv3 | Rule-based; good for screening, not publication |
| Transition metals | Out of scope |

Full failure analysis: run `python3 scripts/bench5k.py ... --detail` and inspect stderr.
"""

    out = ROOT / "docs" / "validation.md"
    out.write_text(md)
    print(f"Written {out}")


if __name__ == "__main__":
    main()
