#!/usr/bin/env python3
"""Run chematic reactions for another process, one JSON request per line.

Published chematic Linux wheels are CPython 3.9 builds while RDKit
2026.03.6 ships no CPython 3.9 wheel, so a harness that compares the two on
a published artifact runs RDKit itself and chematic in this worker, under the
interpreter that has the published wheel:

    {"smirks": "...", "smiles": "..."}  ->  {"status": ..., "detail": ...,
                                             "products": [[product, ...], ...]}

Each product carries ``smiles``, and for products with order-less ``~``
bonds also ``ordered`` (``smiles_with_atom_order()``) and ``molblock``.
Only the standard library and chematic are imported (Python 3.9 syntax).
"""

import json
import sys

import chematic


def product_record(product):
    record = {"smiles": product.smiles}
    if "~" in product.smiles:
        smi, order = product.smiles_with_atom_order()
        record["ordered"] = [smi, list(order)]
        record["molblock"] = product.to_mol_block()
    return record


def run(request):
    try:
        reactant = chematic.from_smiles(request["smiles"])
    except ValueError as exc:
        return {"status": "reactant_parse", "detail": str(exc), "products": None}
    try:
        checked = chematic.run_smirks_checked(request["smirks"], [reactant], rdkit_compat=True)
    except ValueError as exc:
        return {"status": "error", "detail": str(exc), "products": None}
    status = checked["status"]
    if status in {"typed_refusal", "typed_unsupported"}:
        return {"status": status, "detail": checked.get("reason"), "products": None}
    products = [[product_record(p) for p in product_set] for product_set in checked.get("products", [])]
    return {"status": status, "detail": None, "products": products}


def main():
    sys.stdout.write(json.dumps({"version": chematic.__version__, "file": chematic.__file__}) + "\n")
    sys.stdout.flush()
    for line in sys.stdin:
        if not line.strip():
            continue
        sys.stdout.write(json.dumps(run(json.loads(line))) + "\n")
        sys.stdout.flush()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
