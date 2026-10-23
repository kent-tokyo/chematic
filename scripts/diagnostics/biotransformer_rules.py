#!/usr/bin/env python3
"""BioTransformer rule tables and the corpus's reactant requests, without RDKit.

Shared by ``biotransformer_rule_corpus.py`` (RDKit side) and
``biotransformer_chematic_responses.py`` (chematic side), so the two can run
on different hosts: a published wheel on a platform without an RDKit wheel
records chematic's products, and RDKit compares them on Linux.
"""

from __future__ import annotations

import hashlib
from pathlib import Path


def load_rules(paths: list[Path]) -> tuple[list[dict], list[dict]]:
    import json5

    rules, sources = [], []
    for path in paths:
        raw = path.read_bytes()
        sources.append({"path": path.name, "sha256": hashlib.sha256(raw).hexdigest()})
        text = raw.decode("utf-8")
        table = json5.loads(text[text.index("{"):], object_pairs_hook=lambda pairs: pairs)
        top = dict(table)
        entries = top.get("reactions", table)
        for name, body in entries:
            body = dict(body)
            if body.get("smirks"):
                rules.append({
                    "table": path.name,
                    "name": name,
                    "id": body.get("btmrID") or name,
                    "smirks": body["smirks"],
                })
    return rules, sources


def load_requests(path: Path) -> list[tuple[str, str]]:
    """``(implicit, explicit_h)`` reactant SMILES as RDKit 2026.03.6 writes
    them (``validation/biotransformer-requests-400.tsv``), one row per
    reactant that RDKit reads."""
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip() and not line.startswith("#"):
            implicit, explicit_h = line.split("\t")
            rows.append((implicit, explicit_h))
    return rows
