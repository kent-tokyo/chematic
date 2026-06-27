from .chematic import *
from .chematic import __version__, Mol, SimilarityIndex, SdfRecord, SdfIter, bulk
from .chematic import (
    from_smiles,
    from_mol_block,
    from_inchi,
    is_valid_smiles,
    tanimoto,
    smarts_match,
    smarts_find,
    depict_grid,
    run_smirks,
    find_mcs,
    iter_sdf,
    iter_sdf_str,
)


# Task → representation mapping (arXiv 2026: CML/MolJSON outperform SMILES on structural tasks)
_TASK_REPR = {
    "structural_reasoning": "moljson",
    "shortest_path":        "moljson",
    "graph_reasoning":      "moljson",
    "identification":       "inchi",
    "exact_match":          "inchi",
    "property_prediction":  "canonical_smiles",
    "generation":           "canonical_smiles",
    "editing":              "cml",
    "default":              "canonical_smiles",
}


def best_representation(task: str = "default") -> str:
    """Return the recommended molecular text format for an LLM task.

    Based on arXiv 2026 "Rethinking Molecular Text Representations for LLMs":
    CML and MolJSON outperform SMILES on structural reasoning tasks; InChI is
    best for exact identification.

    Args:
        task: one of ``structural_reasoning``, ``shortest_path``,
              ``graph_reasoning``, ``identification``, ``exact_match``,
              ``property_prediction``, ``generation``, ``editing``,
              ``default``

    Returns:
        format string — pass directly to ``mol.to_llm_text(format)``

    Example::

        fmt = chematic.best_representation("structural_reasoning")
        # → "moljson"
        text = mol.to_llm_text(fmt)
    """
    return _TASK_REPR.get(task, "canonical_smiles")


def from_smiles_list(smiles, /, *, skip_invalid=True):
    """Parse a list of SMILES strings into Mol objects.

    Runs in parallel (Rayon). Invalid SMILES are silently dropped by default.

    Args:
        smiles: iterable of SMILES strings
        skip_invalid: if True (default), drop invalid entries; if False, keep None

    Returns:
        list of Mol objects (or list[Mol | None] when skip_invalid=False)

    Example::

        mols = chematic.from_smiles_list(["CCO", "c1ccccc1", "INVALID"])
        # → [<Mol CCO>, <Mol c1ccccc1>]
    """
    parsed = bulk.parse(list(smiles))
    if skip_invalid:
        return [m for m in parsed if m is not None]
    return parsed


class MolContextPack:
    """Assembled molecular context for LLM prompts and RAG pipelines.

    Created via ``chematic.context_pack(mol)``.

    Combines identifiers, physicochemical properties, drug-likeness flags,
    structural alerts, and molecular representations into a single object
    ready to embed in an LLM prompt.

    Example::

        mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
        ctx = chematic.context_pack(mol)
        print(ctx.to_markdown())         # LLM-ready markdown block
        print(ctx.to_json())             # structured JSON for RAG retrieval
        prompt_text = ctx.to_prompt()    # compact single-line summary
    """

    def __init__(self, mol):
        self._mol = mol
        d = mol.descriptors()
        self._data = {
            "identifiers": {
                "smiles": mol.smiles,
                "inchi": mol.inchi,
                "inchikey": mol.inchikey,
                "formula": mol.formula,
                "mw": round(float(d.get("mw", 0)), 2),
                "exact_mass": round(float(d.get("exact_mass", 0)), 4),
            },
            "properties": {
                "logp": round(float(d.get("logp", 0)), 2),
                "tpsa": round(float(d.get("tpsa", 0)), 2),
                "hbd": int(d.get("hbd", 0)),
                "hba": int(d.get("hba", 0)),
                "rotatable_bonds": int(d.get("rotatable_bonds", 0)),
                "aromatic_rings": int(d.get("aromatic_ring_count", 0)),
                "heavy_atoms": int(d.get("heavy_atom_count", 0)),
                "fsp3": round(float(d.get("fsp3", 0)), 3),
                "qed": round(float(d.get("qed", 0)), 3),
                "sa_score": round(float(d.get("sa_score", 0)), 2),
            },
            "drug_likeness": {
                "lipinski_passes": bool(mol.lipinski_passes),
                "pains_passes": bool(mol.pains_passes),
                "brenk_passes": bool(mol.brenk_passes),
                "egan_passes": bool(mol.egan_passes),
                "veber_passes": bool(mol.veber_passes),
                "pains_alerts": list(mol.pains_alerts()),
                "brenk_alerts": list(mol.brenk_alerts()),
            },
            "admet": {
                "bbb_score": round(float(mol.bbb_score), 3),
                "caco2": mol.caco2,
                "cyp3a4_risk": mol.cyp3a4_risk,
                "clearance_class": mol.clearance_class,
            },
            "representations": {
                "canonical_smiles": mol.smiles,
                "inchi": mol.inchi,
                "moljson": mol.to_moljson(),
            },
        }

    def to_dict(self) -> dict:
        """Return the context as a plain Python dict."""
        return self._data

    def to_json(self, indent: int = 2) -> str:
        """Return the context as a JSON string."""
        import json
        return json.dumps(self._data, indent=indent, ensure_ascii=False)

    def to_markdown(self) -> str:
        """Return the context as a Markdown block for LLM prompts."""
        d = self._data
        p = d["properties"]
        dl = d["drug_likeness"]
        ids = d["identifiers"]
        adm = d["admet"]

        alerts = []
        if dl["pains_alerts"]:
            alerts.append(f"PAINS: {', '.join(dl['pains_alerts'])}")
        if dl["brenk_alerts"]:
            alerts.append(f"Brenk: {', '.join(dl['brenk_alerts'])}")
        alerts_str = "; ".join(alerts) if alerts else "none"

        drug_flags = []
        for name, label in [
            ("lipinski_passes", "Lipinski Ro5"),
            ("egan_passes", "Egan"),
            ("veber_passes", "Veber"),
        ]:
            drug_flags.append(f"{'✓' if dl[name] else '✗'} {label}")

        return (
            f"## Molecule\n"
            f"- **SMILES**: {ids['smiles']}\n"
            f"- **Formula**: {ids['formula']}  MW: {ids['mw']} Da\n"
            f"- **InChI**: {ids['inchi']}\n\n"
            f"## Physicochemical Properties\n"
            f"- LogP: {p['logp']}, TPSA: {p['tpsa']} Å², "
            f"HBD: {p['hbd']}, HBA: {p['hba']}\n"
            f"- Rotatable bonds: {p['rotatable_bonds']}, "
            f"Aromatic rings: {p['aromatic_rings']}, "
            f"fsp3: {p['fsp3']}\n"
            f"- QED: {p['qed']} (drug-likeness), "
            f"SA score: {p['sa_score']} (synthesizability)\n\n"
            f"## Drug-Likeness\n"
            f"- {' | '.join(drug_flags)}\n"
            f"- Structural alerts: {alerts_str}\n\n"
            f"## ADMET\n"
            f"- BBB score: {adm['bbb_score']}, "
            f"Caco-2: {adm['caco2']}, "
            f"CYP3A4: {adm['cyp3a4_risk']}, "
            f"Clearance: {adm['clearance_class']}\n\n"
            f"## Structural Representation (MolJSON)\n"
            f"```json\n{d['representations']['moljson']}\n```"
        )

    def to_prompt(self) -> str:
        """Return a compact one-line molecule summary for inline LLM prompts."""
        d = self._data
        p = d["properties"]
        ids = d["identifiers"]
        dl = d["drug_likeness"]
        flags = "Lipinski✓" if dl["lipinski_passes"] else "Lipinski✗"
        return (
            f"{ids['smiles']} | MW={ids['mw']} | "
            f"LogP={p['logp']} | TPSA={p['tpsa']} | "
            f"HBD={p['hbd']} HBA={p['hba']} | "
            f"QED={p['qed']} | {flags}"
        )

    def __repr__(self) -> str:
        return f"MolContextPack({self._data['identifiers']['smiles']!r})"


def context_pack(mol) -> "MolContextPack":
    """Assemble a molecule context pack for LLM prompts and RAG pipelines.

    Combines identifiers, physicochemical properties, drug-likeness flags,
    structural alerts, ADMET profile, and molecular representations.

    Args:
        mol: a ``Mol`` object (from ``from_smiles``, ``from_moljson``, etc.)

    Returns:
        ``MolContextPack`` with ``.to_markdown()``, ``.to_json()``,
        ``.to_dict()``, ``.to_prompt()`` methods.

    Example::

        mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
        ctx = chematic.context_pack(mol)

        # For RAG retrieval:
        vector_store.add(ctx.to_json(), embedding=mol.hdf())

        # For LLM prompt:
        prompt = f"Analyse this molecule:\\n{ctx.to_markdown()}"
    """
    return MolContextPack(mol)


def descriptors_df(smiles):
    """Compute 55+ descriptors for a list of SMILES and return a DataFrame.

    Requires pandas (``pip install pandas``). Runs in parallel via Rayon.

    Args:
        smiles: iterable of SMILES strings (invalid entries are skipped)

    Returns:
        pd.DataFrame with one row per valid molecule and 55+ descriptor columns
        (mw, logp, tpsa, hbd, hba, qed, sa_score, pains_passes, …)

    Example::

        df = chematic.descriptors_df(["CCO", "c1ccccc1", "CC(=O)O"])
        df[["mw", "logp", "tpsa"]].head()
    """
    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas is required: pip install pandas") from None
    return pd.DataFrame(bulk.descriptors(list(smiles)))


def fragment_text(mol, method: str = "brics", fmt: str = "markdown") -> str:
    """Describe a molecule as its structural fragments for LLM prompts.

    Decomposes the molecule using BRICS retrosynthetic rules and returns a
    human-readable description of the scaffold, fragments, and connection
    points — useful for medicinal chemistry reasoning with LLMs.

    Based on: MolLingo (arXiv 2026) — block-based SMILES + fragment names.

    Args:
        mol: a ``Mol`` object
        method: fragmentation method; currently ``"brics"`` (default)
        fmt: output format — ``"markdown"`` (default) or ``"json"``

    Returns:
        str — fragment description for LLM prompts

    Example::

        mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
        print(chematic.fragment_text(mol))
        # ## Fragment Decomposition
        # - **Scaffold**: c1ccccc1
        # - **Fragments** (BRICS):
        #   1. C(=O)(O[*])C  — acetyl ester
        #   2. C(=O)(O)c1ccccc1[*]  — benzoic acid core
        # - **Connection points**: 1 ester linkage
    """
    frags = mol.brics_fragments()
    scaffold_smi = mol.scaffold().smiles
    n_bonds = len(mol.brics_bonds())

    if fmt == "json":
        import json
        return json.dumps({
            "smiles": mol.smiles,
            "scaffold": scaffold_smi,
            "fragments": [f.smiles for f in frags],
            "connection_points": n_bonds,
            "method": method,
        }, indent=2)

    # markdown
    frag_lines = "\n".join(
        f"  {i + 1}. `{f.smiles}`"
        for i, f in enumerate(frags)
    ) if frags else "  (no fragmentation sites)"

    conn = f"{n_bonds} bond{'s' if n_bonds != 1 else ''}" if n_bonds else "none"

    return (
        f"## Fragment Decomposition ({method.upper()})\n"
        f"- **Molecule**: `{mol.smiles}`\n"
        f"- **Scaffold**: `{scaffold_smi}`\n"
        f"- **Fragments**:\n{frag_lines}\n"
        f"- **Connection points**: {conn}\n"
        f"- **Fragment count**: {len(frags)}"
    )
