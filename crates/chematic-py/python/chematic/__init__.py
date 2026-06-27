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
