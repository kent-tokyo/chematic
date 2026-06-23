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
