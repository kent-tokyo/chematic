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
from .chematic import (
    Lattice,
    PeriodicStructure,
    Site,
    PeriodicNeighbor,
    CifSymmetryStatus,
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
    parsed = bulk.parse(_materialize_smiles_batch(smiles))
    if skip_invalid:
        return [m for m in parsed if m is not None]
    return parsed


_FORMAT_ALIASES = {
    "smi": "smiles", "smiles": "smiles", "mol": "mol", "sdf": "mol",
    "mol_v3000": "mol_v3000", "v3000": "mol_v3000", "mol2": "mol2",
    "cml": "cml", "cjson": "cjson", "moljson": "moljson", "cdxml": "cdxml",
    "pdb": "pdb", "xyz": "xyz", "pdbqt": "pdbqt", "gjf": "gjf", "com": "gjf",
}
_CONVERT_MAX_INPUT_BYTES = 16 * 1024 * 1024
_MAX_BATCH_ITEMS = 100_000


def _materialize_smiles_batch(smiles):
    values = list(smiles)
    if len(values) > _MAX_BATCH_ITEMS:
        raise ValueError(
            f"SMILES batch exceeds maximum item count ({len(values)} > "
            f"{_MAX_BATCH_ITEMS})"
        )
    return values


def convert_format(text, input_format, output_format, /, *, coords=None,
                   charges=None, name="LIG", comment=""):
    """Convert a molecule between common interchange formats.

    Supported inputs are SMILES, MOL/SDF (V2000), MOL V3000, MOL2, CML,
    ChemicalJSON, MolJSON, CDXML, PDB, XYZ, PDBQT, and Gaussian input.
    Supported outputs are SMILES, MOL V2000, MOL V3000, MOL2, CML,
    ChemicalJSON, MolJSON, CDXML, PDB, XYZ, and PDBQT. PDB, XYZ, and PDBQT
    output require ``coords=[[x, y, z], ...]`` unless the input provides them.

    The molecular graph and supported stereochemistry are preserved. Format-
    specific metadata is not promised to survive conversion.
    """
    def normalize(value):
        if not isinstance(value, str):
            raise ValueError("format names must be strings")
        key = value.lower().lstrip(".").replace("-", "_")
        try:
            return _FORMAT_ALIASES[key]
        except KeyError as exc:
            supported = ", ".join(sorted(set(_FORMAT_ALIASES.values())))
            raise ValueError(
                f"unsupported molecular format {value!r}; supported formats: {supported}"
            ) from exc

    source, target = normalize(input_format), normalize(output_format)
    if not isinstance(text, str):
        raise ValueError("text must be a string")
    input_bytes = len(text.encode("utf-8"))
    if input_bytes > _CONVERT_MAX_INPUT_BYTES:
        raise ValueError(
            f"format input exceeds maximum input size ({input_bytes} > "
            f"{_CONVERT_MAX_INPUT_BYTES} bytes)"
        )
    if source == "smiles":
        mol = from_smiles(text)
    elif source == "mol":
        mol = from_mol_block(text)
    elif source == "mol_v3000":
        mol = from_mol_v3000(text)
    elif source == "mol2":
        mol = from_mol2(text)
    elif source == "cml":
        mol = from_cml(text)
    elif source == "cjson":
        mol, parsed_coords = from_cjson(text)
        if coords is None and parsed_coords:
            coords = parsed_coords
    elif source == "moljson":
        mol = from_moljson(text)
    elif source == "cdxml":
        mol = from_cdxml(text)
    elif source == "pdb":
        mol, parsed_coords = from_pdb(text)
        if coords is None:
            coords = parsed_coords
    elif source == "xyz":
        mol, parsed_coords = from_xyz(text)
        if coords is None:
            coords = parsed_coords
    elif source == "pdbqt":
        mol = from_pdbqt(text)
    elif source == "gjf":
        mol = from_gjf(text)
    else:  # pragma: no cover
        raise ValueError(f"unsupported input format: {source}")

    if target == "smiles":
        return mol.smiles
    if target == "mol":
        return mol.to_mol_block()
    if target == "mol_v3000":
        return mol.to_mol_v3000([], name=name)
    if target == "mol2":
        return mol.to_mol2()
    if target == "cml":
        return mol.to_cml()
    if target == "cjson":
        return mol.to_cjson(coords or [])
    if target == "moljson":
        return mol.to_moljson()
    if target == "cdxml":
        return mol.to_cdxml()
    if target in {"pdb", "xyz", "pdbqt"}:
        if coords is None:
            raise ValueError(f"{target.upper()} output requires coords=[[x, y, z], ...]")
        if target == "pdb":
            return mol.to_pdb(coords)
        if target == "xyz":
            return mol.to_xyz(coords, comment)
        if charges is None:
            charges = [0.0] * len(coords)
        if len(charges) != len(coords):
            raise ValueError("charges must have the same length as coords")
        return mol.to_pdbqt([tuple(point) for point in coords], charges, name)
    raise ValueError(f"unsupported output format: {target}")


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
    return pd.DataFrame(bulk.descriptors(_materialize_smiles_batch(smiles)))


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



class ParseReport:
    """Result of parse_smiles_report(): mol + warnings instead of raising."""

    __slots__ = ("mol", "ok", "warnings", "error")

    def __init__(self, mol, warnings, error):
        self.mol = mol
        self.ok = mol is not None
        self.warnings: list = warnings
        self.error: str | None = error

    def __repr__(self):
        return f"ParseReport(ok={self.ok}, warnings={self.warnings}, error={self.error!r})"


def parse_smiles_report(smiles: str, *, strict: bool = False) -> "ParseReport":
    """Parse SMILES and return a ParseReport instead of raising on failure.

    Returns (mol, warnings) so batch pipelines can handle errors without
    try/except at every call site. Useful for WASM, MCP, and AI agent pipelines.

    Args:
        smiles: SMILES string to parse
        strict: if True, treat any warning as an error (mol=None)

    Returns:
        ParseReport with .ok, .mol, .warnings, .error

    Example::

        report = chematic.parse_smiles_report("C/C=C/C")
        if report.ok:
            print(report.mol.mw)
        for w in report.warnings:
            print(w)

        # batch
        reports = [chematic.parse_smiles_report(s) for s in smiles_list]
        mols = [r.mol for r in reports if r.ok]
    """
    warnings = []
    try:
        mol = from_smiles(smiles)
    except Exception as exc:
        return ParseReport(None, [], f"W003_PARSE_FAILED: {exc}")

    # Empty molecule (lenient parser swallowed the input without atoms)
    if mol.heavy_atoms == 0 and smiles.strip():
        return ParseReport(None, [], f"W003_PARSE_FAILED: no atoms parsed from {smiles!r}")

    # W002_DROPPED_STEREO — stereo specified but not fully resolved
    sc = mol.stereo_completeness
    if isinstance(sc, dict):
        specified = sc.get("specified", 0)
        unspecified = sc.get("unspecified", 0)
        if specified == 0 and unspecified > 0 and ("@" in smiles or "/" in smiles):
            warnings.append(f"W002_DROPPED_STEREO: {unspecified} stereocenters unresolved")
    elif ("@" in smiles or "/" in smiles) and sc == 0.0:
        warnings.append("W002_DROPPED_STEREO: stereo notation present but none resolved")

    # W001_UNUSUAL_VALENCE — atom with valence outside normal range
    try:
        per_atom = mol.implicit_hcount_per_atom()
        if any(h < 0 for h in per_atom):
            warnings.append("W001_UNUSUAL_VALENCE: one or more atoms have unusual valence")
    except Exception:
        pass

    if strict and warnings:
        return ParseReport(None, warnings, warnings[0])

    return ParseReport(mol, warnings, None)


# ---------------------------------------------------------------------------
# Compound screening
# ---------------------------------------------------------------------------

_SCREEN_PROFILES = {
    "druglike":  ["lipinski", "veber", "pains", "brenk", "qed", "sa_score"],
    "fragment":  ["ro3", "pains", "brenk"],
    "leadlike":  ["lead_like", "pains", "brenk", "qed"],
}


def screen(smiles, profile: str = "druglike", filters=None) -> list:
    """Screen compounds against a preset or custom filter profile.

    Parameters
    ----------
    smiles : list[str] or str
        One or more SMILES strings.
    profile : str
        Preset profile: "druglike" (default), "fragment", or "leadlike".
        Ignored when *filters* is provided.
    filters : list[str] or None
        Explicit list of filters to apply (overrides *profile*).
        Supported: "lipinski", "veber", "pains", "brenk", "egan",
        "ghose", "ro3", "lead_like", "reos", "mcf", "ames",
        "pfizer_3_75", "qed" (>= 0.5), "sa_score" (<= 3.5).

    Returns
    -------
    list[dict]
        One dict per input SMILES with fields:
        - smiles, valid, mw, logp, tpsa, hbd, hba, qed, sa_score
        - one ``<name>_pass`` bool per requested filter
        - overall_pass (True only when all filter passes are True)
    """
    if isinstance(smiles, str):
        smiles = [smiles]
    else:
        smiles = _materialize_smiles_batch(smiles)

    active = filters if filters is not None else _SCREEN_PROFILES.get(profile, _SCREEN_PROFILES["druglike"])

    _FILTER_ATTR = {
        "lipinski":   "lipinski_passes",
        "veber":      "veber_passes",
        "pains":      "pains_passes",
        "brenk":      "brenk_passes",
        "egan":       "egan_passes",
        "ghose":      "ghose_passes",
        "ro3":        "ro3_passes",
        "lead_like":  "lead_like_passes",
        "reos":       "reos_passes",
        "mcf":        "mcf_passes",
        "ames":       "ames_passes",
        "pfizer_3_75": "pfizer_3_75_passes",
    }

    mols = bulk.parse(smiles)
    results = []
    for smi, mol in zip(smiles, mols):
        if mol is None:
            results.append({"smiles": smi, "valid": False, "overall_pass": False})
            continue
        d = mol.descriptors()
        row = {
            "smiles":   smi,
            "valid":    True,
            "mw":       d.get("mw"),
            "logp":     d.get("logp"),
            "tpsa":     d.get("tpsa"),
            "hbd":      d.get("hbd"),
            "hba":      d.get("hba"),
            "qed":      d.get("qed"),
            "sa_score": d.get("sa_score"),
        }
        for f in active:
            if f in _FILTER_ATTR:
                row[f + "_pass"] = bool(getattr(mol, _FILTER_ATTR[f]))
            elif f == "qed":
                row["qed_pass"] = (d.get("qed") or 0.0) >= 0.5
            elif f == "sa_score":
                row["sa_score_pass"] = (d.get("sa_score") or 999.0) <= 3.5
        row["overall_pass"] = all(v for k, v in row.items() if k.endswith("_pass"))
        results.append(row)
    return results
