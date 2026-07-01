//! Molecule format parsers/writers exposed as free `#[pyfunction]`s (SMILES/SDF/MOL/CML/CJSON/MolJSON/CDXML/MOL2/PDBQT/GJF/CIF/InChI/PDB/XYZ/RXN) plus small serialization helpers.

use crate::Mol;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::sync::Arc;

pub(crate) fn bitvec2048_to_bytes(fp: &chematic_fp::bitvec::BitVec2048) -> Vec<u8> {
    (0..256usize)
        .map(|byte_idx| {
            let mut byte = 0u8;
            for bit in 0..8usize {
                if fp.get(byte_idx * 8 + bit) {
                    byte |= 1 << bit;
                }
            }
            byte
        })
        .collect()
}

pub(crate) fn flat_to_coords3d(coords: &[[f64; 3]]) -> chematic_3d::Coords3D {
    chematic_3d::Coords3D {
        points: coords
            .iter()
            .map(|c| chematic_3d::Point3::new(c[0], c[1], c[2]))
            .collect(),
    }
}

// ---------------------------------------------------------------------------
// Helper: convert MoleculeReport to Python dict
// ---------------------------------------------------------------------------

/// Parse a SMILES string and return a Mol.
///
/// Raises ``ValueError`` on invalid SMILES.
#[pyfunction]
fn from_smiles(smiles: &str) -> PyResult<Mol> {
    chematic_smiles::parse(smiles)
        .map(|mol| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a CXSMILES string and return the molecule with CX metadata.
///
/// Returns a 2-tuple ``(mol, cx)`` where ``cx`` is a dict with:
///
/// - ``atom_labels``: list of atom label strings (or ``None`` per atom)
/// - ``atom_props``: list of ``{"atom_idx", "key", "value"}`` dicts
/// - ``atom_radicals``: list of radical class integers (or ``None`` per atom)
///
/// Raises ``ValueError`` on parse failure. CXSMILES without a CX extension
/// block behaves like :func:`from_smiles` (all CX fields are empty).
///
///     mol, cx = chematic.from_cxsmiles("CC |$R1;R2$|")
///     print(cx['atom_labels'])   # ['R1', 'R2']
#[pyfunction]
fn from_cxsmiles<'py>(py: Python<'py>, s: &str) -> PyResult<(Mol, Bound<'py, PyDict>)> {
    let cx =
        chematic_smiles::parse_cxsmiles(s).map_err(|e| PyValueError::new_err(e.to_string()))?;

    let d = PyDict::new(py);
    let labels: Vec<Option<&str>> = cx.atom_labels.iter().map(|l| l.as_deref()).collect();
    d.set_item("atom_labels", labels)?;
    let props: Vec<Bound<'py, PyDict>> = cx
        .atom_props
        .iter()
        .map(|p| {
            let pd = PyDict::new(py);
            pd.set_item("atom_idx", p.atom.0 as usize).unwrap();
            pd.set_item("key", &p.key).unwrap();
            pd.set_item("value", &p.value).unwrap();
            pd
        })
        .collect();
    d.set_item("atom_props", props)?;
    let radicals: Vec<Option<u8>> = cx.atom_radicals.clone();
    d.set_item("atom_radicals", radicals)?;

    Ok((
        Mol {
            inner: Arc::new(cx.mol),
            props: Default::default(),
        },
        d,
    ))
}

/// Parse a MOL/SDF block and return a Mol.
///
/// Parse a condensed molecular formula (e.g., ``"CH3OH"``, ``"C6H12O6"``) into a :class:`Mol`.
///
/// Returns ``None`` if the formula is unknown or ambiguous.
/// Unlike SMILES parsing, condensed formulas may not encode connectivity uniquely;
/// this function uses a built-in formula→structure dictionary.
///
/// Equivalent to chempy's condensed formula support.
///
///     mol = chematic.from_condensed("CH3OH")  # methanol
///     if mol:
///         print(mol.smiles)  # CO
///
///     mol = chematic.from_condensed("C6H12O6")  # glucose
#[pyfunction]
fn from_condensed(formula: &str) -> Option<Mol> {
    chematic_chem::parse_condensed(formula).ok().map(|mol| Mol {
        inner: Arc::new(mol),
        props: Default::default(),
    })
}

/// Raises ``ValueError`` on parse failure.
#[pyfunction]
fn from_mol_block(block: &str) -> PyResult<Mol> {
    chematic_mol::parse_mol(block)
        .map(|(mol, _meta)| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a MDL MOL V2000 block and return the molecule with its 2D layout coordinates.
///
/// Returns a 3-tuple ``(mol, name, coords_2d)`` where:
///
/// - ``mol``: :class:`Mol` object
/// - ``name``: molecule name from the MOL header (may be empty)
/// - ``coords_2d``: list of ``[x, y]`` pairs (one per heavy atom, Å)
///
/// Raises ``ValueError`` on parse failure.
///
/// Use :func:`from_mol_block` if you only need the molecule graph.
/// Use this function when you want to preserve the 2D layout for display or
/// round-trip back to MOL format via :meth:`Mol.to_mol_block_2d`.
///
///     mol, name, coords_2d = chematic.from_mol_block_with_coords(block)
///     new_block = mol.to_mol_block_2d(coords_2d, name=name)
#[pyfunction]
fn from_mol_block_with_coords(block: &str) -> PyResult<(Mol, String, Vec<Vec<f64>>)> {
    chematic_mol::parse_mol_with_coords(block)
        .map(|(mol, meta, coords)| {
            let py_coords: Vec<Vec<f64>> = coords.iter().map(|(x, y)| vec![*x, *y]).collect();
            (
                Mol {
                    inner: Arc::new(mol),
                    props: Default::default(),
                },
                meta.name,
                py_coords,
            )
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a multi-record SDF string and return all molecules with their 2D layout coordinates.
///
/// Returns a list of 3-tuples ``(mol, name, coords_2d)`` — one per SDF record.
/// Invalid records are silently skipped (same behaviour as :func:`iter_sdf`).
///
/// This is the batch equivalent of :func:`from_mol_block_with_coords`.
///
///     with open("library.sdf") as f:
///         records = chematic.parse_sdf_with_coords(f.read())
///     for mol, name, coords_2d in records:
///         new_block = mol.to_mol_block_2d(coords_2d, name=name)
#[pyfunction]
fn parse_sdf_with_coords(text: &str) -> Vec<(Mol, String, Vec<Vec<f64>>)> {
    // Split SDF by $$$$ delimiter and parse each block with parse_mol_with_coords.
    // This avoids the Rust parse_sdf_with_coords leading-blank-line stripping issue.
    let mut results = Vec::new();
    let mut remaining = text;
    loop {
        let (block, rest) = match remaining.find("$$$$") {
            Some(pos) => {
                let after = &remaining[pos + 4..];
                let after = after
                    .strip_prefix("\r\n")
                    .or_else(|| after.strip_prefix('\n'))
                    .unwrap_or(after);
                (&remaining[..pos], after)
            }
            None => (remaining, ""),
        };
        if !block.trim().is_empty() {
            if let Ok((mol, meta, coords)) = chematic_mol::parse_mol_with_coords(block) {
                let py_coords: Vec<Vec<f64>> = coords.iter().map(|(x, y)| vec![*x, *y]).collect();
                results.push((
                    Mol {
                        inner: Arc::new(mol),
                        props: Default::default(),
                    },
                    meta.name,
                    py_coords,
                ));
            }
        }
        if remaining.find("$$$$").is_none() || rest.is_empty() {
            break;
        }
        remaining = rest;
    }
    results
}

/// Parse a Chemical Markup Language (CML) string into a ``Mol`` object.
///
/// Raises ``ValueError`` on parse failure.
///
///     with open("molecule.cml") as f:
///         mol = chematic.from_cml(f.read())
#[pyfunction]
fn from_cml(cml_str: &str) -> PyResult<Mol> {
    chematic_mol::parse_cml(cml_str)
        .map(|(mol, _coords)| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a ChemicalJSON (.cjson) string.
///
/// Returns ``(mol, coords)`` where ``coords`` is a list of ``[x, y, z]``
/// coordinate triples (Å), one per heavy atom.  ``coords`` is empty when
/// the file has no ``atoms.coords.3d`` field.
///
/// ChemicalJSON is the native format of Avogadro 2 and the MolSSI
/// Open Chemistry toolkit.
///
/// Raises ``ValueError`` on parse failure.
///
///     mol, coords = chematic.from_cjson(open("mol.cjson").read())
///     print(mol.smiles)
///     # Round-trip:
///     open("out.cjson", "w").write(mol.to_cjson(coords))
#[pyfunction]
fn from_cjson(cjson_str: &str) -> PyResult<(Mol, Vec<Vec<f64>>)> {
    let (mol, coords) =
        chematic_mol::parse_cjson(cjson_str).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let py_coords = coords.iter().map(|&(x, y, z)| vec![x, y, z]).collect();
    Ok((
        Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        },
        py_coords,
    ))
}

/// Parse a MolJSON string into a ``Mol`` object.
///
/// MolJSON is a JSON-based molecular representation designed for LLM
/// (large language model) compatibility.
///
/// Raises ``ValueError`` on parse failure.
///
///     mol = chematic.from_moljson(open("mol.json").read())
///     # Round-trip:
///     json_str = mol.to_moljson()
///     mol2 = chematic.from_moljson(json_str)
#[pyfunction]
fn from_moljson(json_str: &str) -> PyResult<Mol> {
    chematic_mol::parse_moljson(json_str)
        .map(|mol| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a ChemDraw XML (CDXML) string into a ``Mol`` object.
///
/// Raises ``ValueError`` on parse failure.
///
///     with open("molecule.cdxml") as f:
///         mol = chematic.from_cdxml(f.read())
#[pyfunction]
fn from_cdxml(cdxml_str: &str) -> PyResult<Mol> {
    chematic_mol::parse_cdxml(cdxml_str)
        .map(|(mol, _coords)| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse an MDL MOL V3000 (``V3000``) block into a ``Mol`` object.
///
/// Raises ``ValueError`` on parse failure.
///
///     with open("ligand_v3000.mol") as f:
///         mol = chematic.from_mol_v3000(f.read())
#[pyfunction]
fn from_mol_v3000(block: &str) -> PyResult<Mol> {
    chematic_mol::parse_mol_v3000(block)
        .map(|(mol, _meta)| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a MDL MOL V3000 block and return the molecule with its 2D layout coordinates.
///
/// Returns a 3-tuple ``(mol, name, coords_2d)`` identical to
/// :func:`from_mol_block_with_coords` but for V3000 input.
///
/// Raises ``ValueError`` on parse failure.
///
///     mol, name, coords_2d = chematic.from_mol_v3000_with_coords(block)
///     new_block = mol.to_mol_v3000(coords_2d, name=name)
#[pyfunction]
fn from_mol_v3000_with_coords(block: &str) -> PyResult<(Mol, String, Vec<Vec<f64>>)> {
    chematic_mol::parse_mol_v3000_with_coords(block)
        .map(|(mol, meta, coords)| {
            let py_coords: Vec<Vec<f64>> = coords.iter().map(|(x, y)| vec![*x, *y]).collect();
            (
                Mol {
                    inner: Arc::new(mol),
                    props: Default::default(),
                },
                meta.name,
                py_coords,
            )
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a Tripos MOL2 string into a ``Mol`` object.
///
/// Example::
///
///     with open("ligand.mol2") as f:
///         mol = chematic.from_mol2(f.read())
///     print(mol.mw)
///
/// Raises ``ValueError`` on parse failure.
#[pyfunction]
fn from_mol2(mol2_str: &str) -> PyResult<Mol> {
    chematic_mol::parse_mol2(mol2_str)
        .map(|(mol, _coords)| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse an AutoDock PDBQT string and return a :class:`Mol`.
///
/// Only the molecular graph (elements and bonds) is extracted; 3D coordinates
/// and partial charges are discarded.  To retain them, use the lower-level
/// :func:`chematic_mol.parse_pdbqt` Rust API directly.
///
/// Raises:
///     ValueError: on parse failure.
///
/// Example::
///
///     with open("ligand.pdbqt") as f:
///         mol = chematic.from_pdbqt(f.read())
///     print(mol.mw)
#[pyfunction]
fn from_pdbqt(pdbqt_str: &str) -> PyResult<Mol> {
    chematic_mol::parse_pdbqt(pdbqt_str)
        .map(|(mol, _coords, _charges)| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a Gaussian input file (`.gjf` / `.com`) and return a :class:`Mol`.
///
/// Raises:
///     ValueError: on parse failure.
///
/// Example::
///
///     mol = chematic.from_gjf(open("mol.gjf").read())
#[pyfunction]
fn from_gjf(gjf_str: &str) -> PyResult<Mol> {
    chematic_mol::parse_gjf(gjf_str)
        .map(|(mol, _coords, _charge, _mult)| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a Gaussian output file (`.log` / `.out`) and return a dict with
/// ``mol``, ``coords`` and ``scf_energy`` fields.
///
/// Returns:
///     dict: ``{"mol": Mol, "coords": list[list[float]], "scf_energy": float | None}``
///
/// Raises:
///     ValueError: when no `Standard orientation:` block is found.
#[pyfunction]
fn parse_gaussian_log<'py>(
    py: Python<'py>,
    log_str: &str,
) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    let result = chematic_mol::parse_gaussian_log(log_str)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let mol = Mol {
        inner: Arc::new(result.mol),
        props: Default::default(),
    };
    let coords: Vec<Vec<f64>> = result
        .coords
        .iter()
        .map(|&(x, y, z)| vec![x, y, z])
        .collect();
    let d = pyo3::types::PyDict::new(py);
    d.set_item("mol", mol)?;
    d.set_item("coords", coords)?;
    d.set_item("scf_energy", result.scf_energy)?;
    Ok(d)
}

/// Generate a Gaussian input file (`.gjf`) string from a molecule.
///
/// Args:
///     mol: The molecule to write.
///     coords: Atomic coordinates as ``[[x, y, z], ...]`` in Ångströms.
///     charge: Formal charge (default 0).
///     multiplicity: Spin multiplicity (default 1).
///     method: Route section keywords (default ``"B3LYP/6-31G* opt"``).
///     title: Job title comment (default ``"chematic"``).
///
/// Returns:
///     str: GJF file contents.
#[pyfunction]
#[pyo3(signature = (mol, coords, charge=0, multiplicity=1, method="B3LYP/6-31G* opt", title="chematic"))]
fn write_gjf(
    mol: &Mol,
    coords: Vec<[f64; 3]>,
    charge: i32,
    multiplicity: u32,
    method: &str,
    title: &str,
) -> String {
    let c: Vec<(f64, f64, f64)> = coords.into_iter().map(|[x, y, z]| (x, y, z)).collect();
    chematic_mol::write_gjf(&mol.inner, &c, charge, multiplicity, method, title)
}

/// Parse a CIF (Crystallographic Information File) string and return a dict.
///
/// Returns:
///     dict: ``{"mol": Mol, "coords": list[list[float]], "cell": dict | None}``
///     where ``cell`` has keys ``a, b, c, alpha, beta, gamma``.
///
/// Raises:
///     ValueError: on parse failure.
///
/// Example::
///
///     result = chematic.parse_cif(open("structure.cif").read())
///     mol = result["mol"]
///     coords = result["coords"]
#[pyfunction]
fn parse_cif<'py>(py: Python<'py>, cif_str: &str) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    let result =
        chematic_mol::parse_cif(cif_str).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let mol = Mol {
        inner: Arc::new(result.mol),
        props: Default::default(),
    };
    let coords: Vec<Vec<f64>> = result
        .coords
        .iter()
        .map(|&(x, y, z)| vec![x, y, z])
        .collect();
    let d = pyo3::types::PyDict::new(py);
    d.set_item("mol", mol)?;
    d.set_item("coords", coords)?;
    if let Some(cell) = result.cell {
        let cd = pyo3::types::PyDict::new(py);
        cd.set_item("a", cell.a)?;
        cd.set_item("b", cell.b)?;
        cd.set_item("c", cell.c)?;
        cd.set_item("alpha", cell.alpha)?;
        cd.set_item("beta", cell.beta)?;
        cd.set_item("gamma", cell.gamma)?;
        d.set_item("cell", cd)?;
    } else {
        d.set_item("cell", py.None())?;
    }
    Ok(d)
}

/// Return True if the SMILES can be parsed without error.
#[pyfunction]
fn is_valid_smiles(smiles: &str) -> bool {
    chematic_smiles::parse(smiles).is_ok()
}

/// Return ``True`` if ``smarts`` is a valid SMARTS pattern, ``False`` otherwise.
///
/// Mirrors :func:`is_valid_smiles` for SMARTS pattern validation.
/// Useful for validating user-supplied SMARTS before calling
/// :func:`smarts_match` or :func:`smarts_find`.
///
///     chematic.is_valid_smarts("c1ccccc1")  # True
///     chematic.is_valid_smarts("[invalid")  # False
///     chematic.is_valid_smarts("[#6]-[#7]") # True
#[pyfunction]
fn is_valid_smarts(smarts: &str) -> bool {
    chematic_smarts::parse_smarts(smarts).is_ok()
}

/// Parse an InChI string and return a Mol.
///
/// Raises ``ValueError`` on parse failure.
///
///     mol = chematic.from_inchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
#[pyfunction]
fn from_inchi(inchi: &str) -> PyResult<Mol> {
    chematic_inchi::parse_inchi(inchi)
        .map(|mol| Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse a PDB string and return ``(Mol, coords)`` where coords is a list of ``[x,y,z]``.
///
/// Bond information is inferred from inter-atom distances; only ATOM/HETATM records
/// are used.  Returns ``ValueError`` when no atoms are found.
///
///     mol, coords = chematic.from_pdb(open("ligand.pdb").read())
#[pyfunction]
fn from_pdb(pdb_str: &str) -> PyResult<(Mol, Vec<Vec<f64>>)> {
    let atoms = chematic_3d::parse_pdb_atoms(pdb_str);
    if atoms.is_empty() {
        return Err(PyValueError::new_err(
            "no ATOM/HETATM records found in PDB input",
        ));
    }
    let (mol, c3d) = chematic_3d::pdb_to_molecule(&atoms);
    let coords = c3d.points.iter().map(|p| vec![p.x, p.y, p.z]).collect();
    Ok((
        Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        },
        coords,
    ))
}

/// Parse an XYZ string and return ``(Mol, coords)`` where coords is a list of ``[x,y,z]``.
///
/// Bond information is inferred from inter-atom distances.
///
///     mol, coords = chematic.from_xyz(open("molecule.xyz").read())
#[pyfunction]
fn from_xyz(xyz_str: &str) -> PyResult<(Mol, Vec<Vec<f64>>)> {
    let (mol, c3d) =
        chematic_3d::parse_xyz(xyz_str).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let coords = c3d.points.iter().map(|p| vec![p.x, p.y, p.z]).collect();
    Ok((
        Mol {
            inner: Arc::new(mol),
            props: Default::default(),
        },
        coords,
    ))
}

/// Parse a Hill-notation molecular formula string into an element count dictionary.
///
/// Returns a ``dict[str, int]`` mapping element symbol → atom count.
/// Mirrors the API of PyPI libraries **chemparse** and **chemformula**.
///
/// Supported syntax:
///   - Simple formulas: ``"H2O"``, ``"C6H12O6"``
///   - Parentheses with multipliers: ``"Ca(OH)2"`` → ``{"Ca":1,"O":2,"H":2}``
///   - SMILES-style brackets: ``"[NH4]+"`` → ``{"N":1,"H":4}``
///   - Trailing charge signs are ignored: ``"NH4+"`` → same as ``"NH4"``
///
/// Raises:
///     ValueError: on empty formula or unbalanced parentheses.
///
///     chematic.parse_formula("C6H12O6")  # {"C": 6, "H": 12, "O": 6}
///     chematic.parse_formula("Ca(OH)2")  # {"Ca": 1, "O": 2, "H": 2}
///     chematic.parse_formula("[NH4]+")   # {"N": 1, "H": 4}
#[pyfunction]
fn parse_formula<'py>(formula: &str, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let counts =
        chematic_chem::parse_formula(formula).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let d = PyDict::new(py);
    let mut sorted: Vec<(String, u32)> = counts.into_iter().collect();
    sorted.sort_by(|a, b| a.0.cmp(&b.0));
    for (elem, cnt) in sorted {
        d.set_item(elem, cnt)?;
    }
    Ok(d)
}

/// Compute atom economy of a reaction (green chemistry metric).
/// E-factor (Environmental Factor) — waste-to-product mass ratio.
///
/// E-factor = waste_mass / product_mass.  Lower is greener.
/// Fine chemicals typically E=5–50; pharmaceuticals E=25–100.
///
///     ef = chematic.e_factor(waste_kg=90.0, product_kg=10.0)  # → 9.0
/// Fast structural hash for deduplication — one int per molecule.
///
/// Molecules with the same canonical graph return identical hashes.
/// Use with :func:`are_identical` to confirm true equivalence (no hash collisions).
///
/// Equivalent to RDKit's ``rdMolHash.MolHash()``.
///
///     seen = set()
///     unique = [m for m in mols if (h := chematic.mol_hash(m)) not in seen and not seen.add(h)]
#[pyfunction]
fn mol_hash(mol: &Mol) -> u64 {
    chematic_chem::mol_hash(&mol.inner)
}

/// Check whether two molecules are graph-isomorphic (exact structural identity).
///
/// More reliable than comparing SMILES strings (which depend on canonicalization).
/// Equivalent to RDKit's ``Chem.MolToInchiKey(m1) == Chem.MolToInchiKey(m2)``.
///
///     assert chematic.are_identical(
///         chematic.from_smiles("c1ccccc1"),
///         chematic.from_smiles("C1=CC=CC=C1"),  # kekulé form
///     )
#[pyfunction]
fn are_identical(mol1: &Mol, mol2: &Mol) -> bool {
    chematic_chem::are_identical(&mol1.inner, &mol2.inner)
}

/// Normalize and re-serialize a reaction SMILES.
///
/// Parses the reaction SMILES and writes it back in canonical form.
/// Useful for standardizing reaction data before storing or comparing.
///
///     canon = chematic.write_reaction("CC(=O)Cl.[NH3]>>CC(=O)N.HCl")
#[pyfunction]
fn write_reaction(reaction_smiles: &str) -> PyResult<String> {
    let rxn = chematic_rxn::parse_reaction(reaction_smiles)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(chematic_rxn::write_reaction(&rxn))
}

///
/// atom_economy = MW(desired products) / MW(all reactants) × 100.
/// Parse a MDL RXN V2000 file and return the canonical reaction SMILES.
///
/// Raises ``ValueError`` on parse failure.
///
///     rxn_smiles = chematic.from_rxn_file(text)
///     ae = chematic.atom_economy(rxn_smiles)
#[pyfunction]
fn from_rxn_file(text: &str) -> PyResult<String> {
    let rxn =
        chematic_mol::parse_rxn_file(text).map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(chematic_rxn::write_reaction(&rxn))
}

/// Convert a reaction SMILES string to MDL RXN V2000 format.
///
/// Raises ``ValueError`` on invalid reaction SMILES.
///
///     block = chematic.to_rxn_file("CC(=O)Cl.[NH3]>>CC(=O)N.HCl")
#[pyfunction]
fn to_rxn_file(reaction_smiles: &str) -> PyResult<String> {
    let rxn = chematic_rxn::parse_reaction(reaction_smiles)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(chematic_mol::write_rxn_file(&rxn))
}

/// Find top-K nearest neighbors from precomputed fingerprint byte arrays.
///
/// More efficient than :func:`top_k_similar_fp` when the same ``db_fps`` list is
/// reused across multiple queries (fingerprints computed only once).
///
///     db_fps = [mol.ecfp4() for mol in library]   # compute once
///     for query in queries:
///         hits = chematic.nearest_neighbors_from_fp(query.ecfp4(), db_fps, k=10)
///         for idx, score in hits:
///             print(library_smiles[idx], score)
#[pyfunction]
#[pyo3(signature = (query_fp, db_fps, k = 10))]
fn nearest_neighbors_from_fp(query_fp: &[u8], db_fps: Vec<Vec<u8>>, k: usize) -> Vec<(usize, f64)> {
    let qa: u32 = query_fp.iter().map(|b| b.count_ones()).sum();
    let mut scores: Vec<(usize, f64)> = db_fps
        .iter()
        .enumerate()
        .filter_map(|(i, fp)| {
            if fp.len() != query_fp.len() {
                return None;
            }
            let and: u32 = query_fp
                .iter()
                .zip(fp.iter())
                .map(|(a, b)| (a & b).count_ones())
                .sum();
            let db_cnt: u32 = fp.iter().map(|b| b.count_ones()).sum();
            let or = qa + db_cnt - and;
            if or == 0 {
                return None;
            }
            Some((i, and as f64 / or as f64))
        })
        .collect();
    scores.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    scores.truncate(k);
    scores
}

/// Parse a ``.smi`` file (tab/space-separated SMILES + name) into (Mol, name) pairs.
///
/// Each line is ``SMILES[<tab>name]``. Lines with invalid SMILES are silently skipped.
/// Comment lines starting with ``#`` and blank lines are ignored.
/// Equivalent to RDKit's ``Chem.SmilesMolSupplier``.
///
///     records = chematic.parse_smi_file(open("library.smi").read())
///     for mol, name in records:
///         print(name, mol.mw)
#[pyfunction]
fn parse_smi_file(content: &str) -> Vec<(Mol, String)> {
    chematic_smiles::parse_smi_file(content)
        .into_iter()
        .filter_map(|r| r.ok())
        .map(|(mol, name)| {
            (
                Mol {
                    inner: Arc::new(mol),
                    props: Default::default(),
                },
                name,
            )
        })
        .collect()
}

/// Write (Mol, name) pairs to ``.smi`` format.
///
/// Output format: ``SMILES<TAB>name<NEWLINE>`` per record (name omitted if empty).
/// Equivalent to RDKit's ``Chem.SmilesWriter``.
///
///     text = chematic.write_smi_file([(mol1, "cpd1"), (mol2, "cpd2")])
///     with open("output.smi", "w") as f:
///         f.write(text)
#[pyfunction]
fn write_smi_file(records: Vec<(Mol, String)>) -> String {
    let mut out = String::new();
    for (mol, name) in &records {
        let smiles = chematic_smiles::canonical_smiles(&mol.inner);
        if name.is_empty() {
            out.push_str(&smiles);
        } else {
            out.push_str(&smiles);
            out.push('\t');
            out.push_str(name);
        }
        out.push('\n');
    }
    out
}

// ---------------------------------------------------------------------------
// Register
// ---------------------------------------------------------------------------

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(from_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(from_cxsmiles, m)?)?;
    m.add_function(wrap_pyfunction!(from_condensed, m)?)?;
    m.add_function(wrap_pyfunction!(from_mol_block, m)?)?;
    m.add_function(wrap_pyfunction!(from_mol_block_with_coords, m)?)?;
    m.add_function(wrap_pyfunction!(parse_sdf_with_coords, m)?)?;
    m.add_function(wrap_pyfunction!(from_cml, m)?)?;
    m.add_function(wrap_pyfunction!(from_cjson, m)?)?;
    m.add_function(wrap_pyfunction!(from_moljson, m)?)?;
    m.add_function(wrap_pyfunction!(from_cdxml, m)?)?;
    m.add_function(wrap_pyfunction!(from_mol_v3000, m)?)?;
    m.add_function(wrap_pyfunction!(from_mol_v3000_with_coords, m)?)?;
    m.add_function(wrap_pyfunction!(from_mol2, m)?)?;
    m.add_function(wrap_pyfunction!(from_pdbqt, m)?)?;
    m.add_function(wrap_pyfunction!(from_gjf, m)?)?;
    m.add_function(wrap_pyfunction!(parse_gaussian_log, m)?)?;
    m.add_function(wrap_pyfunction!(write_gjf, m)?)?;
    m.add_function(wrap_pyfunction!(parse_cif, m)?)?;
    m.add_function(wrap_pyfunction!(is_valid_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(is_valid_smarts, m)?)?;
    m.add_function(wrap_pyfunction!(from_inchi, m)?)?;
    m.add_function(wrap_pyfunction!(from_pdb, m)?)?;
    m.add_function(wrap_pyfunction!(from_xyz, m)?)?;
    m.add_function(wrap_pyfunction!(parse_formula, m)?)?;
    m.add_function(wrap_pyfunction!(mol_hash, m)?)?;
    m.add_function(wrap_pyfunction!(are_identical, m)?)?;
    m.add_function(wrap_pyfunction!(write_reaction, m)?)?;
    m.add_function(wrap_pyfunction!(from_rxn_file, m)?)?;
    m.add_function(wrap_pyfunction!(to_rxn_file, m)?)?;
    m.add_function(wrap_pyfunction!(nearest_neighbors_from_fp, m)?)?;
    m.add_function(wrap_pyfunction!(parse_smi_file, m)?)?;
    m.add_function(wrap_pyfunction!(write_smi_file, m)?)?;
    Ok(())
}
