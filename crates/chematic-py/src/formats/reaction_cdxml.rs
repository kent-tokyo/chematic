use pyo3::exceptions::{PyMemoryError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyByteArray, PyBytes, PyString};

fn coordinates(source: &Bound<'_, PyAny>, counts: &[usize]) -> PyResult<Vec<Vec<[f64; 2]>>> {
    let mut rows = Vec::new();
    for row in source.try_iter()? {
        let row = row?;
        let expected = counts
            .get(rows.len())
            .ok_or_else(|| PyValueError::new_err("too many component coordinate rows"))?;
        if row.is_instance_of::<PyString>()
            || row.is_instance_of::<PyBytes>()
            || row.is_instance_of::<PyByteArray>()
        {
            return Err(PyTypeError::new_err(
                "coordinate row must contain numeric pairs",
            ));
        }
        let mut points = Vec::new();
        for point in row.try_iter()? {
            let point = point?;
            if points.len() == *expected {
                return Err(PyValueError::new_err("too many atom coordinates"));
            }
            points.push(point.extract::<[f64; 2]>()?);
        }
        if points.len() != *expected {
            return Err(PyValueError::new_err(
                "one coordinate pair is required per atom",
            ));
        }
        rows.push(points);
    }
    if rows.len() != counts.len() {
        return Err(PyValueError::new_err(
            "coordinate rows must match reaction components",
        ));
    }
    Ok(rows)
}

/// Write a bounded single-step CDXML reaction with explicit roles and arrow.
/// Coordinates follow reactants, agents, products; arrow is [tail, head], in
/// CDXML Y-down units. No layout inference. Non-aromatic, unmapped graphs only;
/// stereo, queries and nonterminal double bonds are refused. Bracket H must
/// match native normal-valence inference; isolated/bridging H is refused. Limits: 16 KiB
/// reaction SMILES, 32 total components, 256 atoms/component, 1024 total atoms.
#[pyfunction]
pub(super) fn write_reaction_cdxml_checked(
    reaction_smiles: &str,
    coords: &Bound<'_, PyAny>,
    arrow: [[f64; 2]; 2],
) -> PyResult<String> {
    let limits = chematic_rxn::ReactionParseLimits {
        max_input_bytes: 16384,
        max_components_per_side: 32,
        max_atoms_per_molecule: 256,
        max_bonds_per_molecule: 512,
    };
    let reaction = chematic_rxn::parse_reaction_with_limits(reaction_smiles, &limits)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let counts: Vec<_> = reaction
        .reactants
        .iter()
        .chain(&reaction.agents)
        .chain(&reaction.products)
        .map(|mol| mol.atom_count())
        .collect();
    if counts.len() > 32 || counts.iter().sum::<usize>() > 1024 {
        return Err(PyValueError::new_err(
            "reaction exceeds component/atom limits",
        ));
    }
    let coords = coordinates(coords, &counts)?;
    chematic_mol::write_reaction_cdxml_checked(&reaction, &coords, arrow).map_err(|e| match e {
        chematic_mol::ReactionCdxmlError::Allocation => PyMemoryError::new_err(e.to_string()),
        _ => PyValueError::new_err(e.to_string()),
    })
}
