//! RDKit's `MurckoScaffold.GetScaffoldForMol` (`MolOps::MurckoDecompose`,
//! `Code/GraphMol/ChemTransforms/ChemTransforms.cpp`) on the RDKit-model
//! molecule.

use super::RdkitSmilesError;
use super::mol::{BondType, ChiralTag, Mol};

/// `MurckoDecompose(mol)` followed by `GetScaffoldForMol`'s
/// `UpdatePropertyCache()` and `GetSymmSSSR`. The result's stereo is not
/// perceived (`_StereochemDone` is cleared).
pub(crate) fn murcko_decompose(mol: &Mol) -> Result<Mol, RdkitSmilesError> {
    let n = mol.atoms.len();
    let mut res = mol.clone();
    if n == 0 {
        return Ok(res);
    }
    let ri = mol.ring_info();
    let mut keep: Vec<bool> = (0..n).map(|i| ri.num_atom_rings(i) > 0).collect();
    // Shortest paths between the first atoms of every pair of rings. A
    // non-ring atom on such a path lies on the unique linker between two
    // ring systems, so any shortest path marks the same atoms.
    let rings = &ri.atom_rings;
    for i in 0..rings.len() {
        for j in (i + 1)..rings.len() {
            let (from, to) = (rings[i][0], rings[j][0]);
            // BFS from `to`: parent[x] = next atom from x toward `to`.
            let mut parent = vec![usize::MAX; n];
            parent[to] = to;
            let mut queue = std::collections::VecDeque::from([to]);
            while let Some(a) = queue.pop_front() {
                for b in mol.nbrs(a) {
                    if parent[b] == usize::MAX {
                        parent[b] = a;
                        queue.push_back(b);
                    }
                }
            }
            let mut at = from;
            while at != to {
                keep[at] = true;
                at = parent[at];
                if at == usize::MAX {
                    break;
                }
            }
        }
    }
    let mut removed = Vec::new();
    for i in 0..n {
        if keep[i] {
            continue;
        }
        let mut remove_it = true;
        let bonds = res.atom_bonds[i].clone();
        for b in bonds {
            let nbr = res.bonds[b].other(i);
            if !keep[nbr] {
                continue;
            }
            let atom = &mut res.atoms[nbr];
            if res.bonds[b].bt == BondType::Double {
                remove_it = false;
                break;
            } else if atom.aromatic && atom.anum != 6 {
                atom.num_explicit_hs = 1;
            } else if atom.aromatic && atom.anum == 6 && atom.charge == 1 {
                atom.num_explicit_hs = 1;
            } else if atom.no_implicit || atom.chiral != ChiralTag::Unspecified {
                atom.no_implicit = false;
                atom.num_explicit_hs = 0;
                atom.chiral = ChiralTag::Unspecified;
            }
        }
        if remove_it {
            removed.push(i);
        }
    }
    // commitBatchEdit: removing a bond clears the stereo atoms that refer
    // across it on the neighbouring bonds.
    let mut dead = vec![false; res.bonds.len()];
    for &a in &removed {
        for &b in &res.atom_bonds[a] {
            dead[b] = true;
        }
    }
    for b in 0..res.bonds.len() {
        if !dead[b] {
            continue;
        }
        let (u, v) = (res.bonds[b].begin, res.bonds[b].end);
        for (x, y) in [(u, v), (v, u)] {
            for &ob in &res.atom_bonds[x].clone() {
                if ob == b || !res.bonds[ob].stereo_atoms.contains(&y) {
                    continue;
                }
                // (The model has no STEREOCIS/STEREOTRANS, which RDKit
                // would also reset to STEREONONE here.)
                res.bonds[ob].stereo_atoms.clear();
            }
        }
    }
    for &a in removed.iter().rev() {
        res.remove_atom(a);
    }
    for a in &mut res.atoms {
        a.cip_code = None;
        a.chirality_possible = false;
    }
    res.update_property_cache(true)?;
    res.find_rings()?;
    Ok(res)
}
