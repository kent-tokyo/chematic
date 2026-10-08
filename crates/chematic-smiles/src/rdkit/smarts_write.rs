//! `SmartsWrite` (RDKit 2026.03.1 `SmartsWrite.cpp`) for molecules without
//! query features: `MolToSmarts` and `MolToCXSmarts` of a molecule read with
//! `MolFromSmiles`.
//!
//! RDKit writes such a molecule with atom indices as ranks (no
//! canonicalization) and ring information cleared, every atom as a bracket
//! `[#n]`/`[Sym]` atom and every bond explicitly.

use super::RdkitSmilesError;
use super::canon::{StackElem, canonicalize_fragment};
use super::mol::{BondDir, BondType, ChiralTag, Mol};
use super::periodic;
use super::write::{fragment, mol_frags};

/// `SmilesWrite::inOrganicSubset`.
fn in_organic_subset(anum: u32) -> bool {
    matches!(anum, 0 | 5 | 6 | 7 | 8 | 9 | 15 | 16 | 17 | 35 | 53)
}

/// `getNonQueryAtomSmarts`.
fn atom_smarts(mol: &Mol, a: usize, isomeric: bool) -> String {
    let atom = &mol.atoms[a];
    let mut res = String::from("[");
    if atom.isotope != 0 {
        res.push_str(&atom.isotope.to_string());
    }
    if in_organic_subset(atom.anum) {
        res.push('#');
        res.push_str(&atom.anum.to_string());
    } else {
        res.push_str(periodic::symbol(atom.anum));
    }
    let mut added_chirality = false;
    if isomeric {
        match atom.chiral {
            ChiralTag::Cw => {
                res.push_str("@@");
                added_chirality = true;
            }
            ChiralTag::Ccw => {
                res.push('@');
                added_chirality = true;
            }
            ChiralTag::Unspecified => {}
        }
    }
    if added_chirality && atom.num_explicit_hs == 1 {
        res.push('H');
    }
    match atom.charge {
        0 => {}
        -1 => res.push('-'),
        1 => res.push('+'),
        c if c < 0 => res.push_str(&c.to_string()),
        c => {
            res.push('+');
            res.push_str(&c.to_string());
        }
    }
    if let Some(m) = atom.map {
        res.push(':');
        res.push_str(&m.to_string());
    }
    res.push(']');
    res
}

/// `getNonQueryBondSmarts`.
fn bond_smarts(mol: &Mol, b: usize, atom_to_left: usize, isomeric: bool, dative: bool) -> String {
    let bond = &mol.bonds[b];
    let dir = |default: &str| -> String {
        if isomeric {
            match bond.dir {
                BondDir::EndDownRight => return "\\".into(),
                BondDir::EndUpRight => return "/".into(),
                BondDir::None => {}
            }
        }
        default.into()
    };
    if bond.aromatic {
        return dir(":");
    }
    match bond.bt {
        BondType::Single => dir("-"),
        BondType::Double => "=".into(),
        BondType::Triple => "#".into(),
        BondType::Quadruple => "$".into(),
        BondType::Aromatic => dir(":"),
        BondType::Dative => {
            if !dative {
                "-".into()
            } else if bond.begin != atom_to_left {
                "<-".into()
            } else {
                "->".into()
            }
        }
    }
}

/// `MolToSmarts(mol, params)` with `isomericSmiles`, `rootedAtAtom` and
/// `includeDativeBonds`; also returns the atom and bond output orders
/// (`_smilesAtomOutputOrder`, `_smilesBondOutputOrder`).
pub(crate) fn mol_to_smarts(
    mol: &Mol,
    isomeric: bool,
    rooted_at_atom: Option<usize>,
    dative: bool,
) -> Result<(String, Vec<usize>, Vec<usize>), RdkitSmilesError> {
    let n = mol.atoms.len();
    if n == 0 {
        return Ok((String::new(), Vec::new(), Vec::new()));
    }
    if let Some(r) = rooted_at_atom
        && r >= n
    {
        return Err(RdkitSmilesError::Unsupported("bad atom index".into()));
    }
    let frags = mol_frags(mol);
    let mut frag_of = vec![0usize; n];
    for (f, atoms) in frags.iter().enumerate() {
        for &a in atoms {
            frag_of[a] = f;
        }
    }
    let mut done = vec![false; frags.len()];
    let mut res = String::new();
    let mut atom_order = Vec::with_capacity(n);
    let mut bond_order = Vec::with_capacity(mol.bonds.len());
    loop {
        let white = |a: usize| !done[frag_of[a]];
        let start = match rooted_at_atom {
            Some(r) if white(r) => r,
            _ => {
                // The first unprocessed non-chiral atom, else the chiral
                // atom with the lowest rank (= index).
                let mut pick = None;
                let mut first_chiral = None;
                for a in 0..n {
                    if !white(a) {
                        continue;
                    }
                    if matches!(mol.atoms[a].chiral, ChiralTag::Cw | ChiralTag::Ccw) {
                        first_chiral.get_or_insert(a);
                    } else {
                        pick = Some(a);
                        break;
                    }
                }
                match pick.or(first_chiral) {
                    Some(a) => a,
                    None => break,
                }
            }
        };
        let f = frag_of[start];
        done[f] = true;
        let atoms = &frags[f];
        let mut sub = if frags.len() == 1 {
            mol.clone()
        } else {
            fragment(mol, atoms)
        };
        // Global bond index of each fragment bond (fragment keeps order).
        let sub_bonds: Vec<usize> = (0..mol.bonds.len())
            .filter(|&b| frag_of[mol.bonds[b].begin] == f)
            .collect();
        // FragmentSmartsConstruct: empty ring information and an updated
        // property cache; ranks are atom indices.
        sub.set_rings(Vec::new());
        sub.update_property_cache(false)?;
        let ranks: Vec<u32> = (0..sub.atoms.len() as u32).collect();
        let local_start = atoms.binary_search(&start).expect("start in fragment");
        let canon = canonicalize_fragment(&mut sub, local_start, &ranks, isomeric)?;
        if !res.is_empty() {
            res.push('.');
        }
        for e in &canon.stack {
            match *e {
                StackElem::Atom(a) => {
                    res.push_str(&atom_smarts(&sub, a, isomeric));
                    atom_order.push(atoms[a]);
                }
                StackElem::Bond(b, left) => {
                    res.push_str(&bond_smarts(&sub, b, left, isomeric, dative));
                    bond_order.push(sub_bonds[b]);
                }
                StackElem::Ring(num) => {
                    if num < 10 {
                        res.push_str(&num.to_string());
                    } else {
                        res.push('%');
                        res.push_str(&num.to_string());
                    }
                }
                StackElem::BranchOpen => res.push('('),
                StackElem::BranchClose => res.push(')'),
            }
        }
    }
    Ok((res, atom_order, bond_order))
}

/// `SmilesWrite::getCXExtensions(mol)` for the fields a `MolFromSmiles`
/// molecule without stereo groups can carry: radicals and coordinate
/// (dative) bonds.
pub(crate) fn cx_extensions(mol: &Mol, atom_order: &[usize], bond_order: &[usize]) -> String {
    let mut res = String::from("|");
    let mut rads: std::collections::BTreeMap<u32, Vec<usize>> = Default::default();
    for (i, &a) in atom_order.iter().enumerate() {
        let nrad = mol.atoms[a].radicals;
        if nrad != 0 {
            rads.entry(nrad).or_default().push(i);
        }
    }
    if !rads.is_empty() {
        let mut block = String::new();
        for (nrad, atoms) in &rads {
            match nrad {
                1 => block.push_str("^1:"),
                2 => block.push_str("^2:"),
                3 => block.push_str("^5:"),
                _ => {}
            }
            for a in atoms {
                block.push_str(&format!("{a},"));
            }
        }
        res.push_str(&block);
        if res.ends_with(',') {
            res.pop();
        }
    }
    let mut coord = String::new();
    for (i, &b) in bond_order.iter().enumerate() {
        let bond = &mol.bonds[b];
        if bond.bt != BondType::Dative {
            continue;
        }
        let beg = atom_order
            .iter()
            .position(|&a| a == bond.begin)
            .unwrap_or(atom_order.len());
        coord.push_str(if coord.is_empty() { "C:" } else { "," });
        coord.push_str(&format!("{beg}.{i}"));
    }
    if !coord.is_empty() {
        if res.len() > 1 {
            res.push(',');
        }
        res.push_str(&coord);
    }
    if res.len() > 1 {
        res.push('|');
        res
    } else {
        String::new()
    }
}
