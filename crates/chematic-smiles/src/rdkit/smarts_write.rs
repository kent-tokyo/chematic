//! `SmartsWrite` (RDKit 2026.03.1 `SmartsWrite.cpp`) for molecules without
//! query features: `MolToSmarts` and `MolToCXSmarts` of a molecule read with
//! `MolFromSmiles`.
//!
//! RDKit writes such a molecule with atom indices as ranks (no
//! canonicalization) and ring information cleared, every atom as a bracket
//! `[#n]`/`[Sym]` atom and every bond explicitly.

use std::borrow::Cow;

use super::RdkitSmilesError;
use super::canon::{StackElem, canonicalize_fragment};
use super::mol::{BondDir, BondType, ChiralTag, Mol};
use super::periodic;
use super::write::{fragment, mol_frags};

/// `SmilesWrite::inOrganicSubset`.
fn in_organic_subset(anum: u32) -> bool {
    matches!(anum, 0 | 5 | 6 | 7 | 8 | 9 | 15 | 16 | 17 | 35 | 53)
}

/// `getNonQueryAtomSmarts`, appended to `res`.
fn atom_smarts(res: &mut String, mol: &Mol, a: usize, isomeric: bool) {
    use std::fmt::Write;
    let atom = &mol.atoms[a];
    res.push('[');
    if atom.isotope != 0 {
        let _ = write!(res, "{}", atom.isotope);
    }
    if in_organic_subset(atom.anum) {
        let _ = write!(res, "#{}", atom.anum);
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
            _ => {}
        }
    }
    if added_chirality && atom.num_explicit_hs == 1 {
        res.push('H');
    }
    match atom.charge {
        0 => {}
        -1 => res.push('-'),
        1 => res.push('+'),
        c if c < 0 => {
            let _ = write!(res, "{c}");
        }
        c => {
            let _ = write!(res, "+{c}");
        }
    }
    if let Some(m) = atom.map {
        let _ = write!(res, ":{m}");
    }
    res.push(']');
}

/// `getNonQueryBondSmarts`.
fn bond_smarts(
    mol: &Mol,
    b: usize,
    atom_to_left: usize,
    isomeric: bool,
    dative: bool,
) -> &'static str {
    let bond = &mol.bonds[b];
    let dir = |default: &'static str| -> &'static str {
        if isomeric {
            match bond.dir {
                BondDir::EndDownRight => return "\\",
                BondDir::EndUpRight => return "/",
                BondDir::None => {}
            }
        }
        default
    };
    if bond.aromatic {
        return dir(":");
    }
    match bond.bt {
        BondType::Single => dir("-"),
        BondType::Double => "=",
        BondType::Triple => "#",
        BondType::Quadruple => "$",
        BondType::Aromatic => dir(":"),
        BondType::Dative => {
            if !dative {
                "-"
            } else if bond.begin != atom_to_left {
                "<-"
            } else {
                "->"
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
    mol_to_smarts_cow(Cow::Borrowed(mol), isomeric, rooted_at_atom, dative)
}

/// [`mol_to_smarts`] on a molecule the caller no longer needs (a
/// single-fragment molecule is then written without a copy).
pub(crate) fn mol_to_smarts_owned(
    mol: Mol,
    isomeric: bool,
    rooted_at_atom: Option<usize>,
    dative: bool,
) -> Result<(String, Vec<usize>, Vec<usize>), RdkitSmilesError> {
    mol_to_smarts_cow(Cow::Owned(mol), isomeric, rooted_at_atom, dative)
}

/// The atom `MolToSmarts` starts the next fragment at among the atoms
/// `white` accepts: the root if allowed, else the first non-chiral atom,
/// else the chiral atom with the lowest rank (= index).
fn pick_start(
    mol: &Mol,
    rooted_at_atom: Option<usize>,
    white: impl Fn(usize) -> bool,
) -> Option<usize> {
    match rooted_at_atom {
        Some(r) if white(r) => Some(r),
        _ => {
            let mut first_chiral = None;
            for a in 0..mol.atoms.len() {
                if !white(a) {
                    continue;
                }
                if matches!(mol.atoms[a].chiral, ChiralTag::Cw | ChiralTag::Ccw) {
                    first_chiral.get_or_insert(a);
                } else {
                    return Some(a);
                }
            }
            first_chiral
        }
    }
}

/// Output of [`mol_to_smarts`] while it is being written.
struct SmartsOut {
    res: String,
    atom_order: Vec<usize>,
    bond_order: Vec<usize>,
}

impl SmartsOut {
    /// `FragmentSmartsConstruct` of fragment `sub` (atoms `atoms` and bonds
    /// `sub_bonds` of the whole molecule) from `local_start`.
    #[allow(clippy::too_many_arguments)]
    fn write_fragment(
        &mut self,
        mut sub: Mol,
        atoms: &[usize],
        sub_bonds: &[usize],
        local_start: usize,
        isomeric: bool,
        dative: bool,
    ) -> Result<(), RdkitSmilesError> {
        // Empty ring information and an updated property cache; ranks are
        // atom indices.
        sub.set_rings(Vec::new());
        sub.update_property_cache(false)?;
        let ranks: Vec<u32> = (0..sub.atoms.len() as u32).collect();
        let canon = canonicalize_fragment(&mut sub, local_start, &ranks, isomeric)?;
        let res = &mut self.res;
        if !res.is_empty() {
            res.push('.');
        }
        for e in &canon.stack {
            match *e {
                StackElem::Atom(a) => {
                    atom_smarts(res, &sub, a, isomeric);
                    self.atom_order.push(atoms[a]);
                }
                StackElem::Bond(b, left) => {
                    res.push_str(bond_smarts(&sub, b, left, isomeric, dative));
                    self.bond_order.push(sub_bonds[b]);
                }
                StackElem::Ring(num) => {
                    use std::fmt::Write;
                    if num >= 10 {
                        res.push('%');
                    }
                    let _ = write!(res, "{num}");
                }
                StackElem::BranchOpen => res.push('('),
                StackElem::BranchClose => res.push(')'),
            }
        }
        Ok(())
    }
}

fn mol_to_smarts_cow(
    mol: Cow<'_, Mol>,
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
    let mut out = SmartsOut {
        res: String::with_capacity(8 * n),
        atom_order: Vec::with_capacity(n),
        bond_order: Vec::with_capacity(mol.bonds.len()),
    };
    let frags = mol_frags(&mol);
    if frags.len() == 1 {
        let start = pick_start(&mol, rooted_at_atom, |_| true).expect("non-empty molecule");
        let sub_bonds: Vec<usize> = (0..mol.bonds.len()).collect();
        out.write_fragment(
            mol.into_owned(),
            &frags[0],
            &sub_bonds,
            start,
            isomeric,
            dative,
        )?;
        return Ok((out.res, out.atom_order, out.bond_order));
    }
    let mol: &Mol = &mol;
    let mut frag_of = vec![0usize; n];
    for (f, atoms) in frags.iter().enumerate() {
        for &a in atoms {
            frag_of[a] = f;
        }
    }
    let mut done = vec![false; frags.len()];
    while let Some(start) = pick_start(mol, rooted_at_atom, |a| !done[frag_of[a]]) {
        let f = frag_of[start];
        done[f] = true;
        let atoms = &frags[f];
        // Global bond index of each fragment bond (fragment keeps order).
        let sub_bonds: Vec<usize> = (0..mol.bonds.len())
            .filter(|&b| frag_of[mol.bonds[b].begin] == f)
            .collect();
        let local_start = atoms.binary_search(&start).expect("start in fragment");
        out.write_fragment(
            fragment(mol, atoms),
            atoms,
            &sub_bonds,
            local_start,
            isomeric,
            dative,
        )?;
    }
    Ok((out.res, out.atom_order, out.bond_order))
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
