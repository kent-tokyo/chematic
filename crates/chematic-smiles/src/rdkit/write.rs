//! `SmilesWrite` (RDKit 2026.03.1 `SmilesWrite.cpp`): `GetAtomSmiles`,
//! `GetBondSmiles`, `FragmentSmilesConstruct` and `MolToSmiles` with the
//! `SmilesWriteParams` of [`RdkitSmilesParams`] (`doRandom=false`).

use std::collections::BTreeMap;

use super::canon::{StackElem, canonicalize_fragment};
use super::kekulize::kekulize_ranked;
use super::mol::{BondDir, BondStereo, BondType, ChiralTag, Mol};
use super::periodic;
use super::rank::rank_mol_atoms_with;
use super::stereo::legacy_stereo_perception;
use super::{RdkitSmilesError, RdkitSmilesParams};

/// `SmilesWrite::inOrganicSubset`.
fn in_organic_subset(anum: u32) -> bool {
    matches!(anum, 0 | 5 | 6 | 7 | 8 | 9 | 15 | 16 | 17 | 35 | 53)
}

/// `atomNeedsBracket`.
fn atom_needs_bracket(mol: &Mol, a: usize, at_string: &str, isomeric: bool) -> bool {
    let atom = &mol.atoms[a];
    if !in_organic_subset(atom.anum) {
        return true;
    }
    if atom.charge != 0 {
        return true;
    }
    if isomeric && (atom.isotope != 0 || !at_string.is_empty()) {
        return true;
    }
    if atom.map.is_some() {
        return true;
    }
    let default_v = periodic::default_valence(atom.anum);
    let non_standard = if atom.radicals != 0 {
        true
    } else if (atom.anum == 7 || atom.anum == 15) && atom.aromatic && atom.num_explicit_hs != 0 {
        true
    } else {
        mol.total_valence(a) != default_v && mol.total_num_hs(a) != 0
    };
    if non_standard {
        return true;
    }
    mol.nbrs(a).any(|nb| periodic::is_metal(mol.atoms[nb].anum))
}

/// `SmilesWrite::GetAtomSmiles`.
fn atom_smiles(mol: &Mol, a: usize, p: &RdkitSmilesParams) -> String {
    let atom = &mol.atoms[a];
    let mut symb = periodic::symbol(atom.anum).to_string();
    let at_string = match atom.chiral {
        _ if !p.isomeric => "",
        ChiralTag::Cw => "@@",
        ChiralTag::Ccw => "@",
        ChiralTag::Unspecified => "",
    };
    let needs_bracket = p.all_hs_explicit || atom_needs_bracket(mol, a, at_string, p.isomeric);
    let mut res = String::new();
    if needs_bracket {
        res.push('[');
    }
    if atom.isotope != 0 && p.isomeric {
        res.push_str(&atom.isotope.to_string());
    }
    if !p.kekule
        && atom.aromatic
        && symb.as_bytes()[0].is_ascii_uppercase()
        && matches!(atom.anum, 5 | 6 | 7 | 8 | 14 | 15 | 16 | 33 | 34 | 52)
    {
        let lower = symb[..1].to_ascii_lowercase();
        symb.replace_range(..1, &lower);
    }
    res.push_str(&symb);
    res.push_str(at_string);
    if needs_bracket {
        let tot_hs = mol.total_num_hs(a);
        if tot_hs > 0 {
            res.push('H');
            if tot_hs > 1 {
                res.push_str(&tot_hs.to_string());
            }
        }
        let fc = atom.charge;
        if fc > 0 {
            res.push('+');
            if fc > 1 {
                res.push_str(&fc.to_string());
            }
        } else if fc < 0 {
            if fc < -1 {
                res.push_str(&fc.to_string());
            } else {
                res.push('-');
            }
        }
        if let Some(m) = atom.map {
            res.push(':');
            res.push_str(&m.to_string());
        }
        res.push(']');
    }
    res
}

/// `SmilesWrite::GetBondSmiles`.
fn bond_smiles(mol: &Mol, b: usize, atom_to_left: usize, p: &RdkitSmilesParams) -> &'static str {
    let bond = &mol.bonds[b];
    let mut aromatic = false;
    if !p.kekule
        && matches!(
            bond.bt,
            BondType::Single | BondType::Double | BondType::Aromatic
        )
    {
        let a1 = &mol.atoms[atom_to_left];
        let a2 = &mol.atoms[bond.other(atom_to_left)];
        if a1.aromatic && a2.aromatic && (a1.anum != 0 || a2.anum != 0) {
            aromatic = true;
        }
    }
    let write_dir = p.all_bonds_explicit || p.isomeric;
    let slash = |dir: BondDir| match dir {
        BondDir::EndDownRight => "\\",
        _ => "/",
    };
    match bond.bt {
        BondType::Single => match bond.dir {
            BondDir::None => {
                if p.all_bonds_explicit || (aromatic && !bond.aromatic) {
                    "-"
                } else {
                    ""
                }
            }
            dir if write_dir => slash(dir),
            _ => "",
        },
        BondType::Double => {
            if !aromatic || !bond.aromatic || p.all_bonds_explicit {
                "="
            } else {
                ""
            }
        }
        BondType::Triple => "#",
        BondType::Quadruple => "$",
        BondType::Aromatic => match bond.dir {
            BondDir::None => {
                if p.all_bonds_explicit || !aromatic {
                    ":"
                } else {
                    ""
                }
            }
            dir if write_dir => slash(dir),
            _ => "",
        },
        BondType::Dative => {
            if bond.begin == atom_to_left {
                "->"
            } else {
                "<-"
            }
        }
    }
}

/// `SmilesWrite::FragmentSmilesConstruct` for a whole fragment.
fn fragment_smiles_construct(
    mol: &mut Mol,
    start: usize,
    ranks: &[u32],
    p: &RdkitSmilesParams,
) -> Result<String, RdkitSmilesError> {
    if p.kekule {
        // `MolOps::Kekulize(mol)`: canonical, with `rankFragmentAtoms`
        // (chirality and isotopes included).
        let kek_ranks = rank_mol_atoms_with(mol, true);
        kekulize_ranked(mol, Some(&kek_ranks))?;
    }
    let canon = canonicalize_fragment(mol, start, ranks, p.isomeric)?;
    let mut res = String::new();
    let mut ring_closure_map: BTreeMap<u32, u32> = BTreeMap::new();
    let mut to_erase: Vec<u32> = Vec::new();
    for e in &canon.stack {
        match *e {
            StackElem::Atom(a) => {
                for r in to_erase.drain(..) {
                    ring_closure_map.remove(&r);
                }
                res.push_str(&atom_smiles(mol, a, p));
            }
            StackElem::Bond(b, left) => res.push_str(bond_smiles(mol, b, left, p)),
            StackElem::Ring(ring_idx) => {
                let closure_val = if let Some(&v) = ring_closure_map.get(&ring_idx) {
                    to_erase.push(ring_idx);
                    v
                } else {
                    let mut v = 1;
                    while ring_closure_map.values().any(|&x| x == v) {
                        v += 1;
                    }
                    ring_closure_map.insert(ring_idx, v);
                    v
                };
                if closure_val < 10 {
                    res.push(char::from(b'0' + closure_val as u8));
                } else if closure_val < 100 {
                    res.push('%');
                    res.push_str(&closure_val.to_string());
                } else {
                    res.push_str("%(");
                    res.push_str(&closure_val.to_string());
                    res.push(')');
                }
            }
            StackElem::BranchOpen => res.push('('),
            StackElem::BranchClose => res.push(')'),
        }
    }
    Ok(res)
}

/// `MolOps::getMolFrags`: atom lists of the connected components, numbered
/// from the lowest atom index, each list ascending.
fn mol_frags(mol: &Mol) -> Vec<Vec<usize>> {
    let n = mol.atoms.len();
    let mut frag_of = vec![usize::MAX; n];
    let mut frags: Vec<Vec<usize>> = Vec::new();
    for start in 0..n {
        if frag_of[start] != usize::MAX {
            continue;
        }
        let id = frags.len();
        frag_of[start] = id;
        let mut members = vec![start];
        let mut i = 0;
        while i < members.len() {
            let a = members[i];
            i += 1;
            for nb in mol.nbrs(a) {
                if frag_of[nb] == usize::MAX {
                    frag_of[nb] = id;
                    members.push(nb);
                }
            }
        }
        members.sort_unstable();
        frags.push(members);
    }
    frags
}

/// The fragment of `mol` on `atoms` (ascending), keeping atom and bond
/// order, properties and the rings inside it.
fn fragment(mol: &Mol, atoms: &[usize]) -> Mol {
    let mut new_idx = vec![usize::MAX; mol.atoms.len()];
    for (k, &a) in atoms.iter().enumerate() {
        new_idx[a] = k;
    }
    let mut out = Mol::default();
    for &a in atoms {
        let mut atom = mol.atoms[a].clone();
        if let Some(rsa) = &atom.ring_stereo_atoms {
            atom.ring_stereo_atoms = Some(
                rsa.iter()
                    .map(|&v| {
                        let idx = new_idx[(v.unsigned_abs() - 1) as usize] as i32 + 1;
                        if v < 0 { -idx } else { idx }
                    })
                    .collect(),
            );
        }
        out.add_atom(atom);
    }
    for bond in &mol.bonds {
        if new_idx[bond.begin] != usize::MAX && new_idx[bond.end] != usize::MAX {
            let mut b = bond.clone();
            b.begin = new_idx[bond.begin];
            b.end = new_idx[bond.end];
            b.stereo_atoms = b.stereo_atoms.iter().map(|&x| new_idx[x]).collect();
            b.requested = b
                .requested
                .map(|(sa, sb, trans)| (new_idx[sa], new_idx[sb], trans));
            out.add_bond(b);
        }
    }
    if let Some(ri) = &mol.rings {
        let rings: Vec<Vec<usize>> = ri
            .atom_rings
            .iter()
            .filter(|r| r.iter().all(|&a| new_idx[a] != usize::MAX))
            .map(|r| r.iter().map(|&a| new_idx[a]).collect())
            .collect();
        out.set_rings(rings);
    }
    out
}

/// `SmilesWrite::detail::MolToSmiles(mol, params)` on a molecule as
/// `MolFromSmiles` leaves it.
pub(crate) fn mol_to_smiles(mol: &Mol, p: &RdkitSmilesParams) -> Result<String, RdkitSmilesError> {
    if mol.atoms.is_empty() {
        return Ok(String::new());
    }
    if let Some(r) = p.rooted_at_atom
        && r >= mol.atoms.len()
    {
        return Err(RdkitSmilesError::Unsupported(
            "rootedAtAtom must be less than the number of atoms".into(),
        ));
    }
    let frags = mol_frags(mol);
    let n_frags = frags.len();
    let mut pieces: Vec<String> = Vec::with_capacity(n_frags);
    for atoms in &frags {
        // RDKit's fragment-local root: the root minus the fragment's first
        // atom index.
        let rooted = p
            .rooted_at_atom
            .filter(|r| atoms.binary_search(r).is_ok())
            .map(|r| r - atoms[0]);
        let mut tmol = if n_frags == 1 {
            mol.clone()
        } else {
            fragment(mol, atoms)
        };
        tmol.update_property_cache(false)?;
        if p.isomeric && n_frags > 1 {
            // The fragment copy lost `_StereochemDone`.
            legacy_stereo_perception(&mut tmol, true, false);
        }
        for b in &mut tmol.bonds {
            if b.stereo == BondStereo::Any {
                b.stereo = BondStereo::None;
            }
        }
        let ranks: Vec<u32> = if p.canonical {
            rank_mol_atoms_with(&tmol, p.isomeric)
        } else {
            (0..tmol.atoms.len() as u32).collect()
        };
        let start = match rooted {
            Some(r) if r < tmol.atoms.len() => r,
            Some(_) => {
                return Err(RdkitSmilesError::Unsupported(
                    "rootedAtAtom maps outside its fragment (RDKit indexes past the fragment)"
                        .into(),
                ));
            }
            None => (0..tmol.atoms.len())
                .min_by_key(|&i| ranks[i])
                .expect("non-empty fragment"),
        };
        pieces.push(fragment_smiles_construct(&mut tmol, start, &ranks, p)?);
    }
    if p.canonical {
        pieces.sort();
    }
    Ok(pieces.join("."))
}
