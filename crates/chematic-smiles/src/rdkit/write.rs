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

/// `SmilesWrite::GetAtomSmiles`, appended to `res`.
fn atom_smiles(res: &mut String, mol: &Mol, a: usize, p: &RdkitSmilesParams) {
    use std::fmt::Write;
    let atom = &mol.atoms[a];
    let symb = periodic::symbol(atom.anum);
    let nontet_string;
    let at_string = match atom.chiral {
        _ if !p.isomeric => "",
        ChiralTag::Cw => "@@",
        ChiralTag::Ccw => "@",
        ChiralTag::Unspecified => "",
        tag => {
            // `getAtomChiralityInfo`: the class, then the permutation when set.
            let class = tag.nontet().expect("non-tetrahedral");
            nontet_string = if atom.chiral_perm == 0 {
                format!("@{}", class.token())
            } else {
                format!("@{}{}", class.token(), atom.chiral_perm)
            };
            nontet_string.as_str()
        }
    };
    let needs_bracket = p.all_hs_explicit || atom_needs_bracket(mol, a, at_string, p.isomeric);
    if needs_bracket {
        res.push('[');
    }
    if atom.isotope != 0 && p.isomeric {
        let _ = write!(res, "{}", atom.isotope);
    }
    if !p.kekule
        && atom.aromatic
        && symb.as_bytes()[0].is_ascii_uppercase()
        && matches!(atom.anum, 5 | 6 | 7 | 8 | 14 | 15 | 16 | 33 | 34 | 52)
    {
        res.push(char::from(symb.as_bytes()[0].to_ascii_lowercase()));
        res.push_str(&symb[1..]);
    } else {
        res.push_str(symb);
    }
    res.push_str(at_string);
    if needs_bracket {
        let tot_hs = mol.total_num_hs(a);
        if tot_hs > 0 {
            res.push('H');
            if tot_hs > 1 {
                let _ = write!(res, "{tot_hs}");
            }
        }
        let fc = atom.charge;
        if fc > 0 {
            res.push('+');
            if fc > 1 {
                let _ = write!(res, "{fc}");
            }
        } else if fc < 0 {
            if fc < -1 {
                let _ = write!(res, "{fc}");
            } else {
                res.push('-');
            }
        }
        if let Some(m) = atom.map {
            let _ = write!(res, ":{m}");
        }
        res.push(']');
    }
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
            BondDir::None if p.all_bonds_explicit || (aromatic && !bond.aromatic) => "-",
            BondDir::None => "",
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
            BondDir::None if p.all_bonds_explicit || !aromatic => ":",
            BondDir::None => "",
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
) -> Result<Piece, RdkitSmilesError> {
    if p.kekule {
        // `MolOps::Kekulize(mol)`: canonical, with `rankFragmentAtoms`
        // (chirality and isotopes included).
        let kek_ranks = rank_mol_atoms_with(mol, true);
        kekulize_ranked(mol, Some(&kek_ranks))?;
    }
    let canon = canonicalize_fragment(mol, start, ranks, p.isomeric)?;
    let mut res = String::with_capacity(2 * canon.stack.len());
    let mut ring_closure_map: BTreeMap<u32, u32> = BTreeMap::new();
    let mut to_erase: Vec<u32> = Vec::new();
    let mut atom_order = Vec::new();
    let mut bond_order = Vec::new();
    for e in &canon.stack {
        match *e {
            StackElem::Atom(a) => {
                for r in to_erase.drain(..) {
                    ring_closure_map.remove(&r);
                }
                atom_smiles(&mut res, mol, a, p);
                atom_order.push(a);
            }
            StackElem::Bond(b, left) => {
                res.push_str(bond_smiles(mol, b, left, p));
                bond_order.push(b);
            }
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
    Ok((res, atom_order, bond_order))
}

/// A fragment's SMILES with its atom and bond output orders
/// (`_smilesAtomOutputOrder`, `_smilesBondOutputOrder`).
type Piece = (String, Vec<usize>, Vec<usize>);

/// `MolOps::getMolFrags`: atom lists of the connected components, numbered
/// from the lowest atom index, each list ascending.
pub(crate) fn mol_frags(mol: &Mol) -> Vec<Vec<usize>> {
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
pub(crate) fn fragment(mol: &Mol, atoms: &[usize]) -> Mol {
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
    mol_to_smiles_cow(std::borrow::Cow::Borrowed(mol), p)
}

/// [`mol_to_smiles`] on a molecule the caller no longer needs (a
/// single-fragment molecule is then written without a copy).
pub(crate) fn mol_to_smiles_owned(
    mol: Mol,
    p: &RdkitSmilesParams,
) -> Result<String, RdkitSmilesError> {
    mol_to_smiles_cow(std::borrow::Cow::Owned(mol), p)
}

fn mol_to_smiles_cow(
    mol: std::borrow::Cow<'_, Mol>,
    p: &RdkitSmilesParams,
) -> Result<String, RdkitSmilesError> {
    Ok(mol_to_smiles_ordered(mol, p, false, false)?.0)
}

/// [`mol_to_smiles`] on a molecule whose `_StereochemDone` is a
/// non-computed property (set by the tautomer enumerator), which edited
/// fragment copies keep.
pub(crate) fn mol_to_smiles_flag_kept(
    mol: &Mol,
    p: &RdkitSmilesParams,
) -> Result<String, RdkitSmilesError> {
    Ok(mol_to_smiles_ordered(std::borrow::Cow::Borrowed(mol), p, false, true)?.0)
}

/// `fragmentHasChallengingFeatures` (no substance or stereo groups here).
fn challenging_fragment(mol: &Mol, atoms: &[usize]) -> bool {
    atoms.iter().any(|&a| {
        !matches!(mol.atoms[a].chiral, ChiralTag::Unspecified)
            || mol.atom_bonds[a].iter().any(|&b| {
                let bond = &mol.bonds[b];
                atoms.binary_search(&bond.other(a)).is_ok()
                    && !matches!(bond.stereo, BondStereo::None | BondStereo::Any)
            })
    })
}

/// `SmilesWrite::detail::MolToSmiles(mol, p, doingCXSmiles)` with the
/// output orders it stores on the molecule (indices into `mol`).
fn mol_to_smiles_ordered(
    mol: std::borrow::Cow<'_, Mol>,
    p: &RdkitSmilesParams,
    cx: bool,
    flag_kept: bool,
) -> Result<Piece, RdkitSmilesError> {
    if mol.atoms.is_empty() {
        return Ok(Default::default());
    }
    if let Some(r) = p.rooted_at_atom
        && r >= mol.atoms.len()
    {
        return Err(RdkitSmilesError::Unsupported(
            "rootedAtAtom must be less than the number of atoms".into(),
        ));
    }
    let frags = mol_frags(&mol);
    let n_frags = frags.len();
    if n_frags == 1 {
        // RDKit's fragment-local root: the root minus the fragment's first
        // atom index (0 here).
        let rooted = p.rooted_at_atom;
        return fragment_piece(mol.into_owned(), rooted, false, p, cx);
    }
    let mut pieces: Vec<Piece> = Vec::with_capacity(n_frags);
    for atoms in &frags {
        let rooted = p
            .rooted_at_atom
            .filter(|r| atoms.binary_search(r).is_ok())
            .map(|r| r - atoms[0]);
        // `getMolFrags` copies a fragment with `copyMolSubset` (no molecule
        // properties) when it is a single atom or one of more than three
        // simple fragments, and as an edited copy of the molecule otherwise
        // (keeping a non-computed `_StereochemDone`).
        let subset_copy = atoms.len() == 1 || (n_frags > 3 && !challenging_fragment(&mol, atoms));
        let copied = !flag_kept || subset_copy;
        let (smi, mut atom_order, mut bond_order) =
            fragment_piece(fragment(&mol, atoms), rooted, copied, p, cx)?;
        // `fragment` keeps the bonds inside the fragment in their order.
        let frag_bonds: Vec<usize> = (0..mol.bonds.len())
            .filter(|&b| atoms.binary_search(&mol.bonds[b].begin).is_ok())
            .collect();
        for a in &mut atom_order {
            *a = atoms[*a];
        }
        for b in &mut bond_order {
            *b = frag_bonds[*b];
        }
        pieces.push((smi, atom_order, bond_order));
    }
    if p.canonical {
        pieces.sort();
    }
    let mut res = (String::new(), Vec::new(), Vec::new());
    for (i, (smi, atoms, bonds)) in pieces.into_iter().enumerate() {
        if i > 0 {
            res.0.push('.');
        }
        res.0.push_str(&smi);
        res.1.extend(atoms);
        res.2.extend(bonds);
    }
    Ok(res)
}

/// One fragment's SMILES; `copied`: the fragment is a copy that lost
/// `_StereochemDone`.
fn fragment_piece(
    mut tmol: Mol,
    rooted: Option<usize>,
    copied: bool,
    p: &RdkitSmilesParams,
    cx: bool,
) -> Result<Piece, RdkitSmilesError> {
    tmol.update_property_cache(false)?;
    if p.isomeric && copied {
        legacy_stereo_perception(&mut tmol, true, false);
    }
    if cx {
        // Coordinate bonds go to the extension; the begin atom's explicit
        // valence is recomputed.
        for b in 0..tmol.bonds.len() {
            if tmol.bonds[b].bt == BondType::Dative {
                tmol.bonds[b].bt = BondType::Single;
                let begin = tmol.bonds[b].begin;
                tmol.calc_explicit_valence(begin, false)?;
            }
        }
    } else {
        for b in &mut tmol.bonds {
            if b.stereo == BondStereo::Any {
                b.stereo = BondStereo::None;
            }
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
                "rootedAtAtom maps outside its fragment (RDKit indexes past the fragment)".into(),
            ));
        }
        None => (0..tmol.atoms.len())
            .min_by_key(|&i| ranks[i])
            .expect("non-empty fragment"),
    };
    fragment_smiles_construct(&mut tmol, start, &ranks, p)
}

/// `SmilesWrite::detail::MolToSmiles(mol, p, doingCXSmiles=true)` followed
/// by MolHash's `addCXExtensions` (`skip_radicals`: without the radical
/// field).
pub(crate) fn mol_to_cx_smiles_for_hash(
    mol: Mol,
    p: &RdkitSmilesParams,
    skip_radicals: bool,
) -> Result<String, RdkitSmilesError> {
    if mol.atoms.is_empty() {
        // `MolToSmiles` stores no output order for an empty molecule, and
        // `getCXExtensions` raises on the missing property.
        return Err(RdkitSmilesError::Unsupported(
            "CXSMILES of an empty molecule (RDKit raises KeyError '_smilesAtomOutputOrder')".into(),
        ));
    }
    let (mut res, atoms, bonds) =
        mol_to_smiles_ordered(std::borrow::Cow::Borrowed(&mol), p, true, false)?;
    let ext = hash_cx_extensions(&mol, &atoms, &bonds, skip_radicals);
    if !ext.is_empty() {
        res.push(' ');
        res.push_str(&ext);
    }
    Ok(res)
}

/// `getCXExtensions(mol, CX_ALL ^ skipped)` for a hash molecule read from
/// SMILES without stereo groups: radicals, ring double bonds of unknown
/// configuration (`ctu`) and coordinate bonds.
fn hash_cx_extensions(mol: &Mol, atoms: &[usize], bonds: &[usize], skip_radicals: bool) -> String {
    let mut res = String::from("|");
    if !skip_radicals {
        let mut rads: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
        for (i, &a) in atoms.iter().enumerate() {
            let nrad = mol.atoms[a].radicals;
            if nrad != 0 {
                rads.entry(nrad).or_default().push(i);
            }
        }
        for (nrad, idxs) in &rads {
            res.push_str(match nrad {
                1 => "^1:",
                2 => "^2:",
                3 => "^5:",
                _ => continue,
            });
            for i in idxs {
                res.push_str(&format!("{i},"));
            }
        }
        if res.ends_with(',') {
            res.pop();
        }
    }
    let append = |res: &mut String, block: &str| {
        if block.is_empty() {
            return;
        }
        if res.len() > 1 {
            res.push(',');
        }
        res.push_str(block);
    };
    // `get_ringbond_cistrans_block`: STEREOANY ring double bonds in rings of
    // at least `minRingSizeForDoubleBondStereo` (8).
    if let Some(ri) = &mol.rings {
        let mut ctu = String::new();
        for (i, &b) in bonds.iter().enumerate() {
            if ri.num_bond_rings(b) == 0 || ri.min_bond_ring_size(b) < 8 {
                continue;
            }
            let bond = &mol.bonds[b];
            if !matches!(bond.bt, BondType::Double | BondType::Aromatic)
                || bond.stereo != BondStereo::Any
            {
                continue;
            }
            ctu.push_str(if ctu.is_empty() { "ctu:" } else { "," });
            ctu.push_str(&i.to_string());
        }
        append(&mut res, &ctu);
    }
    let mut coord = String::new();
    for (i, &b) in bonds.iter().enumerate() {
        let bond = &mol.bonds[b];
        if bond.bt != BondType::Dative {
            continue;
        }
        let beg = atoms
            .iter()
            .position(|&a| a == bond.begin)
            .unwrap_or(atoms.len());
        coord.push_str(if coord.is_empty() { "C:" } else { "," });
        coord.push_str(&format!("{beg}.{i}"));
    }
    append(&mut res, &coord);
    if res.len() > 1 {
        res.push('|');
        res
    } else {
        String::new()
    }
}
