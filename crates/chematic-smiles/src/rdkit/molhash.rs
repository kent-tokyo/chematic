//! RDKit's `MolHash::MolHash` (`Code/GraphMol/MolHash/hashfunctions.cpp`,
//! NextMove Software's molhash) on the RDKit-model molecule.

use std::collections::VecDeque;

use super::mol::{Atom, Bond, BondStereo, BondType, Mol};
use super::{RdkitSmilesError, RdkitSmilesParams, sanitize, smarts_match, stereo, write};

/// RDKit's `rdMolHash.HashFunction`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RdkitHashFunction {
    AnonymousGraph = 1,
    ElementGraph = 2,
    CanonicalSmiles = 3,
    MurckoScaffold = 4,
    ExtendedMurcko = 5,
    MolFormula = 6,
    AtomBondCounts = 7,
    DegreeVector = 8,
    Mesomer = 9,
    HetAtomTautomer = 10,
    HetAtomProtomer = 11,
    RedoxPair = 12,
    Regioisomer = 13,
    NetCharge = 14,
    SmallWorldIndexBR = 15,
    SmallWorldIndexBRL = 16,
    ArthorSubstructureOrder = 17,
    HetAtomTautomerv2 = 18,
    HetAtomProtomerv2 = 19,
}

impl RdkitHashFunction {
    /// Every hash function, in RDKit's enum order.
    pub const ALL: [RdkitHashFunction; 19] = [
        Self::AnonymousGraph,
        Self::ElementGraph,
        Self::CanonicalSmiles,
        Self::MurckoScaffold,
        Self::ExtendedMurcko,
        Self::MolFormula,
        Self::AtomBondCounts,
        Self::DegreeVector,
        Self::Mesomer,
        Self::HetAtomTautomer,
        Self::HetAtomProtomer,
        Self::RedoxPair,
        Self::Regioisomer,
        Self::NetCharge,
        Self::SmallWorldIndexBR,
        Self::SmallWorldIndexBRL,
        Self::ArthorSubstructureOrder,
        Self::HetAtomTautomerv2,
        Self::HetAtomProtomerv2,
    ];

    /// The enum member's name in `rdMolHash.HashFunction`.
    pub fn name(self) -> &'static str {
        match self {
            Self::AnonymousGraph => "AnonymousGraph",
            Self::ElementGraph => "ElementGraph",
            Self::CanonicalSmiles => "CanonicalSmiles",
            Self::MurckoScaffold => "MurckoScaffold",
            Self::ExtendedMurcko => "ExtendedMurcko",
            Self::MolFormula => "MolFormula",
            Self::AtomBondCounts => "AtomBondCounts",
            Self::DegreeVector => "DegreeVector",
            Self::Mesomer => "Mesomer",
            Self::HetAtomTautomer => "HetAtomTautomer",
            Self::HetAtomProtomer => "HetAtomProtomer",
            Self::RedoxPair => "RedoxPair",
            Self::Regioisomer => "Regioisomer",
            Self::NetCharge => "NetCharge",
            Self::SmallWorldIndexBR => "SmallWorldIndexBR",
            Self::SmallWorldIndexBRL => "SmallWorldIndexBRL",
            Self::ArthorSubstructureOrder => "ArthorSubstructureOrder",
            Self::HetAtomTautomerv2 => "HetAtomTautomerv2",
            Self::HetAtomProtomerv2 => "HetAtomProtomerv2",
        }
    }

    /// The function named `name` (case-insensitive).
    pub fn from_name(name: &str) -> Option<Self> {
        Self::ALL
            .into_iter()
            .find(|f| f.name().eq_ignore_ascii_case(name))
    }
}

/// `MolHash(mol, func, useCXSmiles)` on a molecule as `MolFromSmiles`
/// leaves it.
pub(crate) fn mol_hash(
    mut mol: Mol,
    func: RdkitHashFunction,
    use_cx: bool,
) -> Result<String, RdkitSmilesError> {
    use RdkitHashFunction as F;
    sanitize_hydrogens(&mut mol)?;
    // `convertToSmilesWithCXFlags` and, with `useCXSmiles`,
    // `addCXExtensions` (most functions skip the radicals there).
    let skip_radicals = !matches!(func, F::CanonicalSmiles | F::Regioisomer);
    let smiles = |mol: Mol, p: &RdkitSmilesParams| -> Result<String, RdkitSmilesError> {
        if use_cx {
            write::mol_to_cx_smiles_for_hash(mol, p, skip_radicals)
        } else {
            write::mol_to_smiles_owned(mol, p)
        }
    };
    let default = RdkitSmilesParams::default();
    Ok(match func {
        F::AnonymousGraph | F::ElementGraph => {
            anonymous_graph(&mut mol, func == F::ElementGraph)?;
            smiles(mol, &default)?
        }
        F::CanonicalSmiles => smiles(mol, &default)?,
        F::MurckoScaffold => {
            murcko_scaffold(&mut mol)?;
            smiles(mol, &default)?
        }
        F::ExtendedMurcko => {
            extended_murcko(&mut mol)?;
            smiles(mol, &default)?
        }
        F::Mesomer | F::RedoxPair => {
            let charge = mesomer(&mut mol)?;
            let mut res = smiles(mol, &default)?;
            if func == F::Mesomer {
                insert_suffix(&mut res, &format!("_{charge}"), use_cx);
            }
            res
        }
        F::HetAtomTautomer | F::HetAtomProtomer => {
            let (hcount, charge) = tautomer(&mut mol)?;
            let mut res = smiles(mol, &default)?;
            let suffix = if func == F::HetAtomTautomer {
                format!("_{hcount}_{charge}")
            } else {
                format!("_{}", hcount - charge)
            };
            insert_suffix(&mut res, &suffix, use_cx);
            res
        }
        F::HetAtomTautomerv2 | F::HetAtomProtomerv2 => {
            let (hcount, charge) = tautomer_v2(&mut mol)?;
            let p = RdkitSmilesParams {
                all_bonds_explicit: true,
                all_hs_explicit: true,
                ..default
            };
            let mut res = smiles(mol, &p)?;
            // `hcount` is unsigned in RDKit, printed with `%d`.
            let suffix = if func == F::HetAtomTautomerv2 {
                format!("_{}_{charge}", hcount as i32)
            } else {
                format!("_{}", (hcount as i32).wrapping_sub(charge))
            };
            insert_suffix(&mut res, &suffix, use_cx);
            res
        }
        F::Regioisomer => {
            regioisomer(&mut mol)?;
            smiles(mol, &default)?
        }
        F::MolFormula => molecular_formula(&mol),
        F::AtomBondCounts => format!("{},{}", mol.atoms.len(), mol.bonds.len()),
        F::NetCharge => mol.atoms.iter().map(|a| a.charge).sum::<i32>().to_string(),
        F::SmallWorldIndexBR | F::SmallWorldIndexBRL => {
            let acount = mol.atoms.len() as u32;
            let bcount = mol.bonds.len() as u32;
            let rcount = (bcount + 1).wrapping_sub(acount);
            if func == F::SmallWorldIndexBRL {
                let lcount = (0..mol.atoms.len()).filter(|&a| mol.degree(a) == 2).count();
                format!("B{bcount}R{rcount}L{lcount}")
            } else {
                format!("B{bcount}R{rcount}")
            }
        }
        F::DegreeVector => {
            let mut v = [0u32; 4];
            for a in 0..mol.atoms.len() {
                match mol.degree(a) {
                    4 => v[0] += 1,
                    3 => v[1] += 1,
                    2 => v[2] += 1,
                    1 => v[3] += 1,
                    _ => {}
                }
            }
            format!("{},{},{},{}", v[0], v[1], v[2], v[3])
        }
        F::ArthorSubstructureOrder => arthor_sub_order(&mol),
    })
}

/// The suffix goes after the SMILES, before a CXSMILES extension.
fn insert_suffix(res: &mut String, suffix: &str, use_cx: bool) {
    match res.find(' ').filter(|_| use_cx) {
        Some(p) => res.insert_str(p, suffix),
        None => res.push_str(suffix),
    }
}

/// `NMRDKitSanitizeHydrogens`: every atom's hydrogens become explicit.
fn sanitize_hydrogens(mol: &mut Mol) -> Result<(), RdkitSmilesError> {
    for a in 0..mol.atoms.len() {
        let hcount = mol.total_num_hs(a);
        mol.atoms[a].no_implicit = true;
        mol.atoms[a].num_explicit_hs = hcount;
        mol.update_atom_property_cache(a, false)?;
    }
    Ok(())
}

/// `MolOps::assignStereochemistry(mol, cleanIt=true, force=true)` after
/// the hash's edits (`legacyStereoPerception`, which finds rings first
/// when an edit reset them).
fn assign_stereochemistry(mol: &mut Mol) -> Result<(), RdkitSmilesError> {
    if (0..mol.atoms.len()).any(|a| mol.needs_update_property_cache(a)) {
        mol.update_property_cache(false)?;
    }
    if mol.rings.is_none() {
        mol.find_rings()?;
    }
    stereo::legacy_stereo_perception(mol, true, false);
    Ok(())
}

/// `NMRDKitBondGetOrder`.
fn bond_order(bt: BondType) -> u32 {
    match bt {
        BondType::Aromatic | BondType::Single => 1,
        BondType::Double => 2,
        BondType::Triple => 3,
        BondType::Quadruple => 4,
        _ => 0,
    }
}

/// `RWMol::removeBond`'s cleanup of the stereo atoms that refer across
/// bond `b` on the neighbouring bonds.
fn clear_stereo_atoms_across(mol: &mut Mol, b: usize) {
    let (u, v) = (mol.bonds[b].begin, mol.bonds[b].end);
    for (x, y) in [(u, v), (v, u)] {
        for k in 0..mol.atom_bonds[x].len() {
            let ob = mol.atom_bonds[x][k];
            if ob != b && mol.bonds[ob].stereo_atoms.contains(&y) {
                mol.bonds[ob].stereo_atoms.clear();
            }
        }
    }
}

/// `beginBatchEdit`, `removeAtom` on each of `removed` (ascending),
/// `commitBatchEdit`.
fn remove_atoms_batch(mol: &mut Mol, removed: &[usize]) {
    if removed.is_empty() {
        return;
    }
    let mut dead = vec![false; mol.bonds.len()];
    for &a in removed {
        for &b in &mol.atom_bonds[a] {
            dead[b] = true;
        }
    }
    for b in 0..mol.bonds.len() {
        if dead[b] {
            clear_stereo_atoms_across(mol, b);
        }
    }
    mol.remove_atoms(removed);
}

/// `AnonymousGraph(mol, elem)` up to the SMILES.
fn anonymous_graph(mol: &mut Mol, elem: bool) -> Result<(), RdkitSmilesError> {
    for a in 0..mol.atoms.len() {
        let atom = &mut mol.atoms[a];
        atom.aromatic = false;
        atom.charge = 0;
        if !elem {
            atom.num_explicit_hs = 0;
            atom.no_implicit = true;
            atom.anum = 0;
            atom.isotope = 0;
        } else {
            normalize_h_count(mol, a);
        }
    }
    for b in &mut mol.bonds {
        b.bt = BondType::Single;
        b.aromatic = false;
    }
    sanitize::assign_radicals(mol);
    assign_stereochemistry(mol)
}

/// `NormalizeHCount`.
fn normalize_h_count(mol: &mut Mol, a: usize) {
    let deg = mol.degree(a) as u32;
    let hcount = match mol.atoms[a].anum {
        9 | 17 | 35 | 53 => 1u32.saturating_sub(deg),
        8 | 16 => 2u32.saturating_sub(deg),
        5 => 3u32.saturating_sub(deg),
        7 | 15 => {
            if deg < 3 {
                3 - deg
            } else if deg == 4 {
                1
            } else {
                0
            }
        }
        6 => 4u32.saturating_sub(deg),
        _ => 0,
    };
    mol.atoms[a].no_implicit = true;
    mol.atoms[a].num_explicit_hs = hcount;
}

/// `MesomerHash` up to the SMILES; returns the net charge.
fn mesomer(mol: &mut Mol) -> Result<i32, RdkitSmilesError> {
    let mut charge = 0;
    for atom in &mut mol.atoms {
        charge += atom.charge;
        atom.aromatic = false;
        atom.charge = 0;
    }
    for b in &mut mol.bonds {
        b.bt = BondType::Single;
    }
    sanitize::assign_radicals(mol);
    assign_stereochemistry(mol)?;
    Ok(charge)
}

/// `TautomerHash` up to the SMILES; returns `(hcount, charge)`.
fn tautomer(mol: &mut Mol) -> Result<(i32, i32), RdkitSmilesError> {
    let mut hcount = 0i32;
    let mut charge = 0i32;
    for a in 0..mol.atoms.len() {
        charge += mol.atoms[a].charge;
        mol.atoms[a].aromatic = false;
        mol.atoms[a].charge = 0;
        if mol.atoms[a].anum != 6 {
            hcount += mol.total_num_hs(a) as i32;
            mol.atoms[a].no_implicit = true;
            mol.atoms[a].num_explicit_hs = 0;
        }
    }
    for b in 0..mol.bonds.len() {
        let bond = &mol.bonds[b];
        if bond.bt != BondType::Single
            && (bond.conjugated || mol.atoms[bond.begin].anum != 6 || mol.atoms[bond.end].anum != 6)
        {
            let bond = &mut mol.bonds[b];
            bond.aromatic = false;
            bond.bt = BondType::Single;
            bond.stereo = BondStereo::None;
        }
    }
    sanitize::assign_radicals(mol);
    assign_stereochemistry(mol)?;
    Ok((hcount, charge))
}

/// `details::getBondFlags`: bonds in carboxyl-like groups.
fn bond_flags(mol: &Mol) -> Vec<bool> {
    const PATTERNS: [&str; 12] = [
        "[C;!$(C-C(=[NH])-[NH2])]-[C;!$(C(-C)(=[NH])-[NH2])](=[O,N,S])-[O,N,S]",
        "[A;!$(C=[O,N,S])]-[O,N,S]-C=[O,N,S]",
        "[OH0,SH0]-C=[O,N,S]",
        "C-[N;$([N+]=C(N)(N)),$(N-C(N)=N)]",
        "[C]-[c](:[o,n,s]):[o,n,s]",
        "*[SD4](=*)=*",
        "*=[SD4;$(S(=*)=*)](*)*",
        "*-[H0]=,#[C,N]=,#*",
        "*[#7+]-[O-]",
        "[O,S;H]-P=O",
        "[#6]-P=O",
        "[#6]-N=[SD4]=*",
    ];
    thread_local! {
        static QUERIES: Vec<smarts_match::Query> =
            PATTERNS.iter().map(|p| smarts_match::parse(p)).collect();
    }
    let mut flags = vec![false; mol.bonds.len()];
    QUERIES.with(|queries| {
        for q in queries {
            for m in smarts_match::substruct_match(q, mol) {
                let b = mol.bond_between(m[0], m[1]).expect("matched bond");
                flags[b] = true;
            }
        }
    });
    flags
}

/// `queryAtomUnsaturated`.
fn unsaturated(mol: &Mol, a: usize) -> bool {
    (mol.total_degree(a) as i32) < mol.total_valence(a)
}

fn is_hetero(mol: &Mol, a: usize) -> bool {
    let n = mol.atoms[a].anum;
    n != 6 && n > 1
}

/// `isUnsaturatedBond`.
fn unsaturated_bond(bond: &Bond) -> bool {
    bond.aromatic
        || matches!(
            bond.bt,
            BondType::Aromatic | BondType::Double | BondType::Triple
        )
}

/// `TautomerHashv2` up to the SMILES; returns `(hcount, charge)`.
fn tautomer_v2(mol: &mut Mol) -> Result<(u32, i32), RdkitSmilesError> {
    let n_atoms = mol.atoms.len();
    let n_bonds = mol.bonds.len();
    let bond_flags = bond_flags(mol);
    // `atomFlags` are all zero.
    let candidate = |mol: &Mol, a: usize| mol.total_num_hs(a) != 0 || unsaturated(mol, a);
    let mut start_bonds = vec![false; n_bonds];
    for b in 0..n_bonds {
        if bond_flags[b] {
            continue;
        }
        let (beg, end) = (mol.bonds[b].begin, mol.bonds[b].end);
        let hetero_beg = is_hetero(mol, beg) && candidate(mol, beg);
        let hetero_end = is_hetero(mol, end) && candidate(mol, end);
        if !hetero_beg && !hetero_end {
            continue;
        }
        let unsat_beg = unsaturated(mol, beg);
        let unsat_end = unsaturated(mol, end);
        start_bonds[b] = (hetero_beg && unsat_end) || (hetero_end && unsat_beg);
    }
    let has_start_bond = |mol: &Mol, a: usize| mol.atom_bonds[a].iter().any(|&b| start_bonds[b]);
    let num_conj_nbrs: Vec<u32> = (0..n_atoms)
        .map(|a| {
            mol.nbrs(a)
                .filter(|&nb| {
                    mol.atom_bonds[nb]
                        .iter()
                        .any(|&nbb| mol.bonds[nbb].conjugated || start_bonds[nbb])
                })
                .count() as u32
        })
        .collect();
    let skip_neighbor_bond = |mol: &Mol, atom: usize, other: usize, nb: usize| {
        bond_flags[nb]
            || (!candidate(mol, other) && !has_start_bond(mol, other))
            || (!unsaturated_bond(&mol.bonds[nb])
                && !mol.bonds[nb].conjugated
                && !has_start_bond(mol, atom))
    };
    let check_for_overreach = |mol: &Mol, atom: usize, oatom: usize, bond: usize, nb: usize| {
        (start_bonds[bond] || has_start_bond(mol, atom))
            && !start_bonds[nb]
            && is_hetero(mol, atom)
            && !unsaturated_bond(&mol.bonds[nb])
            && num_conj_nbrs[oatom] < 2
    };

    let mut to_modify = vec![false; n_bonds];
    let mut considered = vec![false; n_bonds];
    for b in 0..n_bonds {
        if to_modify[b] || considered[b] || !start_bonds[b] {
            continue;
        }
        let mut conj_system = vec![false; n_bonds];
        let mut conj_atoms = vec![false; n_atoms];
        let mut atoms_in_system = vec![false; n_atoms];
        let mut num_donor_cs = 0u32;
        let mut active_hetero_hs = 0u32;
        let mut bq: VecDeque<usize> = VecDeque::new();
        for atm in [mol.bonds[b].begin, mol.bonds[b].end] {
            if mol.atoms[atm].anum == 6 {
                if mol.total_num_hs(atm) != 0 {
                    num_donor_cs += 1;
                }
            } else if is_hetero(mol, atm) {
                active_hetero_hs += mol.total_num_hs(atm);
            }
            for &nb in &mol.atom_bonds[atm] {
                if nb == b || considered[nb] {
                    continue;
                }
                let oatom = mol.bonds[nb].other(atm);
                if check_for_overreach(mol, atm, oatom, b, nb) {
                    continue;
                }
                if skip_neighbor_bond(mol, atm, oatom, nb)
                    && skip_neighbor_bond(mol, oatom, atm, nb)
                {
                    continue;
                }
                if !bq.contains(&nb) {
                    bq.push_back(nb);
                }
                considered[b] = true;
                conj_system[b] = true;
                conj_atoms[mol.bonds[b].begin] = true;
                conj_atoms[mol.bonds[b].end] = true;
            }
        }
        while let Some(bnd) = bq.pop_front() {
            if considered[bnd] {
                continue;
            }
            considered[bnd] = true;
            conj_system[bnd] = true;
            conj_atoms[mol.bonds[bnd].begin] = true;
            conj_atoms[mol.bonds[bnd].end] = true;
            for atm in [mol.bonds[bnd].begin, mol.bonds[bnd].end] {
                if atoms_in_system[atm] {
                    continue;
                }
                if mol.atoms[atm].anum == 6 {
                    if mol.total_num_hs(atm) != 0 {
                        num_donor_cs += 1;
                        atoms_in_system[atm] = true;
                    }
                } else if mol.atoms[atm].anum > 1 {
                    active_hetero_hs += mol.total_num_hs(atm);
                    atoms_in_system[atm] = true;
                }
                for &nb in &mol.atom_bonds[atm] {
                    if considered[nb] || bq.contains(&nb) {
                        continue;
                    }
                    let oatom = mol.bonds[nb].other(atm);
                    if check_for_overreach(mol, atm, oatom, bnd, nb) {
                        continue;
                    }
                    if skip_neighbor_bond(mol, atm, oatom, nb)
                        && skip_neighbor_bond(mol, oatom, atm, nb)
                    {
                        continue;
                    }
                    bq.push_back(nb);
                }
            }
        }
        let system_size = conj_system.iter().filter(|&&x| x).count();
        if system_size > 1 && (active_hetero_hs != 0 || num_donor_cs != 0) {
            for i in 0..n_atoms {
                if !conj_atoms[i] {
                    continue;
                }
                for &nb in &mol.atom_bonds[i] {
                    if conj_atoms[mol.bonds[nb].other(i)] {
                        to_modify[nb] = true;
                    }
                }
            }
        } else {
            for i in 0..n_bonds {
                if conj_system[i] {
                    considered[i] = false;
                }
            }
        }
    }

    for b in 0..n_bonds {
        if !start_bonds[b] {
            continue;
        }
        for atm in [mol.bonds[b].begin, mol.bonds[b].end] {
            for k in 0..mol.atom_bonds[atm].len() {
                let nb = mol.atom_bonds[atm][k];
                if nb == b || to_modify[nb] {
                    continue;
                }
                let oatom = mol.bonds[nb].other(atm);
                if mol.total_num_hs(oatom) == 0 {
                    continue;
                }
                let modified_nbr = mol.nbrs(oatom).filter(|&x| x != atm).any(|x| {
                    mol.atom_bonds[x]
                        .iter()
                        .any(|&nbnd| to_modify[nbnd] || start_bonds[nbnd])
                });
                if modified_nbr {
                    to_modify[nb] = true;
                }
            }
        }
    }

    let mut atoms_to_modify = vec![false; n_atoms];
    for b in 0..n_bonds {
        if !to_modify[b] {
            continue;
        }
        let bond = &mut mol.bonds[b];
        bond.aromatic = false;
        bond.bt = BondType::Aromatic;
        bond.stereo = BondStereo::None;
        atoms_to_modify[bond.begin] = true;
        atoms_to_modify[bond.end] = true;
    }
    let mut hcount = 0u32;
    let mut charge = 0i32;
    for a in 0..n_atoms {
        if !atoms_to_modify[a] {
            continue;
        }
        charge += mol.atoms[a].charge;
        hcount += mol.total_num_hs(a);
        let atom = &mut mol.atoms[a];
        atom.aromatic = false;
        atom.charge = 0;
        atom.no_implicit = true;
        atom.num_explicit_hs = 0;
    }
    // `!bondsToModify.empty() || !atomsToModify.empty()`: the bitsets'
    // sizes, not their contents.
    if n_bonds != 0 || n_atoms != 0 {
        sanitize::assign_radicals(mol);
        assign_stereochemistry(mol)?;
    }
    Ok((hcount, charge))
}

/// `IsInScaffold`: a ring atom, or one with two or more neighbours that
/// lead (without passing through it) to a ring atom.
fn is_in_scaffold(mol: &Mol, a: usize, in_ring: &[bool]) -> bool {
    if in_ring[a] {
        return true;
    }
    let n = mol.atoms.len();
    let mut count = 0;
    for nbor in mol.nbrs(a) {
        // `DepthFirstSearchForRing(a, nbor)`: a ring atom adjacent to an
        // atom reachable from `nbor` without passing through `a`.
        let mut visit = vec![false; n];
        visit[a] = true;
        visit[nbor] = true;
        let mut stack = vec![nbor];
        let mut found = false;
        'dfs: while let Some(x) = stack.pop() {
            for y in mol.nbrs(x) {
                if !visit[y] {
                    if in_ring[y] {
                        found = true;
                        break 'dfs;
                    }
                    visit[y] = true;
                    stack.push(y);
                }
            }
        }
        if found {
            count += 1;
        }
    }
    count > 1
}

/// `ExtendedMurckoScaffold` up to the SMILES.
fn extended_murcko(mol: &mut Mol) -> Result<(), RdkitSmilesError> {
    if mol.rings.is_none() {
        mol.find_rings()?;
    }
    let n = mol.atoms.len();
    let in_ring: Vec<bool> = (0..n).map(|a| mol.num_atom_rings(a) > 0).collect();
    let scaffold: Vec<bool> = (0..n).map(|a| is_in_scaffold(mol, a, &in_ring)).collect();
    let mut removed = Vec::new();
    for a in 0..n {
        if scaffold[a] {
            continue;
        }
        if mol.nbrs(a).any(|nb| scaffold[nb]) {
            let atom = &mut mol.atoms[a];
            atom.anum = 0;
            atom.charge = 0;
            atom.no_implicit = true;
            atom.num_explicit_hs = 0;
            atom.isotope = 0;
        } else {
            removed.push(a);
        }
    }
    remove_atoms_batch(mol, &removed);
    sanitize::assign_radicals(mol);
    assign_stereochemistry(mol)
}

/// `MurckoScaffoldHash` up to the SMILES: strip degree-0/1 atoms until none
/// is left.
fn murcko_scaffold(mol: &mut Mol) -> Result<(), RdkitSmilesError> {
    loop {
        let mut removed = Vec::new();
        for a in 0..mol.atoms.len() {
            let deg = mol.degree(a);
            if deg < 2 {
                if deg == 1 {
                    let b = mol.atom_bonds[a][0];
                    let nbr = mol.bonds[b].other(a);
                    let hcount = mol.total_num_hs(nbr);
                    mol.atoms[nbr].num_explicit_hs = hcount + bond_order(mol.bonds[b].bt);
                    mol.atoms[nbr].no_implicit = true;
                }
                removed.push(a);
            }
        }
        if removed.is_empty() {
            break;
        }
        remove_atoms_batch(mol, &removed);
    }
    sanitize::assign_radicals(mol);
    assign_stereochemistry(mol)
}

/// `RegioisomerBond`: -1 keeps the bond; otherwise bit 0 puts a `*` on the
/// begin atom, bit 1 one on the end atom (an H where the bit is clear).
fn regioisomer_bond(mol: &Mol, b: usize, in_ring: &[bool], bond_in_ring: &[bool]) -> i32 {
    let bond = &mol.bonds[b];
    if bond_order(bond.bt) != 1 || bond_in_ring[b] {
        return -1;
    }
    let (beg, end) = (bond.begin, bond.end);
    let (beg_elem, end_elem) = (mol.atoms[beg].anum, mol.atoms[end].anum);
    if beg_elem == 0 || end_elem == 0 {
        return -1;
    }
    let has_double = |a: usize| {
        mol.atom_bonds[a]
            .iter()
            .any(|&x| bond_order(mol.bonds[x].bt) == 2)
    };
    if in_ring[beg] {
        return if in_ring[end] { 0 } else { 2 };
    }
    if in_ring[end] {
        return 1;
    }
    if beg_elem != 6 && end_elem == 6 && !has_double(end) {
        return 1;
    }
    if beg_elem == 6 && end_elem != 6 && !has_double(beg) {
        return 2;
    }
    -1
}

/// `RegioisomerHash` up to the SMILES.
fn regioisomer(mol: &mut Mol) -> Result<(), RdkitSmilesError> {
    let molcpy = mol.clone();
    // RDKit runs `fastFindRings` on the copy (its test is inverted, so it
    // does whenever the molecule has rings perceived): unlike the SSSR it
    // also follows dative bonds.
    let (in_ring, bond_in_ring) = fast_ring_membership(&molcpy);
    for i in (0..molcpy.bonds.len()).rev() {
        let split = regioisomer_bond(&molcpy, i, &in_ring, &bond_in_ring);
        if split < 0 {
            continue;
        }
        let (beg, end) = (mol.bonds[i].begin, mol.bonds[i].end);
        clear_stereo_atoms_across(mol, i);
        mol.remove_bond(i);
        for a in [beg, end] {
            for k in 0..mol.atom_bonds[a].len() {
                let b = mol.atom_bonds[a][k];
                if matches!(mol.bonds[b].stereo, BondStereo::Z | BondStereo::E) {
                    mol.bonds[b].stereo = BondStereo::Any;
                }
            }
        }
        for (a, bit) in [(beg, 1), (end, 2)] {
            if split & bit != 0 {
                let mut star = Atom::new(0);
                star.no_implicit = true;
                let s = mol.add_atom(star);
                mol.add_bond(Bond::new(a, s, BondType::Single));
            } else {
                let hcount = mol.total_num_hs(a);
                mol.atoms[a].num_explicit_hs = hcount + 1;
                mol.atoms[a].no_implicit = true;
            }
        }
    }
    assign_stereochemistry(mol)
}

/// Ring membership as `fastFindRings` perceives it (every bond type): the
/// bonds on a cycle (not bridges) and their atoms.
fn fast_ring_membership(mol: &Mol) -> (Vec<bool>, Vec<bool>) {
    let n = mol.atoms.len();
    let mut disc = vec![usize::MAX; n];
    let mut low = vec![0usize; n];
    let mut bond_in_ring = vec![true; mol.bonds.len()];
    let mut time = 0;
    for root in 0..n {
        if disc[root] != usize::MAX {
            continue;
        }
        // (atom, bond it was entered by, next position in its bond list)
        let mut stack: Vec<(usize, usize, usize)> = vec![(root, usize::MAX, 0)];
        disc[root] = time;
        low[root] = time;
        time += 1;
        while let Some(top) = stack.last_mut() {
            let (a, via) = (top.0, top.1);
            if top.2 < mol.atom_bonds[a].len() {
                let b = mol.atom_bonds[a][top.2];
                top.2 += 1;
                if b == via {
                    continue;
                }
                let o = mol.bonds[b].other(a);
                if disc[o] == usize::MAX {
                    disc[o] = time;
                    low[o] = time;
                    time += 1;
                    stack.push((o, b, 0));
                } else {
                    low[a] = low[a].min(disc[o]);
                }
            } else {
                stack.pop();
                if let Some(&(p, _, _)) = stack.last() {
                    low[p] = low[p].min(low[a]);
                    if low[a] > disc[p] {
                        bond_in_ring[via] = false;
                    }
                }
            }
        }
    }
    let mut atom_in_ring = vec![false; n];
    for (b, bond) in mol.bonds.iter().enumerate() {
        if bond_in_ring[b] {
            atom_in_ring[bond.begin] = true;
            atom_in_ring[bond.end] = true;
        }
    }
    (atom_in_ring, bond_in_ring)
}

/// `NMMolecularFormula(mol)`.
fn molecular_formula(mol: &Mol) -> String {
    let mut hist = [0u32; 256];
    let mut charge = 0i32;
    for a in 0..mol.atoms.len() {
        let elem = mol.atoms[a].anum as usize;
        if elem < 256 {
            hist[elem] += 1;
        } else {
            hist[0] += 1;
        }
        hist[1] += mol.total_num_hs(a);
        charge += mol.atoms[a].charge;
    }
    // `OrganicHillOrder` (C, H, then by symbol) or `InorganicHillOrder` (by
    // symbol); element 0 is "X".
    let symbol = |e: usize| -> &'static str {
        if e == 0 {
            "X"
        } else {
            super::periodic::symbol(e as u32)
        }
    };
    let mut order: Vec<usize> = (0..119).collect();
    if hist[6] != 0 {
        order.retain(|&e| e != 6 && e != 1);
        order.sort_by_key(|&e| symbol(e));
        order.splice(0..0, [6, 1]);
    } else {
        order.sort_by_key(|&e| symbol(e));
    }
    let mut res = String::new();
    for e in order {
        if hist[e] > 0 {
            res.push_str(symbol(e));
            if hist[e] > 1 {
                res.push_str(&hist[e].to_string());
            }
        }
    }
    if charge > 0 {
        res.push('+');
        if charge > 1 {
            res.push_str(&charge.to_string());
        }
    } else if charge < 0 {
        res.push('-');
        if charge < -1 {
            res.push_str(&(-charge).to_string());
        }
    }
    res
}

/// `ArthorSubOrderHash`.
fn arthor_sub_order(mol: &Mol) -> String {
    let acount = mol.atoms.len() as u32;
    let bcount = mol.bonds.len() as u32;
    let pcount = write::mol_frags(mol).len() as u32;
    let (mut ccount, mut ocount, mut zcount, mut icount, mut qcount, mut rcount) =
        (0u32, 0u32, 0u32, 0u32, 0u32, 0u32);
    for a in 0..mol.atoms.len() {
        let atom = &mol.atoms[a];
        let elem = atom.anum;
        let charge = atom.charge;
        let valence = mol.total_valence(a);
        match elem {
            6 => {
                ccount += 1;
                if charge == 0 && valence != 4 {
                    rcount += 1;
                }
            }
            7 | 15 => {
                ocount += 1;
                if charge == 0 && valence != 3 && valence != 5 {
                    rcount += 1;
                }
            }
            8 => {
                ocount += 1;
                if charge != 0 && valence != 2 {
                    rcount += 1;
                }
            }
            9 => {
                ocount += 1;
                if charge != 0 && valence != 1 {
                    rcount += 1;
                }
            }
            17 | 35 | 53 => {
                ocount += 1;
                if charge == 0 && !matches!(valence, 1 | 3 | 5 | 7) {
                    rcount += 1;
                }
            }
            16 => {
                ocount += 1;
                if charge == 0 && !matches!(valence, 2 | 4 | 6) {
                    rcount += 1;
                }
            }
            _ => {}
        }
        zcount += elem;
        if atom.isotope != 0 {
            icount += 1;
        }
        if charge != 0 {
            qcount += 1;
        }
    }
    format!(
        "{:04x}{:04x}{:02x}{:04x}{:04x}{:06x}{:02x}{:02x}{:02x}",
        acount.min(0xffff),
        bcount.min(0xffff),
        pcount.min(0xff),
        ccount.min(0xffff),
        ocount.min(0xffff),
        zcount.min(0xffffff),
        rcount.min(0xff),
        qcount.min(0xff),
        icount.min(0xff)
    )
}
