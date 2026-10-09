//! `MolStandardize::TautomerEnumerator` (RDKit 2026.03.1 `Tautomer.cpp`,
//! `TautomerCatalog/tautomerTransforms.in`) with the default
//! `CleanupParameters`: `enumerate`, `canonicalize` (`pickCanonical` with
//! `TautomerScoringFunctions::scoreTautomer`) and the scoring functions.
//!
//! The enumeration walks RDKit's `std::map<std::string, Tautomer>` (keyed by
//! canonical SMILES) the way RDKit does, including entries inserted ahead
//! of or behind the iterator, so the transforms run in the same order and
//! the `maxTransforms`/`maxTautomers` limits cut at the same place.

use std::collections::BTreeMap;
use std::ops::Bound::{Excluded, Unbounded};
use std::sync::OnceLock;

use super::kekulize::{CANT_KEKULIZE, kekulize_full};
use super::mol::{BondStereo, BondType, ChiralTag, Mol};
use super::rank::rank_mol_atoms;
use super::sanitize::{adjust_hs, set_conjugation, set_hybridization};
use super::smarts_match::{self, Query};
use super::stereo::legacy_stereo_perception;
use super::{RdkitSmilesError, RdkitSmilesParams, aromaticity, write};

/// `defaults::defaultTautomerTransforms`: name, SMARTS, bond types,
/// charges.
const TRANSFORMS: [(&str, &str, &str, &str); 37] = [
    (
        "1,3 (thio)keto/enol f",
        "[CX4!H0R{0-2}]-[C;z{1-2}]=[O,S,Se,Te;X1]",
        "",
        "",
    ),
    (
        "1,3 (thio)keto/enol r",
        "[O,S,Se,Te;X2!H0]-[#6;z{1-2}]=[C,cz{0-1}R{0-1}]",
        "",
        "",
    ),
    (
        "1,5 (thio)keto/enol f",
        "[CX4z0,NX3;!H0]-[C]=[C][Cz1H0]=[O,S,Se,Te;X1]",
        "",
        "",
    ),
    (
        "1,5 (thio)keto/enol r",
        "[O,S,Se,Te;X2!H0]-[Cz1H0]=[C]-[C]=[Cz0,N]",
        "",
        "",
    ),
    (
        "aliphatic imine f",
        "[CX4R{0-2}!H0]-[Cz1]=[NX2;N!H0,$(N-a)]",
        "",
        "",
    ),
    (
        "aliphatic imine r",
        "[NX3;NH2,$([NH1]-a)]-[C;z{1-2}]=[CX3]",
        "",
        "",
    ),
    ("special imine f", "[Nz0!H0]-[C]=[Cz0X3R0]", "", ""),
    ("special imine r", "[Cz0R0X4!H0]-[c]=[nz0]", "", ""),
    (
        "special imine r2",
        "[Cz0R0X4!H0]-[cr6](=[cr6])-[nr6H0z0]",
        "==-",
        "",
    ),
    (
        "1,3 aromatic heteroatom H shift f",
        "[#7+0!H0]-[#6R1]=[O,#7X2+0]",
        "",
        "",
    ),
    (
        "1,3 aromatic heteroatom H shift ",
        "[O,#7+0;!H0]-[#6R1]=[#7+0X2]",
        "",
        "",
    ),
    (
        "1,3 heteroatom H shift",
        "[#7+0,S,O,Se,Te;!H0]-[#7X2,#6,#15X3H0]=[#7+0,#16,#8,Se,Te]",
        "",
        "",
    ),
    (
        "1,5 aromatic heteroatom H shift",
        "[#7+0,#16,#8;!H0]-[#6,#7]=[#6]-[#6,#7]=[#7+0,#16,#8;H0]",
        "",
        "",
    ),
    (
        "1,5 aromatic heteroatom H shift ",
        "[#7+0,#16,#8,Se,Te;!H0]-[#6,nX2]=[#6,nX2]-[#6,#7X2]=[#7X2+0,S,O,Se,Te]",
        "",
        "",
    ),
    (
        "1,5 aromatic heteroatom H shift r",
        "[#7+0,S,O,Se,Te;!H0]-[#6,#7X2]=[#6,nX2]-[#6,nX2]=[#7+0,#16,#8,Se,Te]",
        "",
        "",
    ),
    (
        "1,7 aromatic heteroatom H shift f",
        "[#7+0,#8,#16,Se,Te;!H0]-[#6,#7X2]=[#6,#7X2]-[#6,#7X2]=[#6]-[#6,#7X2]=[#7X2+0,S,O,Se,Te,Cz0X3]",
        "",
        "",
    ),
    (
        "1,7 aromatic heteroatom H shift r",
        "[#7+0,S,O,Se,Te,Cz0X4;!H0]-[#6,#7X2]=[#6]-[#6,#7X2]=[#6,#7X2]-[#6,#7X2]=[NX2,S,O,Se,Te]",
        "",
        "",
    ),
    (
        "1,9 aromatic heteroatom H shift f",
        "[#7+0,O;!H0]-[#6,#7X2]=[#6,#7X2]-[#6,#7X2]=[#6,#7X2]-[#6,#7X2]=[#6,#7X2]-[#6,#7X2]=[#7+0,O]",
        "",
        "",
    ),
    (
        "1,11 aromatic heteroatom H shift f",
        "[#7+0,O;!H0]-[#6,nX2]=[#6,nX2]-[#6,nX2]=[#6,nX2]-[#6,nX2]=[#6,nX2]-[#6,nX2]=[#6,nX2]-[#6,nX2]=[#7X2+0,O]",
        "",
        "",
    ),
    (
        "furanone f",
        "[O,S,N;!H0]-[#6z2r5]=,:[#6X3r5;$([#6]([#6;r5])=,:[#6X3r5])]",
        "",
        "",
    ),
    (
        "furanone r",
        "[#6r5!H0;$([#6]([#6r5])[#6r5])][#6z2r5]=[O,S,N]",
        "",
        "",
    ),
    ("keten/ynol f", "[C!H0]=[C]=[O,S,Se,Te;X1]", "#-", ""),
    ("keten/ynol r", "[O,S,Se,Te;!H0X2]-[C]#[C]", "==", ""),
    (
        "ionic nitro/aci-nitro f",
        "[C!H0]-[N+;$([N][O-])]=[O]",
        "",
        "",
    ),
    (
        "ionic nitro/aci-nitro r",
        "[O!H0]-[N+;$([N][O-])]=[C]",
        "",
        "",
    ),
    ("oxim/nitroso f", "[O!H0]-[Nz1]=[C]", "", ""),
    ("oxim/nitroso r", "[C!H0]-[Nz1]=[O]", "", ""),
    (
        "oxim/nitroso via phenol f",
        "[O!H0]-[N]=[C]-[C]=[C]-[C]=[OH0]",
        "",
        "",
    ),
    (
        "oxim/nitroso via phenol r",
        "[O!H0]-[c]=,:[c][c]=,:[c]-[N]=[OH0]",
        "",
        "",
    ),
    ("cyano/iso-cyanic acid f", "[O!H0]-[C]#[N]", "==", ""),
    ("cyano/iso-cyanic acid r", "[N!H0]=[C]=[O]", "#-", ""),
    (
        "formamidinesulfinic acid f",
        "[O,N;!H0]-[C]=[S,Se,Te;v6]=[O]",
        "=--",
        "",
    ),
    (
        "formamidinesulfinic acid r",
        "[O!H0]-[S,Se,Te;v4]-[C]=[O,N]",
        "==-",
        "",
    ),
    ("isocyanide f", "[C-0!H0]#[N+0]", "#", "-+"),
    ("isocyanide r", "[N+!H0]#[C-]", "#", "-+"),
    ("phosphonic acid f", "[OH]-[PX3H0]", "=", ""),
    ("phosphonic acid r", "[PX4H]=[O]", "-", ""),
];

/// `TautomerTransform`.
struct Transform {
    query: Query,
    bond_types: Vec<BondType>,
    charges: Vec<i32>,
}

fn transforms() -> &'static [Transform] {
    static T: OnceLock<Vec<Transform>> = OnceLock::new();
    T.get_or_init(|| {
        TRANSFORMS
            .iter()
            .map(|&(_, smarts, bonds, charges)| Transform {
                query: smarts_match::parse(smarts),
                // `stringToBondType`
                bond_types: bonds
                    .chars()
                    .filter_map(|c| match c {
                        '-' => Some(BondType::Single),
                        '=' => Some(BondType::Double),
                        '#' => Some(BondType::Triple),
                        ':' => Some(BondType::Aromatic),
                        _ => None,
                    })
                    .collect(),
                // `stringToCharge`
                charges: charges
                    .chars()
                    .map(|c| match c {
                        '+' => 1,
                        '0' => 0,
                        '-' => -1,
                        _ => panic!("Charge symbol not recognised."),
                    })
                    .collect(),
            })
            .collect()
    })
}

/// `getDefaultTautomerScoreSubstructs`: SMARTS and score.
const SCORE_TERMS: [(&str, i32); 12] = [
    ("[#6]1([#6]=[#6][#6]([#6]=[#6]1)=,:[N,S,O])=,:[N,S,O]", 25),
    ("[#6]=[N][OH]", 4),
    ("[#6]=,:[#8]", 2),
    ("[#7]=,:[#8]", 2),
    ("[#15]=,:[#8]", 2),
    ("[C]=[!#1;!#6]", 1),
    ("[C](=[!#1;!#6])[!#1;!#6]", 2),
    ("[c]=!@[N]", -1),
    ("[CX4H3]", 1),
    ("[#7]C(=[NR0])[#7H0]", 1),
    ("[#7;R][#6;R]([N])=[#7;R]", 2),
    ("[#6]=[N+]([O-])[OH]", -4),
];

fn score_terms() -> &'static [(Query, i32)] {
    static T: OnceLock<Vec<(Query, i32)>> = OnceLock::new();
    T.get_or_init(|| {
        SCORE_TERMS
            .iter()
            .map(|&(s, v)| (smarts_match::parse(s), v))
            .collect()
    })
}

/// `TautomerScoringFunctions::scoreRings`.
pub(crate) fn score_rings(mol: &Mol) -> i32 {
    let mut score = 0;
    let Some(ri) = mol.rings.as_ref() else {
        return 0;
    };
    for bring in &ri.bond_rings {
        let mut all_c = true;
        let mut all_aromatic = true;
        for &b in bring {
            let bond = &mol.bonds[b];
            if !bond.aromatic {
                all_aromatic = false;
                break;
            }
            if mol.atoms[bond.begin].anum != 6 || mol.atoms[bond.end].anum != 6 {
                all_c = false;
            }
        }
        if all_aromatic {
            score += 100;
            if all_c {
                score += 150;
            }
        }
    }
    score
}

/// `TautomerScoringFunctions::scoreSubstructs` with the default terms.
pub(crate) fn score_substructs(mol: &Mol) -> i32 {
    score_terms()
        .iter()
        .map(|(q, v)| smarts_match::substruct_match(q, mol).len() as i32 * v)
        .sum()
}

/// `TautomerScoringFunctions::scoreHeteroHs`.
pub(crate) fn score_hetero_hs(mol: &Mol) -> i32 {
    let mut score = 0;
    for a in 0..mol.atoms.len() {
        if matches!(mol.atoms[a].anum, 15 | 16 | 34 | 52) {
            score -= mol.total_num_hs(a) as i32;
        }
    }
    score
}

/// `TautomerScoringFunctions::scoreTautomer`.
pub(crate) fn score_tautomer(mol: &Mol) -> i32 {
    score_rings(mol) + score_substructs(mol) + score_hetero_hs(mol)
}

/// `TautomerEnumeratorStatus`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RdkitTautomerStatus {
    Completed = 0,
    MaxTautomersReached = 1,
    MaxTransformsReached = 2,
    Canceled = 3,
}

/// `Tautomer`.
pub(crate) struct Tautomer {
    pub mol: Mol,
    kekulized: Mol,
    n_modified_atoms: usize,
    n_modified_bonds: usize,
    done: bool,
    /// `_StereochemDone` is a non-computed property (see [`smiles`]).
    flag_kept: bool,
}

/// `TautomerEnumeratorResult`.
pub(crate) struct EnumerateResult {
    /// `d_tautomers`: canonical SMILES (as stored) to tautomer.
    pub tautomers: BTreeMap<String, Tautomer>,
    pub status: RdkitTautomerStatus,
    pub modified_atoms: Vec<bool>,
    pub modified_bonds: Vec<bool>,
}

/// The enumerator's settings (`CleanupParameters` defaults).
#[derive(Clone, Copy, Debug)]
pub(crate) struct Settings {
    pub max_tautomers: usize,
    pub max_transforms: usize,
    pub reassign_stereo: bool,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            max_tautomers: 1000,
            max_transforms: 1000,
            reassign_stereo: true,
        }
    }
}

/// `MolToSmiles(mol)`; `flag_kept`: the molecule's `_StereochemDone` is a
/// non-computed property (set by `setTautomerStereoAndIsoHs` without
/// `reassignStereo`), which edited fragment copies keep.
fn smiles(mol: &Mol, flag_kept: bool) -> Result<String, RdkitSmilesError> {
    if flag_kept {
        write::mol_to_smiles_flag_kept(mol, &RdkitSmilesParams::default())
    } else {
        write::mol_to_smiles(mol, &RdkitSmilesParams::default())
    }
}

/// `MolOps::Kekulize(mol, markAtomsBonds=false, canonical=true)`.
fn kekulize_canonical(mol: &mut Mol) -> Result<(), RdkitSmilesError> {
    let ranks = rank_mol_atoms(mol);
    kekulize_full(mol, Some(&ranks), false)
}

/// Whether bond `b` of `taut` counts as a ring bond for the bond stereo
/// contract.
fn is_ring_bond(taut: &Mol, b: usize) -> bool {
    let bond = &taut.bonds[b];
    taut.num_bond_rings(b) != 0
        || (taut.num_atom_rings(bond.begin) != 0 && taut.num_atom_rings(bond.end) != 0)
}

/// `TautomerEnumerator::setTautomerStereoAndIsoHs` with
/// `removeSp3Stereo=removeBondStereo=removeIsotopicHs=true`.
fn set_tautomer_stereo(
    mol: &Mol,
    taut: &mut Mol,
    modified_atoms: &[bool],
    modified_bonds: &[bool],
    reassign_stereo: bool,
) -> bool {
    let mut modified = false;
    for a in 0..mol.atoms.len() {
        if !modified_atoms[a] {
            continue;
        }
        let t = &mut taut.atoms[a];
        modified |= t.chiral != ChiralTag::Unspecified;
        t.chiral = ChiralTag::Unspecified;
        t.chiral_perm = 0;
        t.cip_code = None;
    }
    for b in 0..mol.bonds.len() {
        if !modified_bonds[b] {
            continue;
        }
        let mut clear_dirs = Vec::new();
        let bond = &mol.bonds[b];
        if bond.bt == BondType::Double && bond.stereo > BondStereo::Any {
            for atom in [bond.begin, bond.end] {
                for &ob in &mol.atom_bonds[atom] {
                    if mol.bonds[ob].dir.is_set() {
                        clear_dirs.push(ob);
                    }
                }
            }
        }
        let double = taut.bonds[b].bt == BondType::Double;
        let ring = double && is_ring_bond(taut, b);
        let target = if double && !ring {
            BondStereo::Any
        } else {
            BondStereo::None
        };
        modified |= taut.bonds[b].stereo != target;
        taut.bonds[b].stereo = target;
        taut.bonds[b].stereo_atoms.clear();
        for bi in clear_dirs {
            taut.bonds[bi].dir = super::mol::BondDir::None;
        }
    }
    if reassign_stereo {
        legacy_stereo_perception(taut, true, false);
        for b in 0..taut.bonds.len() {
            if !modified_bonds[b] {
                continue;
            }
            if taut.bonds[b].bt != BondType::Double {
                taut.bonds[b].stereo = BondStereo::None;
                taut.bonds[b].stereo_atoms.clear();
                continue;
            }
            taut.bonds[b].stereo = if is_ring_bond(taut, b) {
                BondStereo::None
            } else {
                BondStereo::Any
            };
            taut.bonds[b].stereo_atoms.clear();
        }
    }
    modified
}

/// The first key after `key` (`std::map` iteration).
fn next_key(map: &BTreeMap<String, Tautomer>, key: &str) -> Option<String> {
    map.range::<str, _>((Excluded(key), Unbounded))
        .next()
        .map(|(k, _)| k.clone())
}

/// The partial `sanitizeMol` RDKit runs on each product:
/// `SANITIZE_KEKULIZE | SANITIZE_SETAROMATICITY | SANITIZE_SETCONJUGATION |
/// SANITIZE_SETHYBRIDIZATION | SANITIZE_ADJUSTHS`. `Ok(false)`: a
/// `KekulizeException` (the product is skipped).
fn sanitize_product(product: &mut Mol) -> Result<bool, RdkitSmilesError> {
    product.update_property_cache(false)?;
    match kekulize_full(product, None, true) {
        Ok(()) => {}
        Err(RdkitSmilesError::Sanitization(msg)) if msg == CANT_KEKULIZE => return Ok(false),
        Err(e) => return Err(e),
    }
    aromaticity::set_aromaticity(product);
    set_conjugation(product);
    set_hybridization(product);
    adjust_hs(product)?;
    Ok(true)
}

/// `TautomerEnumerator::enumerate(mol)`; `mol` as `MolFromSmiles` leaves
/// it.
pub(crate) fn enumerate(
    mol: &Mol,
    settings: &Settings,
) -> Result<EnumerateResult, RdkitSmilesError> {
    let transforms = transforms();
    let smi = smiles(mol, false)?;
    let mut taut = mol.clone();
    if (0..taut.atoms.len()).any(|a| taut.needs_update_property_cache(a)) {
        taut.update_property_cache(false)?;
    }
    if taut.rings.is_none() {
        taut.find_rings()?;
    }
    let mut kekulized = taut.clone();
    kekulize_canonical(&mut kekulized)?;
    let mut res = EnumerateResult {
        tautomers: BTreeMap::new(),
        status: RdkitTautomerStatus::Completed,
        modified_atoms: vec![false; mol.atoms.len()],
        modified_bonds: vec![false; mol.bonds.len()],
    };
    res.tautomers.insert(
        smi,
        Tautomer {
            mol: taut,
            kekulized,
            n_modified_atoms: 0,
            n_modified_bonds: 0,
            done: false,
            flag_kept: false,
        },
    );
    let mut completed = false;
    let mut bail_out = false;
    let mut n_transforms = 0usize;
    while !completed && !bail_out {
        let mut cursor = res.tautomers.keys().next().cloned();
        while let Some(key) = cursor {
            if !res.tautomers[&key].done {
                let kmol = res.tautomers[&key].kekulized.clone();
                for transform in transforms {
                    if bail_out {
                        break;
                    }
                    let matches = smarts_match::substruct_match(&transform.query, &kmol);
                    if matches.is_empty() {
                        continue;
                    }
                    n_transforms += 1;
                    for m in &matches {
                        if n_transforms >= settings.max_transforms {
                            res.status = RdkitTautomerStatus::MaxTransformsReached;
                            bail_out = true;
                        } else if res.tautomers.len() >= settings.max_tautomers {
                            res.status = RdkitTautomerStatus::MaxTautomersReached;
                            bail_out = true;
                        }
                        if bail_out {
                            break;
                        }
                        let mut product = kmol.clone();
                        let first = m[0];
                        let last = *m.last().expect("non-empty match");
                        res.modified_atoms[first] = true;
                        res.modified_atoms[last] = true;
                        let first_hs = product.total_num_hs(first);
                        product.atoms[first].num_explicit_hs = first_hs.saturating_sub(1);
                        let last_hs = product.total_num_hs(last);
                        product.atoms[last].num_explicit_hs = last_hs + 1;
                        product.atoms[first].no_implicit = true;
                        product.atoms[last].no_implicit = true;
                        for (bi, qb) in (0..transform.query.num_bonds()).enumerate() {
                            let (qa, qz) = transform.query.bond_atoms(qb);
                            let b = product
                                .bond_between(m[qa], m[qz])
                                .expect("required bond not found");
                            if !transform.bond_types.is_empty() {
                                product.bonds[b].bt = transform.bond_types[bi];
                            } else {
                                match product.bonds[b].bt {
                                    BondType::Single => product.bonds[b].bt = BondType::Double,
                                    BondType::Double => product.bonds[b].bt = BondType::Single,
                                    _ => {}
                                }
                            }
                            res.modified_bonds[b] = true;
                        }
                        if !transform.charges.is_empty() {
                            for (ci, &a) in m.iter().enumerate() {
                                product.atoms[a].charge += transform.charges[ci];
                            }
                        }
                        if !sanitize_product(&mut product)? {
                            continue;
                        }
                        set_tautomer_stereo(
                            mol,
                            &mut product,
                            &res.modified_atoms,
                            &res.modified_bonds,
                            settings.reassign_stereo,
                        );
                        let tsmiles = smiles(&product, !settings.reassign_stereo)?;
                        if res.tautomers.contains_key(&tsmiles) {
                            continue;
                        }
                        for i in 0..mol.bonds.len() {
                            if mol.bonds[i].bt != product.bonds[i].bt && !res.modified_bonds[i] {
                                res.modified_bonds[i] = true;
                            }
                        }
                        let mut kekulized_product = product.clone();
                        kekulize_canonical(&mut kekulized_product)?;
                        let n_atoms = res.modified_atoms.iter().filter(|&&x| x).count();
                        let n_bonds = res.modified_bonds.iter().filter(|&&x| x).count();
                        res.tautomers.insert(
                            tsmiles,
                            Tautomer {
                                mol: product,
                                kekulized: kekulized_product,
                                n_modified_atoms: n_atoms,
                                n_modified_bonds: n_bonds,
                                done: false,
                                flag_kept: !settings.reassign_stereo,
                            },
                        );
                    }
                }
                res.tautomers.get_mut(&key).expect("current tautomer").done = true;
            }
            cursor = next_key(&res.tautomers, &key);
        }
        completed = true;
        let max_atoms = res.modified_atoms.iter().filter(|&&x| x).count();
        let max_bonds = res.modified_bonds.iter().filter(|&&x| x).count();
        let mut cursor = res.tautomers.keys().next().cloned();
        while let Some(key) = cursor {
            let t = res.tautomers.get_mut(&key).expect("current tautomer");
            if !t.done {
                completed = false;
            }
            if (t.n_modified_atoms < max_atoms || t.n_modified_bonds < max_bonds)
                && set_tautomer_stereo(
                    mol,
                    &mut t.mol,
                    &res.modified_atoms,
                    &res.modified_bonds,
                    settings.reassign_stereo,
                )
            {
                let mut stored = res.tautomers.remove(&key).expect("current tautomer");
                let after = next_key(&res.tautomers, &key);
                stored.n_modified_atoms = max_atoms;
                stored.n_modified_bonds = max_bonds;
                let new_key = smiles(&stored.mol, stored.flag_kept)?;
                if res.tautomers.contains_key(&new_key) {
                    cursor = after;
                } else {
                    res.tautomers.insert(new_key.clone(), stored);
                    cursor = Some(new_key);
                }
            } else {
                cursor = next_key(&res.tautomers, &key);
            }
        }
        if bail_out
            && res.tautomers.len() < settings.max_tautomers
            && res.status == RdkitTautomerStatus::MaxTautomersReached
        {
            res.status = RdkitTautomerStatus::Completed;
            bail_out = false;
        }
    }
    Ok(res)
}

/// `TautomerEnumerator::canonicalize(mol)` (default scoring): the
/// canonical tautomer, with `assignStereochemistry(cleanIt, force)` run on
/// it, and its score.
pub(crate) fn canonicalize(mol: &Mol) -> Result<(Mol, i32), RdkitSmilesError> {
    let settings = Settings {
        reassign_stereo: false,
        ..Settings::default()
    };
    let res = enumerate(mol, &settings)?;
    let mut best: Option<(&String, &Tautomer, i32)> = None;
    if res.tautomers.len() == 1 {
        let (k, t) = res.tautomers.iter().next().expect("one tautomer");
        best = Some((k, t, score_tautomer(&t.mol)));
    } else {
        for (k, t) in &res.tautomers {
            let score = score_tautomer(&t.mol);
            match best {
                Some((bk, _, bs)) if score < bs || (score == bs && k >= bk) => {}
                _ => best = Some((k, t, score)),
            }
        }
    }
    let (_, t, score) = best.expect("at least the input tautomer");
    let mut out = t.mol.clone();
    legacy_stereo_perception(&mut out, true, false);
    Ok((out, score))
}
