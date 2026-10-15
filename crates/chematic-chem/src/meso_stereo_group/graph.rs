use super::MesoStereoGroupError as Error;
use chematic_core::ez_markers::{EzFact, read_ez_facts};
use chematic_core::{AtomIdx, BondIdx, BondOrder, Chirality, Molecule, STEREO_H_SENTINEL};

#[derive(Clone, PartialEq, Eq)]
pub(super) struct Color {
    element: u8,
    isotope: Option<u16>,
    charge: i8,
    hydrogen: u8,
    aromatic: bool,
    atom_map: Option<u16>,
    tag: Option<u16>,
    group: Option<usize>,
    tetra: bool,
    degree: usize,
}
#[derive(Clone, Copy)]
pub(super) struct Edge {
    pub index: BondIdx,
    order: u8,
    group: Option<usize>,
    geometric: bool,
}
impl Edge {
    pub fn signature(self) -> (u8, Option<usize>, bool) {
        (self.order, self.group, self.geometric)
    }
}
pub(super) struct Graph<'a> {
    pub mol: &'a Molecule,
    pub colors: Vec<Color>,
    pub adj: Vec<Vec<Option<Edge>>>,
    pub frames: Vec<Option<[u32; 4]>>,
    pub facts: Vec<EzFact>,
    pub fact_by_bond: Vec<Option<usize>>,
}
fn frame(mol: &Molecule, idx: AtomIdx) -> Result<Option<[u32; 4]>, Error> {
    let a = mol.atom(idx);
    if a.wildcard
        || mol.r_group_label(idx).is_some()
        || a.chirality.nontetrahedral().is_some()
        || mol.is_allene_stereo(idx)
    {
        return Err(Error::UnsupportedAtom(idx));
    }
    let h = chematic_core::implicit_hcount(mol, idx);
    if a.chirality == Chirality::None {
        if !a.aromatic
            && h <= 1
            && mol.degree(idx) + h as usize == 4
            && mol.neighbors(idx).all(|(_, b)| {
                matches!(
                    mol.bond(b).order,
                    BondOrder::Single | BondOrder::Up | BondOrder::Down
                )
            })
        {
            return Err(Error::UnspecifiedTetrahedralAtom(idx));
        }
        return Ok(None);
    }
    let stored = mol.stereo_neighbor_order(idx).unwrap_or_default();
    let mut expected: Vec<u32> = mol.neighbors(idx).map(|(a, _)| a.0).collect();
    expected.extend(std::iter::repeat_n(STEREO_H_SENTINEL, h as usize));
    let mut sorted = stored.to_vec();
    sorted.sort_unstable();
    expected.sort_unstable();
    if expected.len() != 4 || sorted != expected || sorted.windows(2).any(|v| v[0] == v[1]) {
        return Err(Error::IncompleteFrame(idx));
    }
    Ok(Some(stored.try_into().expect("four slots validated")))
}
fn order(order: BondOrder, index: BondIdx) -> Result<u8, Error> {
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down => Ok(1),
        BondOrder::Double => Ok(2),
        BondOrder::Triple => Ok(3),
        BondOrder::Aromatic => Ok(4),
        _ => Err(Error::UnsupportedBond(index)),
    }
}
fn geometry(mol: &Molecule) -> Result<Vec<EzFact>, Error> {
    if mol.has_bond_directions() {
        return Err(Error::UnsupportedDirections);
    }
    let Some((facts, false)) = read_ez_facts(mol) else {
        return Err(Error::UnsupportedDirections);
    };
    // No marker may be ignored as a lone wedge or incomplete E/Z declaration.
    for (_, bond) in mol.bonds() {
        if !matches!(bond.order, BondOrder::Up | BondOrder::Down) {
            continue;
        }
        let accounted = facts.iter().any(|f| {
            f.ends
                .iter()
                .any(|e| bond.atom1 == e.atom || bond.atom2 == e.atom)
        });
        if !accounted {
            return Err(Error::UnsupportedDirections);
        }
    }
    Ok(facts)
}
impl<'a> Graph<'a> {
    pub fn new(mol: &'a Molecule) -> Result<Self, Error> {
        let facts = geometry(mol)?;
        let n = mol.atom_count();
        let mut atom_groups = vec![None; n];
        let mut bond_groups = vec![None; mol.bond_count()];
        for (i, g) in mol.stereo_groups().iter().enumerate() {
            for a in &g.atom_indices {
                atom_groups[a.0 as usize] = Some(i);
            }
            for b in mol.stereo_group_bond_indices(i).unwrap_or_default() {
                bond_groups[b.0 as usize] = Some(i);
            }
        }
        let mut fact_by_bond = vec![None; mol.bond_count()];
        for (i, f) in facts.iter().enumerate() {
            fact_by_bond[f.double.0 as usize] = Some(i);
        }
        let mut colors = Vec::new();
        let mut frames = Vec::new();
        for (idx, a) in mol.atoms() {
            let f = frame(mol, idx)?;
            colors.push(Color {
                element: a.element.atomic_number(),
                isotope: a.isotope,
                charge: a.charge,
                hydrogen: chematic_core::implicit_hcount(mol, idx),
                aromatic: a.aromatic,
                atom_map: a.atom_map,
                tag: mol.atom_tag(idx).map(|v| v.get()),
                group: atom_groups[idx.0 as usize],
                tetra: f.is_some(),
                degree: mol.degree(idx),
            });
            frames.push(f);
        }
        let mut adj = vec![vec![None; n]; n];
        for (idx, bond) in mol.bonds() {
            let a = bond.atom1.0 as usize;
            let b = bond.atom2.0 as usize;
            if a == b || adj[a][b].is_some() {
                return Err(Error::UnsupportedBond(idx));
            }
            let edge = Edge {
                index: idx,
                order: order(bond.order, idx)?,
                group: bond_groups[idx.0 as usize],
                geometric: fact_by_bond[idx.0 as usize].is_some(),
            };
            adj[a][b] = Some(edge);
            adj[b][a] = Some(edge);
        }
        Ok(Self {
            mol,
            colors,
            adj,
            frames,
            facts,
            fact_by_bond,
        })
    }
}
