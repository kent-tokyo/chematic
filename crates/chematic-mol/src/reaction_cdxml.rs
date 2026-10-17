//! Single-step CDXML output with explicit component roles and caller layout.
use std::fmt::Write;

use chematic_core::{BondOrder, Chirality, Molecule};
use chematic_rxn::Reaction;

pub const MAX_REACTION_CDXML_COMPONENTS: usize = 32;
pub const MAX_REACTION_CDXML_ATOMS: usize = 1024;

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ReactionCdxmlError {
    InvalidInput(String),
    Unsupported {
        component: usize,
        reason: &'static str,
    },
    Allocation,
}
impl std::fmt::Display for ReactionCdxmlError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidInput(s) => write!(f, "reaction CDXML: {s}"),
            Self::Unsupported { component, reason } => {
                write!(f, "reaction CDXML component {component}: {reason}")
            }
            Self::Allocation => f.write_str("reaction CDXML allocation failed"),
        }
    }
}
impl std::error::Error for ReactionCdxmlError {}
fn invalid(message: &str) -> ReactionCdxmlError {
    ReactionCdxmlError::InvalidInput(message.into())
}
fn components(r: &Reaction) -> impl Iterator<Item = &Molecule> {
    r.reactants.iter().chain(&r.agents).chain(&r.products)
}
fn point_valid(point: [f64; 2]) -> bool {
    point
        .into_iter()
        .all(|value| value.is_finite() && value.abs() <= 1_000_000.)
}
fn check_graph(mol: &Molecule, component: usize) -> Result<(), ReactionCdxmlError> {
    let unsupported = |reason| ReactionCdxmlError::Unsupported { component, reason };
    if !mol.stereo_groups().is_empty() {
        return Err(unsupported("enhanced stereo groups are unsupported"));
    }
    for (idx, atom) in mol.atoms() {
        if atom.wildcard || mol.r_group_label(idx).is_some() || atom.atom_map.is_some() {
            return Err(unsupported(
                "query atoms, R groups and atom maps are unsupported",
            ));
        }
        if atom.aromatic
            || atom.chirality != Chirality::None
            || atom.cip_code.is_some()
            || mol.is_allene_stereo(idx)
        {
            return Err(unsupported(
                "aromatic or stereochemical atom annotations are unsupported",
            ));
        }
    }
    for (idx, atom) in mol.atoms() {
        if atom.element.atomic_number() == 1 && mol.degree(idx) != 1 {
            return Err(unsupported(
                "isolated or bridging hydrogen atoms are unsupported",
            ));
        }
        if atom
            .hydrogen_count
            .is_some_and(|h| h != chematic_core::valence::valence_inferred_hcount(mol, idx))
        {
            return Err(unsupported(
                "bracket hydrogen count differs from native normal-valence inference",
            ));
        }
    }
    for (_, bond) in mol.bonds() {
        match bond.order {
            BondOrder::Single | BondOrder::Triple => {}
            BondOrder::Double if mol.degree(bond.atom1) == 1 || mol.degree(bond.atom2) == 1 => {}
            _ => {
                return Err(unsupported(
                    "only single/triple and terminal double bonds are supported",
                ));
            }
        }
    }
    Ok(())
}
fn validate(
    r: &Reaction,
    coords: &[Vec<[f64; 2]>],
    arrow: [[f64; 2]; 2],
) -> Result<(), ReactionCdxmlError> {
    let count = components(r).count();
    if r.reactants.is_empty() || r.products.is_empty() || count > MAX_REACTION_CDXML_COMPONENTS {
        return Err(invalid(
            "requires reactants and products, with at most 32 total components",
        ));
    }
    if coords.len() != count {
        return Err(invalid(
            "coordinate rows must match reactants, then agents, then products",
        ));
    }
    if !arrow.into_iter().all(point_valid) || arrow[0] == arrow[1] {
        return Err(invalid(
            "arrow requires distinct finite endpoints within +/-1000000",
        ));
    }
    let mut atoms = 0;
    let mut bonds = 0;
    for (i, (mol, points)) in components(r).zip(coords).enumerate() {
        atoms += mol.atom_count();
        bonds += mol.bond_count();
        if mol.atom_count() == 0
            || mol.atom_count() > 256
            || atoms > MAX_REACTION_CDXML_ATOMS
            || bonds > 4096
        {
            return Err(invalid(
                "supports nonempty components of at most 256 atoms, 1024 total atoms and 4096 bonds",
            ));
        }
        if points.len() != mol.atom_count() || !points.iter().copied().all(point_valid) {
            return Err(invalid(
                "coordinates require one finite pair per atom within +/-1000000",
            ));
        }
        check_graph(mol, i)?;
    }
    Ok(())
}
fn fragment(out: &mut String, mol: &Molecule, coords: &[[f64; 2]], id: usize, next: &mut usize) {
    writeln!(out, "<fragment id=\"{id}\">").unwrap();
    let first = *next;
    for ((idx, atom), [x, y]) in mol.atoms().zip(coords) {
        let atom_id = first + idx.0 as usize;
        write!(
            out,
            "<n id=\"{atom_id}\" Element=\"{}\" p=\"{x:e} {y:e}\"",
            atom.element.atomic_number()
        )
        .unwrap();
        if let Some(h) = atom.hydrogen_count {
            write!(out, " NumHydrogens=\"{h}\"").unwrap();
        }
        if atom.charge != 0 {
            write!(out, " Charge=\"{}\"", atom.charge).unwrap();
        }
        if let Some(isotope) = atom.isotope {
            write!(out, " Isotope=\"{isotope}\"").unwrap();
        }
        out.push_str("/>\n");
    }
    *next += mol.atom_count();
    for (_, bond) in mol.bonds() {
        let (begin, end) = (first + bond.atom1.0 as usize, first + bond.atom2.0 as usize);
        writeln!(
            out,
            "<b id=\"{}\" B=\"{begin}\" E=\"{end}\" Order=\"{}\"/>",
            *next,
            bond.order.order_int()
        )
        .unwrap();
        *next += 1;
    }
    out.push_str("</fragment>\n");
}
fn ids(start: usize, count: usize) -> String {
    (start..start + count)
        .map(|id| id.to_string())
        .collect::<Vec<_>>()
        .join(" ")
}

/// Write one reaction step with globally unique fragment/node/bond/arrow IDs.
/// Coordinates are in CDXML Y-down units, ordered reactants, agents, products.
/// The arrow uses tail then head; no layout or role inference is performed.
///
/// Profile: non-aromatic unmapped non-query graphs without stereo, with single,
/// triple or terminal double bonds. Bracket H must match native normal-valence
/// inference; isolated/bridging H is refused. Other profiles fail instead of losing labels
/// or deriving E/Z stereo from caller coordinates. At most 32 components, 256
/// atoms per component, 1024 total atoms, 4096 bonds; coordinates are finite and
/// bounded to +/-1e6. At least one reactant and product are required. Agents are
/// listed as ReactionStepObjectsAboveArrow. Only graph/role/layout data is taken;
/// this does not serialize ReactionDocument conditions or source metadata.
/// ChemDraw application acceptance and automatic depiction are not promised.
pub fn write_reaction_cdxml_checked(
    r: &Reaction,
    coordinates: &[Vec<[f64; 2]>],
    arrow: [[f64; 2]; 2],
) -> Result<String, ReactionCdxmlError> {
    validate(r, coordinates, arrow)?;
    let n = components(r).count();
    let mut out = String::new();
    // All numeric attributes use bounded-width scientific notation. This cap
    // exceeds the validated maximum serialized atom/bond payload.
    out.try_reserve(2 * 1024 * 1024)
        .map_err(|_| ReactionCdxmlError::Allocation)?;
    out.push_str(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<CDXML BondLength=\"30\">\n<page id=\"1\">\n",
    );
    let mut next = 2 + n;
    for (i, (mol, coords)) in components(r).zip(coordinates).enumerate() {
        fragment(&mut out, mol, coords, i + 2, &mut next);
    }
    let [[tx, ty], [hx, hy]] = arrow;
    let arrow_id = next;
    writeln!(out, "<arrow id=\"{arrow_id}\" ArrowheadHead=\"Full\" ArrowheadType=\"Solid\" Head3D=\"{hx:e} {hy:e} 0\" Tail3D=\"{tx:e} {ty:e} 0\" BoundingBox=\"{:e} {:e} {:e} {:e}\"/>", tx.min(hx), ty.min(hy), tx.max(hx), ty.max(hy)).unwrap();
    let reactants = ids(2, r.reactants.len());
    let agents = ids(2 + r.reactants.len(), r.agents.len());
    let products = ids(2 + r.reactants.len() + r.agents.len(), r.products.len());
    writeln!(out, "<scheme id=\"{}\">\n<step id=\"{}\" ReactionStepReactants=\"{reactants}\" ReactionStepProducts=\"{products}\" ReactionStepObjectsAboveArrow=\"{agents}\" ReactionStepArrows=\"{arrow_id}\"/>\n</scheme>\n</page>\n</CDXML>", next+1, next+2).unwrap();
    Ok(out)
}
