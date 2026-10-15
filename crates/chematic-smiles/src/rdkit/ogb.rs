//! OGB molecular graph features, using the same sanitized model as the
//! RDKit-compatible SMILES implementation. Atom indices are never reordered.

use super::mol::{BondStereo, BondType, ChiralTag, Hybridization};
use super::{RdkitSmilesError, parse, sanitize, stereo};
use chematic_core::Molecule;
use std::borrow::Cow;

/// The categorical graph schema used by `ogb.utils.mol.smiles2graph`.
///
/// By default RDKit-removable hydrogen graph atoms are removed, as in OGB's
/// `smiles2graph`. Other atoms retain their relative order. The optional
/// variant can preserve every input atom instead. Bonds are emitted in RDKit
/// bond order, in both directions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OgbGraphFeatures {
    pub node_feat: Vec<[i64; 9]>,
    pub edge_index: Vec<[i64; 2]>,
    pub edge_feat: Vec<[i64; 3]>,
}

/// Produce OGB's nine atom and three bond category indices from one molecule.
/// Hydrogen graph atoms follow OGB's default removal policy.
///
/// Refuses inputs the RDKit-compatible sanitizer cannot model. Derived CIP
/// labels are ignored: OGB reads chiral tags, not the `_CIPCode` property.
/// A `Molecule` carries no source SDF coordinates, so this function cannot
/// derive tetrahedral chirality from 3D geometry.
pub fn ogb_graph_features(mol: &Molecule) -> Result<OgbGraphFeatures, RdkitSmilesError> {
    ogb_graph_features_with_options(mol, true)
}

/// OGB graph features with a choice to retain or remove hydrogen graph atoms.
/// `remove_hs=false` retains the graph's current atoms; it does not add new
/// hydrogen atoms like OGB's `smiles2graph(removeHs=False)` does.
pub fn ogb_graph_features_with_options(
    mol: &Molecule,
    remove_hs: bool,
) -> Result<OgbGraphFeatures, RdkitSmilesError> {
    let source = if mol.atoms().any(|(_, atom)| atom.cip_code.is_some()) {
        let mut copy = mol.clone();
        for (idx, atom) in mol.atoms() {
            if atom.cip_code.is_some() {
                copy.set_cip_code(idx, None);
            }
        }
        Cow::Owned(copy)
    } else {
        Cow::Borrowed(mol)
    };
    let mut model = parse::from_chematic(&source)?;
    if remove_hs {
        sanitize::remove_hs(&mut model, true)?;
    } else {
        sanitize::sanitize_keeping_hs(&mut model)?;
    }
    stereo::legacy_stereo_perception_unflagged(&mut model);

    let mut node_feat = Vec::with_capacity(model.atoms.len());
    for (idx, atom) in model.atoms.iter().enumerate() {
        let atomic_num = if (1..=118).contains(&atom.anum) {
            i64::from(atom.anum - 1)
        } else {
            118 // misc
        };
        let chirality = match atom.chiral {
            ChiralTag::Unspecified => 0,
            ChiralTag::Cw => 1,
            ChiralTag::Ccw => 2,
            _ => 4, // unknown OGB chirality category (misc)
        };
        let degree = model.degree(idx) as u32 + model.total_num_hs(idx);
        let degree = if degree <= 10 { i64::from(degree) } else { 11 };
        let charge = if (-5..=5).contains(&atom.charge) {
            i64::from(atom.charge + 5)
        } else {
            11
        };
        let num_hs = model.total_num_hs(idx);
        let num_hs = if num_hs <= 8 { i64::from(num_hs) } else { 9 };
        let radicals = if atom.radicals <= 4 {
            i64::from(atom.radicals)
        } else {
            5
        };
        let hybridization = match atom.hybrid {
            Hybridization::Sp => 0,
            Hybridization::Sp2 => 1,
            Hybridization::Sp3 => 2,
            Hybridization::Sp3d => 3,
            Hybridization::Sp3d2 => 4,
            _ => 5, // S, SP2D, unspecified, or another OGB misc value
        };
        node_feat.push([
            atomic_num,
            chirality,
            degree,
            charge,
            num_hs,
            radicals,
            hybridization,
            i64::from(atom.aromatic),
            i64::from(model.num_atom_rings(idx) != 0),
        ]);
    }

    let mut edge_index = Vec::with_capacity(model.bonds.len() * 2);
    let mut edge_feat = Vec::with_capacity(model.bonds.len() * 2);
    for bond in &model.bonds {
        let bond_type = match bond.bt {
            BondType::Single => 0,
            BondType::Double => 1,
            BondType::Triple => 2,
            BondType::Aromatic => 3,
            _ => 4, // misc
        };
        let bond_stereo = match bond.stereo {
            BondStereo::None => 0,
            BondStereo::Z => 1,
            BondStereo::E => 2,
            BondStereo::Any => 5,
        };
        let feature = [bond_type, bond_stereo, i64::from(bond.conjugated)];
        edge_index.push([bond.begin as i64, bond.end as i64]);
        edge_feat.push(feature);
        edge_index.push([bond.end as i64, bond.begin as i64]);
        edge_feat.push(feature);
    }
    Ok(OgbGraphFeatures {
        node_feat,
        edge_index,
        edge_feat,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse;

    #[test]
    fn ethanol_matches_ogb_categories_and_directed_edge_order() {
        let graph = ogb_graph_features(&parse("CCO").unwrap()).unwrap();
        assert_eq!(
            graph.node_feat,
            vec![
                [5, 0, 4, 5, 3, 0, 2, 0, 0],
                [5, 0, 4, 5, 2, 0, 2, 0, 0],
                [7, 0, 2, 5, 1, 0, 2, 0, 0],
            ]
        );
        assert_eq!(graph.edge_index, [[0, 1], [1, 0], [1, 2], [2, 1]]);
        assert_eq!(graph.edge_feat, [[0, 0, 0]; 4]);
    }

    #[test]
    fn aromatic_ring_and_no_bond_have_defined_shapes() {
        let ring = ogb_graph_features(&parse("c1ccccc1").unwrap()).unwrap();
        assert_eq!(ring.node_feat, vec![[5, 0, 3, 5, 1, 0, 1, 1, 1]; 6]);
        assert_eq!(ring.edge_index.len(), 12);
        assert_eq!(ring.edge_feat, vec![[3, 0, 1]; 12]);

        let ammonium = ogb_graph_features(&parse("[NH4+]").unwrap()).unwrap();
        assert_eq!(ammonium.node_feat, [[6, 0, 4, 6, 4, 0, 2, 0, 0]]);
        assert!(ammonium.edge_index.is_empty());
        assert!(ammonium.edge_feat.is_empty());
    }

    #[test]
    fn radicals_tetrahedral_stereo_and_alkene_stereo_match_rdkit() {
        let methyl = ogb_graph_features(&parse("[CH3]").unwrap()).unwrap();
        assert_eq!(methyl.node_feat, [[5, 0, 3, 5, 3, 1, 2, 0, 0]]);

        let chiral = ogb_graph_features(&parse("F[C@](Cl)(Br)I").unwrap()).unwrap();
        assert_eq!(chiral.node_feat[1][1], 2); // CHI_TETRAHEDRAL_CCW

        let alkene = ogb_graph_features(&parse("C/C=C/C").unwrap()).unwrap();
        assert_eq!(alkene.edge_feat[2], [1, 2, 0]); // DOUBLE, STEREOE
        assert_eq!(alkene.edge_feat[3], alkene.edge_feat[2]);
    }

    #[test]
    fn graph_hydrogen_atoms_keep_their_indices() {
        let mol = parse("[H]C([H])([H])[H]").unwrap();
        let default = ogb_graph_features(&mol).unwrap();
        assert_eq!(default.node_feat, [[5, 0, 4, 5, 4, 0, 2, 0, 0]]);
        assert!(default.edge_index.is_empty());

        let graph = ogb_graph_features_with_options(&mol, false).unwrap();
        assert_eq!(graph.node_feat.len(), mol.atom_count());
        assert_eq!(graph.edge_index.len(), mol.bond_count() * 2);
    }
}
