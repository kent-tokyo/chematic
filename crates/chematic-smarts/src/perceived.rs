//! SMARTS matching against the target's perceived aromaticity (issue #635).
//!
//! Daylight SMARTS semantics assume the target molecule's aromaticity has
//! been perceived: `c` matches every aromatic carbon, `C` only aliphatic ones,
//! and an unspecified bond matches single or aromatic bonds. The core
//! [`crate::find_matches`] family matches the molecule exactly as given, so a
//! benzene written in Kekulé form (`C1=CC=CC=C1`) has no `c` atom there.
//!
//! These functions match against the molecule's RDKit-parity aromatic view
//! ([`chematic_perception::apply_aromaticity_rdkit_parity_shared`]) — the
//! same memoized view the RDKit-compatible descriptors and fingerprints use.
//! That view keeps every atom and bond index, so the returned matches index
//! the molecule that was passed in. When perception fails (a molecule the
//! RDKit-parity engine cannot kekulize), matching falls back to the molecule
//! as given, which is exactly the core behavior.
//!
//! These are the SMARTS semantics of the language bindings (Python, WASM);
//! the core functions are unchanged.

use rustc_hash::FxHashMap;

use chematic_core::{AtomIdx, Molecule};

use crate::match_vf2::{MatchConfig, find_matches_with_config, has_match_with_config};
use crate::query::QueryMolecule;

/// Run `f` on the molecule SMARTS matching should see: the RDKit-parity
/// aromatic view of `mol`, or `mol` itself if perception fails.
pub fn with_perceived_target<R>(mol: &Molecule, f: impl FnOnce(&Molecule) -> R) -> R {
    // When the view would be an exact copy of `mol`, match `mol` itself and
    // skip building (and copying) the view.
    if chematic_perception::rdkit_parity_view_is_identity(mol) {
        return f(mol);
    }
    let view = chematic_perception::apply_aromaticity_rdkit_parity_shared(mol);
    match view.as_ref() {
        Ok(perceived) => f(perceived),
        Err(_) => f(mol),
    }
}

/// [`find_matches_with_config`] against the perceived aromatic view of `mol`
/// (see the module docs). Atom indices refer to `mol`.
pub fn find_matches_perceived(
    query: &QueryMolecule,
    mol: &Molecule,
    config: &MatchConfig,
) -> Vec<FxHashMap<usize, AtomIdx>> {
    with_perceived_target(mol, |target| {
        find_matches_with_config(query, target, config)
    })
}

/// Whether `query` matches the perceived aromatic view of `mol` at least once.
/// `max_matches` and `uniquify` in `config` are ignored.
pub fn has_match_perceived(query: &QueryMolecule, mol: &Molecule, config: &MatchConfig) -> bool {
    with_perceived_target(mol, |target| has_match_with_config(query, target, config))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse_smarts;
    use chematic_smiles::parse;

    fn sets(query: &str, smiles: &str) -> Vec<Vec<u32>> {
        let q = parse_smarts(query).unwrap();
        let mol = parse(smiles).unwrap();
        let mut out: Vec<Vec<u32>> = find_matches_perceived(&q, &mol, &MatchConfig::default())
            .into_iter()
            .map(|m| {
                let mut v: Vec<u32> = m.values().map(|a| a.0).collect();
                v.sort_unstable();
                v
            })
            .collect();
        out.sort();
        out
    }

    #[test]
    fn kekule_and_aromatic_spellings_match_alike() {
        for q in [
            "c",
            "c1ccccc1",
            "[#6]=[#6]",
            "C",
            "[N,O]",
            "[nH]",
            "c:c",
            "C=O",
        ] {
            assert_eq!(sets(q, "C1=CC=CC=C1"), sets(q, "c1ccccc1"), "{q} benzene");
            assert_eq!(sets(q, "C1=CNC=C1"), sets(q, "c1c[nH]cc1"), "{q} pyrrole");
            assert_eq!(
                sets(q, "O=C1C=CC(=O)C=C1"),
                sets(q, "O=C1C=CC(=O)C=C1"),
                "{q} quinone"
            );
        }
        // Positive and negative truth table on a Kekulé input.
        assert_eq!(sets("c", "C1=CC=CC=C1").len(), 6);
        assert!(sets("[#6]=[#6]", "C1=CC=CC=C1").is_empty());
        assert!(sets("C", "C1=CC=CC=C1").is_empty());
        assert_eq!(sets("[nH]", "C1=CC=CN1"), vec![vec![4]]);
        // Non-aromatic rings keep aliphatic semantics.
        assert_eq!(sets("[#6]=[#6]", "C1=CCCCC1").len(), 1);
        assert_eq!(sets("C", "C1=CCCCC1").len(), 6);
    }

    #[test]
    fn indices_refer_to_the_input_molecule() {
        // Pyridine written Kekulé with the N last: the aromatic n is atom 5.
        assert_eq!(sets("n", "C1=CC=CC=N1"), vec![vec![5]]);
        assert!(has_match_perceived(
            &parse_smarts("c1ccncc1").unwrap(),
            &parse("C1=CC=CC=N1").unwrap(),
            &MatchConfig::default()
        ));
    }
}
