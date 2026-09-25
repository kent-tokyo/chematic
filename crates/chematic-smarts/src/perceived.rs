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
use crate::query::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, QueryMolecule};

/// Whether every primitive of `query` evaluates identically on a molecule
/// and on its RDKit-parity aromatic view.
///
/// The view keeps every atom (element, charge, isotope, chirality) and every
/// bond's endpoints, hence each atom's degree and its membership in some
/// ring (a graph property); it may change aromatic flags, bond orders,
/// implicit hydrogen counts and, through the ranks that break SSSR ties,
/// which rings are chosen. Queries built only from atomic number, element,
/// charge, isotope, chirality, wildcard, degree and ring membership atom
/// primitives (and recursive SMARTS of such queries) joined by `~` bonds
/// therefore match both identically, embedding for embedding.
pub fn is_aromaticity_insensitive(query: &QueryMolecule) -> bool {
    fn atom_ok(q: &AtomQuery) -> bool {
        match q {
            AtomQuery::Primitive(p) => match p {
                AtomPrimitive::AtomicNum(_)
                | AtomPrimitive::Symbol(_)
                | AtomPrimitive::Charge(_)
                | AtomPrimitive::Isotope(_)
                | AtomPrimitive::Chirality(_)
                | AtomPrimitive::Wildcard
                | AtomPrimitive::Degree(_)
                | AtomPrimitive::RingMembership(_) => true,
                AtomPrimitive::Recursive(sub) => is_aromaticity_insensitive(sub),
                _ => false,
            },
            AtomQuery::And(a, b) | AtomQuery::Or(a, b) => atom_ok(a) && atom_ok(b),
            AtomQuery::Not(a) => atom_ok(a),
        }
    }
    fn bond_ok(q: &BondQuery) -> bool {
        match q {
            BondQuery::Primitive(BondPrimitive::Any) => true,
            BondQuery::Primitive(_) | BondQuery::Any => false,
            BondQuery::And(a, b) | BondQuery::Or(a, b) => bond_ok(a) && bond_ok(b),
            BondQuery::Not(a) => bond_ok(a),
        }
    }
    query.atoms.iter().all(|a| atom_ok(&a.query)) && query.bonds.iter().all(|b| bond_ok(&b.query))
}

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
    if is_aromaticity_insensitive(query) {
        return find_matches_with_config(query, mol, config);
    }
    with_perceived_target(mol, |target| {
        find_matches_with_config(query, target, config)
    })
}

/// The matches of [`find_matches_perceived`] as target-atom sets: each match's
/// atom indices sorted ascending, in the same order as the matches. Built
/// directly from the embeddings (no per-match maps) when `config` has no
/// `max_matches`; otherwise identical to mapping [`find_matches_perceived`].
pub fn find_match_atom_sets_perceived(
    query: &QueryMolecule,
    mol: &Molecule,
    config: &MatchConfig,
) -> Vec<Vec<usize>> {
    if config.max_matches.is_some() {
        return find_matches_perceived(query, mol, config)
            .into_iter()
            .map(|m| {
                let mut v: Vec<usize> = m.values().map(|idx| idx.0 as usize).collect();
                v.sort_unstable();
                v
            })
            .collect();
    }
    let run = |target: &Molecule| {
        let mut out: Vec<Vec<usize>> = Vec::new();
        let mut seen: rustc_hash::FxHashSet<Vec<usize>> = rustc_hash::FxHashSet::default();
        let _ = crate::match_vf2::for_each_embedding(query, target, config, |m| {
            let mut v: Vec<usize> = m.iter().map(|&t| t as usize).collect();
            v.sort_unstable();
            if !config.uniquify || seen.insert(v.clone()) {
                out.push(v);
            }
        });
        out
    };
    if is_aromaticity_insensitive(query) {
        return run(mol);
    }
    with_perceived_target(mol, run)
}

/// Whether `query` matches the perceived aromatic view of `mol` at least once.
/// `max_matches` and `uniquify` in `config` are ignored.
pub fn has_match_perceived(query: &QueryMolecule, mol: &Molecule, config: &MatchConfig) -> bool {
    if is_aromaticity_insensitive(query) {
        return has_match_with_config(query, mol, config);
    }
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
    fn aromaticity_insensitive_queries_are_classified_conservatively() {
        for (q, expected) in [
            ("[#7;R]", true),
            ("*~*~*~*~*~*", true),
            ("[#6]~[#7]", true),
            ("[!#6;!#1]~*~[!#6;!#1]", true),
            ("[$([#8]~[#6])]", true),
            ("[#6;D3;+0]", true),
            ("c1ccccc1", false),
            ("[OH]", false),
            ("C(=O)N", false),
            ("[#6]-[#7]", false),
            ("[#6]:[#7]", false),
            ("[#6][#7]", false),
            ("[R2]", false),
            ("[r5]", false),
            ("[#6]@[#6]", false),
            ("[$(C=O)]", false),
            ("[#6;H1]", false),
            ("[#6;X3]", false),
        ] {
            let query = parse_smarts(q).unwrap();
            assert_eq!(is_aromaticity_insensitive(&query), expected, "{q}");
        }
    }

    #[test]
    fn insensitive_queries_match_the_view_like_the_molecule() {
        let cfg = MatchConfig::default();
        for q in ["[#7;R]", "*~*~*~*~*~*", "[#6]~[#7]", "[!#6;!#1]~*~[!#6;!#1]"] {
            let query = parse_smarts(q).unwrap();
            for smi in ["C1=CC=CC=C1N", "O=C1c2ccccc2C(=O)N1C", "C1=CNC=C1", "c1ccc2[nH]ccc2c1"] {
                let mol = parse(smi).unwrap();
                let via_view = with_perceived_target(&mol, |t| find_matches_with_config(&query, t, &cfg));
                assert_eq!(find_matches_perceived(&query, &mol, &cfg), via_view, "{q} {smi}");
            }
        }
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
