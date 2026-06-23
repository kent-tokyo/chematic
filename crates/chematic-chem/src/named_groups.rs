//! Named functional group detection via SMARTS patterns.
//!
//! Returns human-readable group names (e.g. "hydroxyl", "carbonyl") and the
//! atom indices that constitute each match. Multiple matches of the same group
//! (e.g. a dicarboxylic acid) each produce a separate entry.

use std::sync::OnceLock;

use chematic_core::{AtomIdx, Molecule};
use chematic_perception::find_sssr;
use chematic_smarts::{
    MatchConfig, QueryMolecule, find_matches_with_rings_and_config, parse_smarts,
};

/// A single detected functional group with its name and matched atom indices.
pub struct NamedGroup {
    pub name: &'static str,
    pub atoms: Vec<AtomIdx>,
}

static NAMED_GROUPS: &[(&str, &str)] = &[
    ("hydroxyl", "[OX2H1]"),
    ("carbonyl", "[CX3]=[OX1]"),
    ("carboxyl", "[CX3](=[OX1])[OX2H1]"),
    ("aldehyde", "[CX3H1]=[OX1]"),
    ("ketone", "[#6][CX3](=[OX1])[#6]"),
    ("ester", "[#6][CX3](=[OX1])[OX2H0][#6]"),
    ("amide", "[NX3][CX3]=[OX1]"),
    ("primary_amine", "[NX3H2]"),
    ("secondary_amine", "[NX3H1]"),
    ("tertiary_amine", "[NX3H0;!R]"),
    ("nitro", "[#7](=[OX1])=[OX1]"),
    ("nitrile", "[NX1]#[CX2]"),
    ("thiol", "[SX2H1]"),
    ("sulfide", "[SX2H0]"),
    ("sulfoxide", "[SX3]=[OX1]"),
    ("sulfone", "[SX4](=[OX1])=[OX1]"),
    ("halogen", "[F,Cl,Br,I]"),
    ("ether", "[OX2H0]"),
    ("epoxide", "[OX2r3]"),
    ("aromatic", "[a]"),
];

fn parsed_patterns() -> &'static [(&'static str, QueryMolecule)] {
    static CACHE: OnceLock<Vec<(&'static str, QueryMolecule)>> = OnceLock::new();
    CACHE.get_or_init(|| {
        NAMED_GROUPS
            .iter()
            .filter_map(|&(name, smarts)| parse_smarts(smarts).ok().map(|q| (name, q)))
            .collect()
    })
}

/// Detect named functional groups in `mol`.
///
/// Each SMARTS pattern is matched independently; overlapping groups (e.g. a
/// carboxylic acid is detected as both "carboxyl" and "hydroxyl") are returned
/// as separate entries.
pub fn detect_named_functional_groups(mol: &Molecule) -> Vec<NamedGroup> {
    let rings = find_sssr(mol);
    let config = MatchConfig::default();
    let patterns = parsed_patterns();
    let mut result = Vec::new();
    for (name, query) in patterns {
        for m in find_matches_with_rings_and_config(query, mol, &rings, &config) {
            let mut atoms: Vec<AtomIdx> = m.values().copied().collect();
            atoms.sort_unstable();
            result.push(NamedGroup { name, atoms });
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_all_named_group_patterns_parse() {
        for &(name, smarts) in NAMED_GROUPS {
            assert!(
                parse_smarts(smarts).is_ok(),
                "named_groups pattern for '{name}' failed to parse: {smarts}"
            );
        }
    }

    #[test]
    fn test_hydroxyl_ethanol() {
        let mol = chematic_smiles::parse("CCO").unwrap();
        let groups = detect_named_functional_groups(&mol);
        let names: Vec<&str> = groups.iter().map(|g| g.name).collect();
        assert!(
            names.contains(&"hydroxyl"),
            "expected hydroxyl in {names:?}"
        );
    }

    #[test]
    fn test_carbonyl_acetone() {
        let mol = chematic_smiles::parse("CC(=O)C").unwrap();
        let groups = detect_named_functional_groups(&mol);
        let names: Vec<&str> = groups.iter().map(|g| g.name).collect();
        assert!(
            names.contains(&"carbonyl"),
            "expected carbonyl in {names:?}"
        );
        assert!(names.contains(&"ketone"), "expected ketone in {names:?}");
        assert!(
            !names.contains(&"aldehyde"),
            "aldehyde should NOT match acetone: {names:?}"
        );
        assert!(
            !names.contains(&"hydroxyl"),
            "hydroxyl should NOT match acetone: {names:?}"
        );
    }

    #[test]
    fn test_carboxyl_acetic_acid() {
        let mol = chematic_smiles::parse("CC(=O)O").unwrap();
        let groups = detect_named_functional_groups(&mol);
        let names: Vec<&str> = groups.iter().map(|g| g.name).collect();
        assert!(
            names.contains(&"carboxyl"),
            "expected carboxyl in {names:?}"
        );
        assert!(
            names.contains(&"hydroxyl"),
            "expected hydroxyl in {names:?}"
        );
        assert!(
            !names.contains(&"amine"),
            "amine should NOT match acetic acid: {names:?}"
        );
    }

    #[test]
    fn test_amine_methylamine() {
        let mol = chematic_smiles::parse("CN").unwrap();
        let groups = detect_named_functional_groups(&mol);
        let names: Vec<&str> = groups.iter().map(|g| g.name).collect();
        assert!(
            names.contains(&"primary_amine"),
            "expected primary_amine in {names:?}"
        );
    }

    #[test]
    fn test_halogen_chlorobenzene() {
        let mol = chematic_smiles::parse("c1ccccc1Cl").unwrap();
        let groups = detect_named_functional_groups(&mol);
        let names: Vec<&str> = groups.iter().map(|g| g.name).collect();
        assert!(names.contains(&"halogen"), "expected halogen in {names:?}");
        assert!(
            names.contains(&"aromatic"),
            "expected aromatic in {names:?}"
        );
    }

    #[test]
    fn test_nitrile_acetonitrile() {
        let mol = chematic_smiles::parse("CC#N").unwrap();
        let groups = detect_named_functional_groups(&mol);
        let names: Vec<&str> = groups.iter().map(|g| g.name).collect();
        assert!(names.contains(&"nitrile"), "expected nitrile in {names:?}");
    }
}
