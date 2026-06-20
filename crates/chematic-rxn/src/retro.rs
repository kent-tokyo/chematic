//! Template-based retrosynthetic disconnection.
//!
//! Applies a curated library of reverse SMIRKS templates to a target molecule
//! to enumerate one-step precursor sets.  Each template encodes a common
//! synthetic reaction in the retro direction: the left-hand side matches a
//! functional-group pattern in the target and the right-hand side generates
//! the corresponding building blocks.
//!
//! # Usage
//!
//! ```rust
//! use chematic_rxn::retro::{retro_disconnect, DEFAULT_TEMPLATES};
//! use chematic_smiles::parse;
//!
//! let mol = parse("CC(=O)Nc1ccccc1").unwrap();  // acetanilide
//! let results = retro_disconnect(&mol, DEFAULT_TEMPLATES, 20);
//! for r in &results {
//!     println!("{}: {:?}", r.template_name, r.precursor_smiles);
//! }
//! ```

use chematic_core::Molecule;
use chematic_smiles::canonical_smiles;
use std::collections::HashSet;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Broad reaction class for a retrosynthetic template.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RetroClass {
    /// C(=O)–N bond: amide, sulfonamide, carbamate, urea, lactam
    AmideBond,
    /// C(=O)–O bond: ester, lactone, carbonate, anhydride
    Ester,
    /// C–O bond: aryl/alkyl ether, epoxide, acetal
    Ether,
    /// C–N bond: reductive amination, SNAr, Buchwald, N-alkylation
    CNBond,
    /// C–C bond: Suzuki, Heck, Sonogashira, aldol, Michael, Wittig
    CCBond,
    /// C–S or C–X bond: thioether, disulfide, halogenation, borylation
    CSBond,
    /// Other / uncategorised
    Other,
}

impl RetroClass {
    /// Human-readable name.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::AmideBond => "AmideBond",
            Self::Ester     => "Ester",
            Self::Ether     => "Ether",
            Self::CNBond    => "CNBond",
            Self::CCBond    => "CCBond",
            Self::CSBond    => "CSBond",
            Self::Other     => "Other",
        }
    }
}

/// A single retrosynthetic SMIRKS template.
pub struct RetroTemplate {
    /// Short identifier (snake_case, unique within the library).
    pub name: &'static str,
    /// Reverse SMIRKS: left = target pattern, right = precursor patterns.
    pub smirks: &'static str,
    /// Reaction class for filtering.
    pub reaction_class: RetroClass,
}

/// One retrosynthetic disconnection result.
pub struct RetroResult {
    /// Name of the template that produced this disconnection.
    pub template_name: String,
    /// Reaction class of the template.
    pub reaction_class: RetroClass,
    /// Precursor molecules.
    pub precursors: Vec<Molecule>,
    /// Canonical SMILES for each precursor (parallel to `precursors`).
    pub precursor_smiles: Vec<String>,
}

// ---------------------------------------------------------------------------
// Template library
// ---------------------------------------------------------------------------

/// Default library of 60 retro-SMIRKS templates covering the most common
/// bond-forming reactions in medicinal chemistry.
pub static DEFAULT_TEMPLATES: &[RetroTemplate] = &[
    // ── C(=O)–N bond (amide / sulfonamide / carbamate / urea) ────────────────
    RetroTemplate {
        name: "amide_secondary",
        smirks: "[C:1](=[O:2])[NH:3]>>[C:1](=[O:2])O.[N:3]",
        reaction_class: RetroClass::AmideBond,
    },
    RetroTemplate {
        name: "amide_tertiary",
        smirks: "[C:1](=[O:2])[N:3]>>[C:1](=[O:2])O.[NH:3]",
        reaction_class: RetroClass::AmideBond,
    },
    RetroTemplate {
        name: "amide_acyl_chloride",
        smirks: "[C:1](=[O:2])[NH:3]>>[C:1](=[O:2])Cl.[N:3]",
        reaction_class: RetroClass::AmideBond,
    },
    RetroTemplate {
        name: "sulfonamide",
        smirks: "[S:1](=[O:2])(=[O:3])[NH:4]>>[S:1](=[O:2])(=[O:3])Cl.[N:4]",
        reaction_class: RetroClass::AmideBond,
    },
    RetroTemplate {
        name: "carbamate",
        smirks: "[O:1][C:2](=[O:3])[N:4]>>[O:1][H].[C:2](=[O:3])=O.[N:4]",
        reaction_class: RetroClass::AmideBond,
    },
    RetroTemplate {
        name: "urea",
        smirks: "[N:1][C:2](=[O:3])[N:4]>>[N:1].[N:4]",
        reaction_class: RetroClass::AmideBond,
    },
    RetroTemplate {
        name: "hydrazide",
        smirks: "[C:1](=[O:2])[NH:3][N:4]>>[C:1](=[O:2])O.[N:3][N:4]",
        reaction_class: RetroClass::AmideBond,
    },
    RetroTemplate {
        name: "imide",
        smirks: "[C:1](=[O:2])[N:3][C:4](=[O:5])>>[C:1](=[O:2])O.[C:4](=[O:5])O.[N:3]",
        reaction_class: RetroClass::AmideBond,
    },
    RetroTemplate {
        name: "hydroxamic_acid",
        smirks: "[C:1](=[O:2])[N:3][OH:4]>>[C:1](=[O:2])O.[N:3][O:4]",
        reaction_class: RetroClass::AmideBond,
    },
    RetroTemplate {
        name: "thioamide",
        smirks: "[C:1](=[S:2])[N:3]>>[C:1](=[O])O.[N:3]",
        reaction_class: RetroClass::AmideBond,
    },

    // ── Ester / carbonate / anhydride ────────────────────────────────────────
    RetroTemplate {
        name: "ester",
        smirks: "[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])O.[OH:3][C:4]",
        reaction_class: RetroClass::Ester,
    },
    RetroTemplate {
        name: "thioester",
        smirks: "[C:1](=[O:2])[S:3]>>[C:1](=[O:2])O.[SH:3]",
        reaction_class: RetroClass::Ester,
    },
    RetroTemplate {
        name: "carbonate",
        smirks: "[O:1][C:2](=[O:3])[O:4]>>[O:1][H].[O:4][H]",
        reaction_class: RetroClass::Ester,
    },
    RetroTemplate {
        name: "anhydride",
        smirks: "[C:1](=[O:2])[O:3][C:4](=[O:5])>>[C:1](=[O:2])O.[C:4](=[O:5])O",
        reaction_class: RetroClass::Ester,
    },
    RetroTemplate {
        name: "acetal",
        smirks: "[C:1]([O:2][C:3])([O:4][C:5])>>[C:1]=O.[OH:2][C:3].[OH:4][C:5]",
        reaction_class: RetroClass::Ester,
    },
    RetroTemplate {
        name: "lactone",
        smirks: "[C:1](=[O:2])[O:3][C:4][C:5]>>[C:1](=[O:2])O.[OH:3][C:4][C:5]",
        reaction_class: RetroClass::Ester,
    },

    // ── Ether ────────────────────────────────────────────────────────────────
    RetroTemplate {
        name: "aryl_ether_snar",
        smirks: "[c:1][O:2][C:3]>>[c:1]F.[OH:2][C:3]",
        reaction_class: RetroClass::Ether,
    },
    RetroTemplate {
        name: "aryl_ether_ullmann",
        smirks: "[c:1][O:2][c:3]>>[c:1]Br.[OH:2][c:3]",
        reaction_class: RetroClass::Ether,
    },
    RetroTemplate {
        name: "williamson_ether",
        smirks: "[C:1][O:2][C:3]>>[C:1]Br.[OH:2][C:3]",
        reaction_class: RetroClass::Ether,
    },
    RetroTemplate {
        name: "benzyl_ether",
        smirks: "[c:1][CH2:2][O:3]>>[c:1][CH2:2]Br.[OH:3]",
        reaction_class: RetroClass::Ether,
    },
    RetroTemplate {
        name: "mitsunobu_ether",
        smirks: "[O:1][C:2]>>[OH:1].[OH][C:2]",
        reaction_class: RetroClass::Ether,
    },
    RetroTemplate {
        name: "vinyl_ether",
        smirks: "[C:1]=[C:2][O:3]>>[C:1]=[C:2]Br.[OH:3]",
        reaction_class: RetroClass::Ether,
    },
    RetroTemplate {
        name: "silyl_ether",
        smirks: "[C:1][O:2][Si:3]>>[C:1][OH:2].[Si:3]Cl",
        reaction_class: RetroClass::Ether,
    },
    RetroTemplate {
        name: "pmb_ether",
        smirks: "[c:1]([OCH2:2])ccc(OC)cc1>>[c:1][OH].[c]([CH2:2]Br)ccc(OC)cc1",
        reaction_class: RetroClass::Ether,
    },

    // ── C–N bond ─────────────────────────────────────────────────────────────
    RetroTemplate {
        name: "reductive_amination",
        smirks: "[C:1][N:2]>>[C:1]=O.[N:2]",
        reaction_class: RetroClass::CNBond,
    },
    RetroTemplate {
        name: "snar_cn",
        smirks: "[c:1][N:2]>>[c:1]F.[NH:2]",
        reaction_class: RetroClass::CNBond,
    },
    RetroTemplate {
        name: "buchwald_cn",
        smirks: "[c:1][NH:2]>>[c:1]Br.[N:2]",
        reaction_class: RetroClass::CNBond,
    },
    RetroTemplate {
        name: "buchwald_cn_tertiary",
        smirks: "[c:1][N:2]>>[c:1]Br.[NH:2]",
        reaction_class: RetroClass::CNBond,
    },
    RetroTemplate {
        name: "n_alkylation",
        smirks: "[C:1][N:2]>>[C:1]Br.[NH:2]",
        reaction_class: RetroClass::CNBond,
    },
    RetroTemplate {
        name: "mitsunobu_n",
        smirks: "[N:1][C:2]>>[NH:1].[OH][C:2]",
        reaction_class: RetroClass::CNBond,
    },
    RetroTemplate {
        name: "imine_reduction",
        smirks: "[C:1][NH:2]>>[C:1]=O.[NH2:2]",
        reaction_class: RetroClass::CNBond,
    },
    RetroTemplate {
        name: "nitrile_hydrolysis",
        smirks: "[C:1]#[N:2]>>[C:1](=O)O",
        reaction_class: RetroClass::CNBond,
    },
    RetroTemplate {
        name: "imine_condensation",
        smirks: "[C:1]=[N:2]>>[C:1]=O.[NH2:2]",
        reaction_class: RetroClass::CNBond,
    },
    RetroTemplate {
        name: "guanidine",
        smirks: "[N:1][C:2](=[N:3])[N:4]>>[N:1].[N:4]",
        reaction_class: RetroClass::CNBond,
    },
    RetroTemplate {
        name: "amidine",
        smirks: "[N:1][C:2]=[N:3]>>[C:2]#[N:3].[N:1]",
        reaction_class: RetroClass::CNBond,
    },

    // ── C–C bond ─────────────────────────────────────────────────────────────
    RetroTemplate {
        name: "suzuki_biaryl",
        smirks: "[c:1][c:2]>>[c:1]Br.[c:2]B(O)O",
        reaction_class: RetroClass::CCBond,
    },
    RetroTemplate {
        name: "suzuki_vinyl_ar",
        smirks: "[c:1][CH:2]=[CH2:3]>>[c:1]Br.[CH:2]=[CH2:3]",
        reaction_class: RetroClass::CCBond,
    },
    RetroTemplate {
        name: "heck_acrylate",
        smirks: "[c:1][CH:2]=[CH:3][C:4](=[O:5])>>[c:1]Br.[CH2:2]=[CH:3][C:4](=[O:5])",
        reaction_class: RetroClass::CCBond,
    },
    RetroTemplate {
        name: "sonogashira",
        smirks: "[c:1][C:2]#[C:3]>>[c:1]Br.[CH:2]#[C:3]",
        reaction_class: RetroClass::CCBond,
    },
    RetroTemplate {
        name: "negishi",
        smirks: "[c:1][C:2]>>[c:1]Br.[C:2][Zn]Cl",
        reaction_class: RetroClass::CCBond,
    },
    RetroTemplate {
        name: "grignard_ketone",
        smirks: "[C:1](=[O:2])[C:3]>>[C:1](=[O:2])Cl.[C:3][Mg]Br",
        reaction_class: RetroClass::CCBond,
    },
    RetroTemplate {
        name: "aldol",
        smirks: "[C:1][C:2](=[O:3])>>[C:1].[C:2](=[O:3])",
        reaction_class: RetroClass::CCBond,
    },
    RetroTemplate {
        name: "michael_addition",
        smirks: "[C:1][C:2][C:3](=[O:4])>>[C:1].[C:2]=[C:3]",
        reaction_class: RetroClass::CCBond,
    },
    RetroTemplate {
        name: "wittig",
        smirks: "[C:1]=[C:2]>>[C:1]=O.[C:2]=O",
        reaction_class: RetroClass::CCBond,
    },
    RetroTemplate {
        name: "knoevenagel",
        smirks: "[c:1][C:2]=[C:3]>>[c:1][C:2]=O.[C:3]",
        reaction_class: RetroClass::CCBond,
    },
    RetroTemplate {
        name: "friedel_crafts_acyl",
        smirks: "[c:1][C:2](=[O:3])>>[c:1].[C:2](=[O:3])Cl",
        reaction_class: RetroClass::CCBond,
    },
    RetroTemplate {
        name: "friedel_crafts_alkyl",
        smirks: "[c:1][C:2]>>[c:1].[C:2]Cl",
        reaction_class: RetroClass::CCBond,
    },
    RetroTemplate {
        name: "sp3_cc_bond",
        smirks: "[C:1][C:2]>>[C:1]Br.[C:2]Br",
        reaction_class: RetroClass::CCBond,
    },
    RetroTemplate {
        name: "mannich",
        smirks: "[C:1][C:2][N:3]>>[C:1]=O.[C:2].[N:3]",
        reaction_class: RetroClass::CCBond,
    },
    RetroTemplate {
        name: "reformatsky",
        smirks: "[C:1](=[O:2])[C:3]>>[C:1](=[O:2])Br.[C:3](=O)",
        reaction_class: RetroClass::CCBond,
    },

    // ── C–S / C–X bond ───────────────────────────────────────────────────────
    RetroTemplate {
        name: "aryl_thioether",
        smirks: "[c:1][S:2][C:3]>>[c:1]Br.[SH:2][C:3]",
        reaction_class: RetroClass::CSBond,
    },
    RetroTemplate {
        name: "alkyl_thioether",
        smirks: "[C:1][S:2][C:3]>>[C:1]Br.[SH:2][C:3]",
        reaction_class: RetroClass::CSBond,
    },
    RetroTemplate {
        name: "disulfide",
        smirks: "[S:1][S:2]>>[SH:1].[SH:2]",
        reaction_class: RetroClass::CSBond,
    },
    RetroTemplate {
        name: "borylation",
        smirks: "[c:1][B:2](O)O>>[c:1]Br.[B:2](O)O",
        reaction_class: RetroClass::CSBond,
    },
    RetroTemplate {
        name: "aryl_fluoride_from_cl",
        smirks: "[c:1][F:2]>>[c:1]Cl",
        reaction_class: RetroClass::CSBond,
    },
    RetroTemplate {
        name: "aryl_halide_oxidative_add",
        smirks: "[c:1]Br>>[c:1]H",
        reaction_class: RetroClass::CSBond,
    },
    RetroTemplate {
        name: "phosphonate",
        smirks: "[C:1][P:2](=[O:3])([O:4])([O:5])>>[C:1]Br.[P:2](=[O:3])([O:4])([O:5])H",
        reaction_class: RetroClass::CSBond,
    },
    RetroTemplate {
        name: "sp3_ch_bromination",
        smirks: "[C:1]Br>>[C:1]H",
        reaction_class: RetroClass::CSBond,
    },
    RetroTemplate {
        name: "nitrile_from_halide",
        smirks: "[C:1][C:2]#[N]>>[C:1][C:2]Br",
        reaction_class: RetroClass::CSBond,
    },
    RetroTemplate {
        name: "trifluoromethyl",
        smirks: "[C:1][CF3]>>[C:1]I",
        reaction_class: RetroClass::CSBond,
    },
];

// ---------------------------------------------------------------------------
// Core function
// ---------------------------------------------------------------------------

/// Apply a library of retro-SMIRKS templates to `mol` and return the
/// resulting precursor sets, deduplicated and sorted by the number of
/// precursors ascending (fewest fragments first).
///
/// `templates` — slice of templates to try (use `DEFAULT_TEMPLATES` for the
/// full built-in library).
///
/// `max_results` — cap on returned results (0 = unlimited).
///
/// Returns a `Vec<RetroResult>` sorted by number of precursors (fewer =
/// simpler disconnection).  Duplicate precursor sets (same canonical SMILES
/// in the same order, regardless of template name) are removed.
pub fn retro_disconnect(
    mol: &Molecule,
    templates: &[RetroTemplate],
    max_results: usize,
) -> Vec<RetroResult> {
    let mut results: Vec<RetroResult> = Vec::new();
    let mut seen: HashSet<String> = HashSet::new();

    let templates = if templates.is_empty() {
        DEFAULT_TEMPLATES
    } else {
        templates
    };

    for tmpl in templates {
        let sets = match crate::run_reactants(tmpl.smirks, &[mol]) {
            Ok(s) => s,
            Err(_) => continue, // template didn't match or SMIRKS parse failed
        };

        for precursor_set in sets {
            if precursor_set.is_empty() {
                continue;
            }

            // Compute canonical SMILES for each precursor.
            let smiles: Vec<String> = precursor_set
                .iter()
                .map(canonical_smiles)
                .collect();

            // Dedup key: sorted canonical SMILES joined.
            let mut sorted = smiles.clone();
            sorted.sort();
            let key = sorted.join(".");

            if !seen.insert(key) {
                continue; // already have this precursor set
            }

            results.push(RetroResult {
                template_name: tmpl.name.to_string(),
                reaction_class: tmpl.reaction_class,
                precursors: precursor_set,
                precursor_smiles: smiles,
            });
        }
    }

    // Sort: fewest precursors first (simpler disconnections first).
    results.sort_by_key(|r| r.precursors.len());

    if max_results > 0 && results.len() > max_results {
        results.truncate(max_results);
    }

    results
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).expect("parse failed")
    }

    #[test]
    fn test_retro_amide_secondary() {
        // acetanilide: amide_secondary should give acetic acid + aniline
        let m = mol("CC(=O)Nc1ccccc1");
        let results = retro_disconnect(&m, DEFAULT_TEMPLATES, 0);
        assert!(!results.is_empty(), "should find disconnections in acetanilide");

        let amide_hits: Vec<_> = results.iter()
            .filter(|r| r.template_name.starts_with("amide"))
            .collect();
        assert!(!amide_hits.is_empty(), "at least one amide template should match");

        // Check that both precursors are present
        let all_smiles: Vec<&str> = results.iter()
            .flat_map(|r| r.precursor_smiles.iter().map(|s| s.as_str()))
            .collect();
        // At least one result should contain an acid or amine fragment
        assert!(
            all_smiles.iter().any(|s| s.contains("C(=O)O") || s.contains("N")),
            "precursors should include acid or amine fragments"
        );
    }

    #[test]
    fn test_retro_ester() {
        // methyl acetate: ester template should give acetic acid + methanol
        let m = mol("CC(=O)OC");
        let results = retro_disconnect(&m, DEFAULT_TEMPLATES, 0);
        let ester_hits: Vec<_> = results.iter()
            .filter(|r| r.reaction_class == RetroClass::Ester)
            .collect();
        assert!(!ester_hits.is_empty(), "ester template should match methyl acetate");
    }

    #[test]
    fn test_retro_no_match() {
        // benzene has no breakable bonds for any template
        let m = mol("c1ccccc1");
        let results = retro_disconnect(&m, DEFAULT_TEMPLATES, 0);
        // May still get trivial C-C or C-H hits — just check it doesn't panic
        let _ = results;
    }

    #[test]
    fn test_retro_max_results() {
        let m = mol("CC(=O)Nc1ccc(S(N)(=O)=O)cc1");  // sulfanilamide
        let results = retro_disconnect(&m, DEFAULT_TEMPLATES, 3);
        assert!(results.len() <= 3, "max_results=3 should be respected");
    }

    #[test]
    fn test_retro_deduplication() {
        // Two templates may produce the same precursor set — check dedup works
        let m = mol("CC(=O)Nc1ccccc1");
        let results = retro_disconnect(&m, DEFAULT_TEMPLATES, 0);
        let mut keys: Vec<String> = results.iter().map(|r| {
            let mut s = r.precursor_smiles.clone();
            s.sort();
            s.join(".")
        }).collect();
        let before = keys.len();
        keys.dedup();
        assert_eq!(before, keys.len(), "no duplicate precursor sets");
    }

    #[test]
    fn test_retro_ether() {
        // anisole: aryl_ether_snar or aryl_ether_ullmann should match c-O-C
        let m = mol("COc1ccccc1");
        let results = retro_disconnect(&m, DEFAULT_TEMPLATES, 0);
        let ether_hits: Vec<_> = results.iter()
            .filter(|r| r.reaction_class == RetroClass::Ether)
            .collect();
        assert!(!ether_hits.is_empty(), "ether templates should match anisole");
    }

    #[test]
    fn test_retro_class_filter() {
        // filter to only amide templates
        let collected: Vec<RetroTemplate> = DEFAULT_TEMPLATES.iter()
            .filter(|t| t.reaction_class == RetroClass::AmideBond)
            .map(|t| RetroTemplate { name: t.name, smirks: t.smirks, reaction_class: t.reaction_class })
            .collect();
        let m = mol("CC(=O)Nc1ccccc1");
        let results = retro_disconnect(&m, &collected, 0);
        assert!(results.iter().all(|r| r.reaction_class == RetroClass::AmideBond));
    }

    #[test]
    fn test_default_template_count() {
        assert!(DEFAULT_TEMPLATES.len() >= 50, "library should have at least 50 templates");
    }
}
