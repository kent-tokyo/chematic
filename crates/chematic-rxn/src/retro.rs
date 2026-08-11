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
            Self::Ester => "Ester",
            Self::Ether => "Ether",
            Self::CNBond => "CNBond",
            Self::CCBond => "CCBond",
            Self::CSBond => "CSBond",
            Self::Other => "Other",
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
        // Mitsunobu ether: inverted secondary alcohol + phenol/alcohol.
        // #296: `CX4` (sp3-carbon connectivity query) is SMARTS syntax this
        // crate's SMILES-based template parser cannot express (see
        // `mol_to_query`'s doc) -- broadened to any non-aromatic carbon.
        // Trades "restricted to sp3 C" for "actually parses and fires";
        // documented, not silently dropped.
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
        // PMB (p-methoxybenzyl) ether. #296: the original SMIRKS was
        // malformed (unbalanced ring-closure digit; `OCH2` used as if it
        // were a single valid bracket-atom symbol, which it isn't) --
        // rewritten as a real aryl-O-CH2-aryl(p-OMe) fragment, atom-map :2
        // on the benzylic CH2 carbon that carries through to the PMB-Br
        // product, matching every other ether template's shape.
        name: "pmb_ether",
        smirks: "[c:1]O[CH2:2]c1ccc(OC)cc1>>[c:1][OH].Br[CH2:2]c1ccc(OC)cc1",
        reaction_class: RetroClass::Ether,
    },
    // ── C–N bond ─────────────────────────────────────────────────────────────
    RetroTemplate {
        // #296: `CX4`/`NX3` connectivity queries aren't expressible in this
        // crate's SMILES-based template grammar -- broadened to any
        // non-aromatic C/N. See mitsunobu_ether's comment for the same
        // tradeoff.
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
        // #296: `CX4` broadened to any non-aromatic C (see mitsunobu_ether).
        name: "n_alkylation",
        smirks: "[C:1][N:2]>>[C:1]Br.[NH:2]",
        reaction_class: RetroClass::CNBond,
    },
    RetroTemplate {
        // Mitsunobu N: inverted secondary alcohol + amine.
        // #296: `CX4` broadened to any non-aromatic C (see mitsunobu_ether).
        name: "mitsunobu_n",
        smirks: "[NH:1][C:2]>>[N:1].[OH][C:2]",
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
        // Explicit `-` (not the default implicit-aromatic bond two adjacent
        // aromatic atoms get in this crate's SMILES-based template parser) is
        // load-bearing: a real biaryl axis is *never* an aromatic bond (the two
        // rings' own aromatic systems don't extend across it), only ordinary
        // intra-ring aromatic bonds are. Without the `-`, `[c:1][c:2]` matched
        // every intra-ring aromatic C-C bond and *zero* real biaryl bonds --
        // the opposite of what was intended (issue #294).
        name: "suzuki_biaryl",
        smirks: "[c:1]-[c:2]>>[c:1]Br.[c:2]B(O)O",
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
        // Negishi coupling — aryl + alkyl zinc.
        // #296: `CX4` broadened to any non-aromatic C (see mitsunobu_ether)
        // -- `[c:1][C:2]` still can't match a biaryl bond (that needs two
        // aromatic atoms, `[c:1][c:2]`, covered separately by suzuki_biaryl).
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
        // Alpha carbon must bear at least one H (enolisable position).
        // #296: `H1,H2,H3` (an H-count OR-list) isn't expressible -- a
        // bracket atom's H-count is a single exact value, not a set, in
        // this crate's SMILES-based template grammar. Narrowed to the
        // smallest listed alternative (`H1`) rather than dropped entirely
        // (dropping it would also match non-enolisable H0 alpha carbons,
        // i.e. quaternary centers, which can never actually do an aldol --
        // a real chemistry error, not just reduced coverage). Known,
        // disclosed tradeoff: real CH2/CH3 alpha-carbon aldol precursors
        // are now missed (false negatives), never falsely matched (no
        // false positives).
        name: "aldol",
        smirks: "[CH1:1][C:2](=[O:3])>>[C:1].[C:2](=[O:3])",
        reaction_class: RetroClass::CCBond,
    },
    RetroTemplate {
        // Michael addition: beta carbon must be sp3 with at least one H.
        // #296: same `CX4` + H-count-OR-list issue as aldol above -- `CX4`
        // dropped (see mitsunobu_ether), H-count narrowed to `H1` (see
        // aldol's comment for the false-negative-over-false-positive
        // reasoning).
        name: "michael_addition",
        smirks: "[CH1:1][C:2][C:3](=[O:4])>>[C:1].[C:2]=[C:3]",
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
        // Friedel-Crafts alkylation — restrict to benzylic CH2 to avoid
        // matching every aryl–alkyl bond in the molecule.
        // #296: the `;X4` (AND-combinator + connectivity query) suffix
        // isn't expressible and is dropped -- `CH2` alone (an exact H
        // -count, valid bracket-atom syntax) already does almost all of
        // the intended restricting work here, `X4` was near-redundant
        // given `CH2` on its own.
        name: "friedel_crafts_alkyl",
        smirks: "[c:1][CH2:2]>>[c:1].[CH2:2]Cl",
        reaction_class: RetroClass::CCBond,
    },
    RetroTemplate {
        // Mannich: alpha-aminomethylation — alpha CH2 with enolisable H.
        // #296: `CX4` broadened to any non-aromatic C (see mitsunobu_ether).
        name: "mannich",
        smirks: "[CH2:1][C:2][N:3]>>[C:1]=O.[C:2].[N:3]",
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
        // #296: a bare, unbracketed `H` is not valid SMILES -- `H` is not
        // part of the organic subset this crate's SMILES parser accepts
        // outside brackets (only B/C/N/O/P/S/F/Cl/Br/I), matching the
        // SMILES spec (explicit terminal hydrogen must always be written
        // `[H]`). Not a parser bug; the template's own SMIRKS was invalid
        // SMILES from the start.
        name: "aryl_halide_oxidative_add",
        smirks: "[c:1]Br>>[c:1][H]",
        reaction_class: RetroClass::CSBond,
    },
    RetroTemplate {
        // #296: trailing bare `H` -> `[H]`, same issue as
        // aryl_halide_oxidative_add above.
        name: "phosphonate",
        smirks: "[C:1][P:2](=[O:3])([O:4])([O:5])>>[C:1]Br.[P:2](=[O:3])([O:4])([O:5])[H]",
        reaction_class: RetroClass::CSBond,
    },
    RetroTemplate {
        // #296: bare `H` -> `[H]`, same issue as aryl_halide_oxidative_add.
        name: "sp3_ch_bromination",
        smirks: "[C:1]Br>>[C:1][H]",
        reaction_class: RetroClass::CSBond,
    },
    RetroTemplate {
        name: "nitrile_from_halide",
        smirks: "[C:1][C:2]#[N]>>[C:1][C:2]Br",
        reaction_class: RetroClass::CSBond,
    },
    RetroTemplate {
        // #296: `[CF3]` was never valid bracket-atom syntax (a bracket atom
        // is always exactly one atom; a trifluoromethyl group is four) --
        // rewritten as the real 4-atom substructure `C(F)(F)F`.
        name: "trifluoromethyl",
        smirks: "[C:1]C(F)(F)F>>[C:1]I",
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
            let smiles: Vec<String> = precursor_set.iter().map(canonical_smiles).collect();

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
        assert!(
            !results.is_empty(),
            "should find disconnections in acetanilide"
        );

        let amide_hits: Vec<_> = results
            .iter()
            .filter(|r| r.template_name.starts_with("amide"))
            .collect();
        assert!(
            !amide_hits.is_empty(),
            "at least one amide template should match"
        );

        // Check that both precursors are present
        let all_smiles: Vec<&str> = results
            .iter()
            .flat_map(|r| r.precursor_smiles.iter().map(|s| s.as_str()))
            .collect();
        // At least one result should contain an acid or amine fragment
        assert!(
            all_smiles
                .iter()
                .any(|s| s.contains("C(=O)O") || s.contains("N")),
            "precursors should include acid or amine fragments"
        );
    }

    #[test]
    fn test_retro_ester() {
        // methyl acetate: ester template should give acetic acid + methanol
        let m = mol("CC(=O)OC");
        let results = retro_disconnect(&m, DEFAULT_TEMPLATES, 0);
        let ester_hits: Vec<_> = results
            .iter()
            .filter(|r| r.reaction_class == RetroClass::Ester)
            .collect();
        assert!(
            !ester_hits.is_empty(),
            "ester template should match methyl acetate"
        );
    }

    #[test]
    fn test_retro_no_match() {
        // benzene has no breakable bonds for any template
        let m = mol("c1ccccc1");
        let results = retro_disconnect(&m, DEFAULT_TEMPLATES, 0);
        assert!(
            results.is_empty(),
            "plain benzene has no real disconnection, got: {:?}",
            results.iter().map(|r| &r.template_name).collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_suzuki_biaryl_does_not_match_intra_ring_bonds() {
        // Issue #294: `[c:1][c:2]` with no ring-crossing constraint matched
        // any aromatic C-C bond, including ordinary intra-ring ones.
        for smiles in ["c1ccccc1", "c1ccc2ccccc2c1", "c1ccncc1"] {
            let m = mol(smiles);
            let results = retro_disconnect(&m, DEFAULT_TEMPLATES, 50);
            assert!(
                !results.iter().any(|r| r.template_name == "suzuki_biaryl"),
                "{smiles}: suzuki_biaryl must not fire on a molecule with no biaryl bond"
            );
        }

        // Two rings connected only by an ether oxygen -- no direct ring-to-ring
        // C-C bond exists, so suzuki_biaryl must not fire (aryl_ether_ullmann
        // is the correct match here).
        let diphenyl_ether = mol("c1ccccc1Oc1ccccc1");
        let results = retro_disconnect(&diphenyl_ether, DEFAULT_TEMPLATES, 50);
        assert!(
            !results.iter().any(|r| r.template_name == "suzuki_biaryl"),
            "diphenyl ether has no biaryl C-C bond"
        );
    }

    #[test]
    fn test_suzuki_biaryl_matches_real_biaryl_bond() {
        // biphenyl: the two rings ARE connected by a real (non-ring) C-C bond.
        let biphenyl = mol("c1ccc(-c2ccccc2)cc1");
        let results = retro_disconnect(&biphenyl, DEFAULT_TEMPLATES, 50);
        assert!(
            results.iter().any(|r| r.template_name == "suzuki_biaryl"),
            "suzuki_biaryl must still fire on a genuine biaryl bond"
        );
    }

    #[test]
    fn test_retro_max_results() {
        let m = mol("CC(=O)Nc1ccc(S(N)(=O)=O)cc1"); // sulfanilamide
        let results = retro_disconnect(&m, DEFAULT_TEMPLATES, 3);
        assert!(results.len() <= 3, "max_results=3 should be respected");
    }

    #[test]
    fn test_retro_deduplication() {
        // Two templates may produce the same precursor set — check dedup works
        let m = mol("CC(=O)Nc1ccccc1");
        let results = retro_disconnect(&m, DEFAULT_TEMPLATES, 0);
        let mut keys: Vec<String> = results
            .iter()
            .map(|r| {
                let mut s = r.precursor_smiles.clone();
                s.sort();
                s.join(".")
            })
            .collect();
        let before = keys.len();
        keys.dedup();
        assert_eq!(before, keys.len(), "no duplicate precursor sets");
    }

    #[test]
    fn test_retro_ether() {
        // anisole: aryl_ether_snar or aryl_ether_ullmann should match c-O-C
        let m = mol("COc1ccccc1");
        let results = retro_disconnect(&m, DEFAULT_TEMPLATES, 0);
        let ether_hits: Vec<_> = results
            .iter()
            .filter(|r| r.reaction_class == RetroClass::Ether)
            .collect();
        assert!(
            !ether_hits.is_empty(),
            "ether templates should match anisole"
        );
    }

    #[test]
    fn test_retro_class_filter() {
        // filter to only amide templates
        let collected: Vec<RetroTemplate> = DEFAULT_TEMPLATES
            .iter()
            .filter(|t| t.reaction_class == RetroClass::AmideBond)
            .map(|t| RetroTemplate {
                name: t.name,
                smirks: t.smirks,
                reaction_class: t.reaction_class,
            })
            .collect();
        let m = mol("CC(=O)Nc1ccccc1");
        let results = retro_disconnect(&m, &collected, 0);
        assert!(
            results
                .iter()
                .all(|r| r.reaction_class == RetroClass::AmideBond)
        );
    }

    #[test]
    fn test_default_template_count() {
        assert!(
            DEFAULT_TEMPLATES.len() >= 50,
            "library should have at least 50 templates"
        );
    }

    /// Results from a specific named template only -- reduces boilerplate
    /// in the #296 representative-correctness tests below, which each care
    /// about one template firing (or not) on a specific probe molecule,
    /// not the full cross-template `retro_disconnect` result set.
    fn hits_for(m: &Molecule, template_name: &str) -> Vec<RetroResult> {
        retro_disconnect(m, DEFAULT_TEMPLATES, 0)
            .into_iter()
            .filter(|r| r.template_name == template_name)
            .collect()
    }

    // ── #296: representative positive/negative/false-positive tests for
    // the 14 templates whose SMIRKS previously failed to parse (and so
    // never matched anything, for any caller, ever -- see #296). Parsing
    // successfully and matching chemically-sensibly are different claims;
    // these are checked separately from `all_default_templates_parse`
    // above, which only proves the SMIRKS parses.

    #[test]
    fn mitsunobu_ether_matches_dialkyl_ether() {
        // methyl propyl ether: retro should give methanol + propanol.
        let m = mol("COCCC");
        let hits = hits_for(&m, "mitsunobu_ether");
        assert!(!hits.is_empty(), "should disconnect an aliphatic ether");
        let all_smiles: Vec<&str> = hits
            .iter()
            .flat_map(|r| r.precursor_smiles.iter().map(String::as_str))
            .collect();
        assert!(
            all_smiles.iter().any(|s| s.contains('O')),
            "precursors should be alcohols"
        );
    }

    #[test]
    fn mitsunobu_ether_does_not_fire_on_carbonyl_oxygen() {
        // acetone: the carbonyl O has no single O-C bond in the position
        // this template's [O:1][C:2] pattern needs (it's a C=O double
        // bond, not the ether's single bond) -- confirms broadening CX4
        // away didn't also start matching unrelated C=O oxygens.
        let m = mol("CC(=O)C");
        let hits = hits_for(&m, "mitsunobu_ether");
        assert!(
            hits.is_empty(),
            "should not treat a ketone's C=O as an ether O-C bond"
        );
    }

    #[test]
    fn reductive_amination_matches_secondary_amine() {
        // N-ethylpropan-1-amine: retro should give an aldehyde/imine-side
        // and an amine fragment.
        let m = mol("CCNCCC");
        let hits = hits_for(&m, "reductive_amination");
        assert!(
            !hits.is_empty(),
            "should disconnect a dialkylamine C-N bond"
        );
    }

    #[test]
    fn n_alkylation_matches_alkyl_amine() {
        let m = mol("CCCNCC"); // N-ethylpropan-1-amine, same substrate, different template's own filter
        let hits = hits_for(&m, "n_alkylation");
        assert!(!hits.is_empty(), "should disconnect an N-alkyl bond");
    }

    #[test]
    fn mitsunobu_n_matches_primary_amine_alkylation() {
        let m = mol("CCCNCC");
        let hits = hits_for(&m, "mitsunobu_n");
        assert!(
            !hits.is_empty(),
            "should disconnect a secondary-amine C-N bond"
        );
    }

    #[test]
    fn negishi_matches_aryl_alkyl_bond() {
        // propylbenzene: aryl-C(sp3) bond should disconnect via negishi.
        let m = mol("c1ccccc1CCC");
        let hits = hits_for(&m, "negishi");
        assert!(!hits.is_empty(), "should disconnect an aryl-alkyl C-C bond");
    }

    #[test]
    fn negishi_does_not_fire_on_biaryl_bond() {
        // biphenyl: the aryl-aryl bond must stay suzuki_biaryl's territory
        // (both atoms aromatic), not negishi's ([c:1][C:2] needs C:2
        // non-aromatic).
        let m = mol("c1ccc(-c2ccccc2)cc1");
        let hits = hits_for(&m, "negishi");
        assert!(
            hits.is_empty(),
            "negishi should not match an aryl-aryl bond"
        );
    }

    #[test]
    fn aldol_matches_ch1_alpha_carbon() {
        // 3-methylpentan-2-one: the alpha carbon bonded to the carbonyl
        // (CH(CH3)(CH2CH3)) has exactly 1 H -- matches the template's
        // narrowed [CH1:1] (see aldol's own comment for why H1, not H2).
        let m = mol("CC(=O)C(C)CC");
        let hits = hits_for(&m, "aldol");
        assert!(
            !hits.is_empty(),
            "should disconnect an enolisable alpha-CH1/C=O bond"
        );
    }

    #[test]
    fn aldol_does_not_fire_on_quaternary_alpha_carbon() {
        // 2,2-dimethylpropan-1-one-like center: alpha carbon has zero H
        // (fully substituted) -- chemically cannot enolise, must not match
        // even though the narrowed [CH1:1] test is a stand-in for the
        // original H1,H2,H3 OR-list.
        let m = mol("CC(=O)C(C)(C)C"); // methyl tert-butyl ketone: alpha C (the quaternary one) has 0 H
        let hits = hits_for(&m, "aldol");
        // the *other* alpha carbon (the methyl, CH3) also isn't CH1, so
        // this substrate has no CH1 alpha carbon at all -- correctly no match.
        assert!(
            hits.is_empty(),
            "should not match when no alpha carbon has exactly 1 H"
        );
    }

    #[test]
    fn michael_addition_matches_beta_ch_carbon() {
        // 4-methylpentan-2-one-like chain with a CH beta to carbonyl via an
        // intervening carbon: build a simple enone-derived saturated
        // analog with a CH1 beta carbon two bonds from the carbonyl.
        let m = mol("CC(C)CC(=O)C");
        let hits = hits_for(&m, "michael_addition");
        assert!(
            !hits.is_empty(),
            "should disconnect a beta-CH/alpha-C/carbonyl chain"
        );
    }

    #[test]
    fn friedel_crafts_alkyl_matches_benzylic_ch2() {
        let m = mol("c1ccccc1CCl"); // benzyl chloride's own C-C precursor shape: toluene-like CH2 substrate
        let hits = hits_for(&m, "friedel_crafts_alkyl");
        assert!(!hits.is_empty(), "should disconnect a benzylic CH2");
    }

    #[test]
    fn mannich_matches_alpha_aminomethyl_chain() {
        let m = mol("CC(N(C)C)CC"); // aminomethyl-flanked sp3 chain
        let hits = hits_for(&m, "mannich");
        assert!(
            !hits.is_empty(),
            "should disconnect a CH2-C-N Mannich-type chain"
        );
    }

    #[test]
    fn trifluoromethyl_matches_cf3_group() {
        let m = mol("CC(F)(F)F"); // ethyl-CF3-like substrate (1,1,1-trifluoroethane)
        let hits = hits_for(&m, "trifluoromethyl");
        assert!(!hits.is_empty(), "should disconnect a C-CF3 bond");
    }

    #[test]
    fn trifluoromethyl_does_not_fire_without_three_fluorines() {
        let m = mol("CC(F)F"); // only 2 F -- not a CF3 group
        let hits = hits_for(&m, "trifluoromethyl");
        assert!(
            hits.is_empty(),
            "should require all three fluorines of a true CF3 group"
        );
    }

    #[test]
    fn pmb_ether_matches_pmb_protected_phenol() {
        // Tests pmb_ether's own SMIRKS directly (`run_reactants`, not
        // `retro_disconnect`'s DEFAULT_TEMPLATES sweep): a plain
        // PMB-protected phenol also matches the pre-existing, already
        // -parsing `benzyl_ether` template, which produces the identical
        // canonical precursor set (4-methoxybenzyl bromide + phenol) --
        // `retro_disconnect`'s deliberate cross-template dedup (see its
        // own doc comment) then hides pmb_ether's result behind
        // benzyl_ether's, which is correct, expected behavior of that
        // dedup, not a pmb_ether defect. This test isolates pmb_ether's
        // own matching behavior from that unrelated, already-tested dedup
        // feature.
        let m = mol("c1ccc(OCc2ccc(OC)cc2)cc1"); // phenol protected as its PMB ether
        let pmb_ether_smirks = DEFAULT_TEMPLATES
            .iter()
            .find(|t| t.name == "pmb_ether")
            .expect("pmb_ether template exists")
            .smirks;
        let results =
            crate::run_reactants(pmb_ether_smirks, &[&m]).expect("pmb_ether should parse and run");
        assert!(!results.is_empty(), "should disconnect a PMB ether");
    }

    #[test]
    fn aryl_halide_oxidative_add_matches_aryl_bromide() {
        let m = mol("c1ccccc1Br"); // bromobenzene
        let hits = hits_for(&m, "aryl_halide_oxidative_add");
        assert!(!hits.is_empty(), "should disconnect an aryl C-Br bond");
    }

    #[test]
    fn sp3_ch_bromination_matches_alkyl_bromide() {
        let m = mol("CCCBr"); // 1-bromopropane
        let hits = hits_for(&m, "sp3_ch_bromination");
        assert!(!hits.is_empty(), "should disconnect an sp3 C-Br bond");
    }

    #[test]
    fn phosphonate_matches_dialkyl_phosphonate() {
        let m = mol("CCP(=O)(OCC)OCC"); // diethyl ethylphosphonate
        let hits = hits_for(&m, "phosphonate");
        assert!(
            !hits.is_empty(),
            "should disconnect a C-P(=O) phosphonate bond"
        );
    }

    /// Issue #296's acceptance criterion: every `DEFAULT_TEMPLATES` entry's
    /// SMIRKS must parse successfully. `retro_disconnect`'s `Err(_) =>
    /// continue` (this file, `retro_disconnect`) silently treats a parse
    /// failure exactly like "no match" -- a template that never parses
    /// simply never fires, for any caller, forever, with no signal. This
    /// test is the CI gate: a future built-in template with unparseable
    /// SMIRKS fails the build immediately instead of silently going dark
    /// (which is exactly how the 14 templates fixed by this same commit
    /// went unnoticed since inception -- see #296 for the full audit).
    ///
    /// Deliberately tests parsing only (`parse_reaction`), not matching
    /// against a specific molecule: a template can legitimately return
    /// `Ok(vec![])` for a probe molecule it doesn't happen to match (that's
    /// correct behavior, not a bug), so a matching-based test would either
    /// need a bespoke probe molecule per template (fragile, and this is
    /// what the representative-molecule tests below are for) or risk
    /// false negatives. Parsing success/failure is unambiguous and
    /// molecule-independent, which is exactly what #296 needs gated.
    #[test]
    fn all_default_templates_parse() {
        let mut failures = Vec::new();
        for tmpl in DEFAULT_TEMPLATES {
            if let Err(e) = crate::reaction::parse_reaction(tmpl.smirks) {
                failures.push(format!("{} ({}): {e}", tmpl.name, tmpl.smirks));
            }
        }
        assert!(
            failures.is_empty(),
            "{} of {} DEFAULT_TEMPLATES entries failed to parse:\n{}",
            failures.len(),
            DEFAULT_TEMPLATES.len(),
            failures.join("\n")
        );
    }
}
