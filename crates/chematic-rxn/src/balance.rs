//! Reaction atom-balance checking.
//!
//! [`balance_check`] counts the atoms on each side of a reaction and reports
//! whether they are equal.  Implicit hydrogens are included so that
//! `C>>C` (methane on both sides) counts as balanced.

use std::collections::BTreeMap;

use chematic_core::{Element, Molecule, implicit_hcount};

use crate::reaction::Reaction;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Result of a reaction balance check.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BalanceResult {
    /// `true` when every element has equal counts on both sides.
    pub balanced: bool,
    /// Element counts for all reactants combined (including implicit H).
    pub reactant_formula: BTreeMap<u8, u32>,
    /// Element counts for all products combined (including implicit H).
    pub product_formula: BTreeMap<u8, u32>,
}

impl BalanceResult {
    /// Return a human-readable difference string, e.g. `"C: 1 reactant vs 2 product"`.
    pub fn diff(&self) -> Vec<String> {
        let mut keys: std::collections::BTreeSet<u8> = std::collections::BTreeSet::new();
        keys.extend(self.reactant_formula.keys());
        keys.extend(self.product_formula.keys());

        let mut out = Vec::new();
        for an in keys {
            let r = self.reactant_formula.get(&an).copied().unwrap_or(0);
            let p = self.product_formula.get(&an).copied().unwrap_or(0);
            if r != p {
                let sym = Element::from_atomic_number(an)
                    .map(|e| e.symbol().to_string())
                    .unwrap_or_else(|| format!("#{an}"));
                out.push(format!("{sym}: {r} reactant vs {p} product"));
            }
        }
        out
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Check whether a reaction is atom-balanced.
///
/// Counts every explicit atom plus its implicit hydrogens in the reactants
/// and products, then compares element-by-element.
///
/// Agents are ignored (they are catalysts / solvents, not consumed or formed).
pub fn balance_check(rxn: &Reaction) -> BalanceResult {
    let reactant_formula = sum_formula(&rxn.reactants);
    let product_formula = sum_formula(&rxn.products);
    let balanced = reactant_formula == product_formula;
    BalanceResult {
        balanced,
        reactant_formula,
        product_formula,
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn sum_formula(mols: &[Molecule]) -> BTreeMap<u8, u32> {
    let mut counts: BTreeMap<u8, u32> = BTreeMap::new();
    for mol in mols {
        for (idx, atom) in mol.atoms() {
            let an = atom.element.atomic_number();
            *counts.entry(an).or_insert(0) += 1;
            // Add implicit hydrogens (atomic number 1).
            let h = implicit_hcount(mol, idx) as u32;
            if h > 0 {
                *counts.entry(1).or_insert(0) += h;
            }
        }
    }
    counts
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reaction::parse_reaction;

    #[test]
    fn test_simple_balanced() {
        // CH4 + O2 → CO2 + H2O  (not stoichiometrically balanced but atom-balanced)
        // Use a simpler example: ethane → ethane (trivially balanced)
        let rxn = parse_reaction("CC>>CC").unwrap();
        let result = balance_check(&rxn);
        assert!(result.balanced, "CC>>CC should be balanced");
    }

    #[test]
    fn test_unbalanced() {
        // Methane → ethane: 1C vs 2C — unbalanced
        let rxn = parse_reaction("C>>CC").unwrap();
        let result = balance_check(&rxn);
        assert!(!result.balanced, "C>>CC should be unbalanced");
    }

    #[test]
    fn test_balanced_with_h() {
        // Methanol dehydration: 2 CH3OH → CH3OCH3 + H2O
        // Both sides: 2C 8H 2O → balanced
        let rxn = parse_reaction("CO.CO>>COC.O").unwrap();
        let result = balance_check(&rxn);
        assert!(
            result.balanced,
            "2 methanol → dimethyl ether + water should be balanced"
        );
    }

    #[test]
    fn test_unbalanced_diff() {
        let rxn = parse_reaction("C>>CC").unwrap();
        let result = balance_check(&rxn);
        let diff = result.diff();
        assert!(!diff.is_empty(), "diff should report imbalance");
        assert!(
            diff.iter().any(|s| s.contains('C')),
            "diff should mention C"
        );
    }

    #[test]
    fn test_agents_ignored() {
        // Agent (middle section) should not contribute to balance count
        // CC>O>CC — ethane with water as agent
        let rxn = parse_reaction("CC>O>CC").unwrap();
        let result = balance_check(&rxn);
        assert!(result.balanced, "agents should be ignored in balance check");
    }
}
