//! Chemical formula parsing: convert Hill-notation strings to element count maps.
//!
//! Mirrors the API of PyPI libraries **chemparse** and **chemformula**, which provide
//! `parse_formula("C6H12O6")` → `{"C": 6, "H": 12, "O": 6}` style conversion.
//!
//! # Supported syntax
//!
//! - One- and two-letter element symbols (case-sensitive: `Ca`, `Fe`, `Cl`, …)
//! - Atom counts (omitted = 1): `H2O`, `CH4`, `C6H12O6`
//! - Parentheses with multipliers: `Ca(OH)2` → `{"Ca":1,"O":2,"H":2}`
//! - Square brackets (SMILES-style): `[NH4]+` → `{"N":1,"H":4}`
//! - Trailing charge notation (`+`, `-`, `+2`, `-1`) is ignored
//!
//! # Examples
//!
//! ```
//! # use chematic_chem::parse_formula;
//! let counts = parse_formula("C6H12O6").unwrap();
//! assert_eq!(counts["C"], 6);
//! assert_eq!(counts["H"], 12);
//! assert_eq!(counts["O"], 6);
//!
//! let ca_oh2 = parse_formula("Ca(OH)2").unwrap();
//! assert_eq!(ca_oh2["Ca"], 1);
//! assert_eq!(ca_oh2["O"], 2);
//! assert_eq!(ca_oh2["H"], 2);
//! ```

use std::collections::HashMap;

/// Error type for chemical formula parsing.
#[derive(Debug, Clone, PartialEq)]
pub enum FormulaParseError {
    /// The formula string is empty.
    EmptyFormula,
    /// An unknown or malformed element symbol was encountered.
    UnknownElement(String),
    /// Parentheses are mismatched.
    UnbalancedParentheses,
    /// A malformed number was encountered.
    InvalidCount(String),
}

impl std::fmt::Display for FormulaParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyFormula => write!(f, "formula string is empty"),
            Self::UnknownElement(s) => write!(f, "unknown element symbol: {s}"),
            Self::UnbalancedParentheses => write!(f, "unbalanced parentheses in formula"),
            Self::InvalidCount(s) => write!(f, "invalid atom count: {s}"),
        }
    }
}

impl std::error::Error for FormulaParseError {}

/// Parse a Hill-notation molecular formula string into an element count map.
///
/// Returns a `HashMap<String, u32>` mapping element symbol to total atom count.
/// Charge annotations (`+`, `-`, `+2`, `-1`) at the end of the formula are ignored.
///
/// # Examples
///
/// ```
/// # use chematic_chem::parse_formula;
/// assert_eq!(parse_formula("H2O").unwrap()["O"], 1);
/// assert_eq!(parse_formula("Ca(OH)2").unwrap()["H"], 2);
/// ```
pub fn parse_formula(formula: &str) -> Result<HashMap<String, u32>, FormulaParseError> {
    let formula = formula.trim();
    if formula.is_empty() {
        return Err(FormulaParseError::EmptyFormula);
    }
    // Strip square brackets and their contents for SMILES-like notation.
    // [NH4]+ → NH4+ → then strip trailing +/-
    let formula = strip_brackets(formula);
    // Strip trailing charge: +, -, +2, -3, 2+, 2-, etc.
    let formula = strip_charge(&formula);

    let mut counts: HashMap<String, u32> = HashMap::new();
    let chars: Vec<char> = formula.chars().collect();
    parse_segment(&chars, 0, &mut counts)?;
    Ok(counts)
}

/// Strip SMILES-style square brackets, keeping inner content.
/// "[NH4]" → "NH4", "[Ca]2+" → "Ca2+"
fn strip_brackets(s: &str) -> String {
    s.chars().filter(|&c| c != '[' && c != ']').collect()
}

/// Strip trailing charge sign(s) (`+` or `-`) from a formula string.
///
/// Only the sign characters themselves are stripped; preceding atom-count digits
/// are left untouched.  Charge-magnitude digits before the sign (e.g. the "2"
/// in "Fe2+") are handled by `parse_segment` which skips unknown characters.
///
/// If no sign is present the string is returned unchanged.
fn strip_charge(s: &str) -> &str {
    let bytes = s.as_bytes();
    let mut end = bytes.len();
    // Strip any trailing + or - characters (there may be more than one in exotic notations).
    while end > 0 && (bytes[end - 1] == b'+' || bytes[end - 1] == b'-') {
        end -= 1;
    }
    &s[..end]
}

/// Parse a flat (non-nested) or nested formula segment starting at `pos`.
/// Returns the position after the last parsed character.
fn parse_segment(
    chars: &[char],
    start: usize,
    counts: &mut HashMap<String, u32>,
) -> Result<usize, FormulaParseError> {
    let mut i = start;
    while i < chars.len() {
        let c = chars[i];
        if c == '(' {
            // Parse sub-group recursively.
            let mut sub_counts: HashMap<String, u32> = HashMap::new();
            let end = parse_segment(chars, i + 1, &mut sub_counts)?;
            // After the closing ')' there may be a multiplier.
            let (mult, next) = parse_count(chars, end + 1)?; // end points to ')'
            for (elem, cnt) in sub_counts {
                let total = cnt
                    .checked_mul(mult)
                    .ok_or_else(|| FormulaParseError::InvalidCount(mult.to_string()))?;
                let entry = counts.entry(elem).or_insert(0);
                *entry = entry
                    .checked_add(total)
                    .ok_or_else(|| FormulaParseError::InvalidCount(total.to_string()))?;
            }
            i = next;
        } else if c == ')' {
            // Signal end of sub-group to caller.
            return Ok(i);
        } else if c.is_uppercase() {
            // Start of an element symbol.
            let mut sym = String::new();
            sym.push(c);
            i += 1;
            // Consume following lowercase letters (e.g. "Ca", "Fe", "Cl").
            while i < chars.len() && chars[i].is_lowercase() {
                sym.push(chars[i]);
                i += 1;
            }
            // Validate: symbol must be non-empty (always true here).
            let (cnt, next) = parse_count(chars, i)?;
            let entry = counts.entry(sym).or_insert(0);
            *entry = entry
                .checked_add(cnt)
                .ok_or_else(|| FormulaParseError::InvalidCount(cnt.to_string()))?;
            i = next;
        } else if c.is_ascii_digit() {
            // Stray digit (shouldn't appear at top level after stripping charge).
            i += 1;
        } else if c == '+' || c == '-' {
            // Charge character (already stripped from end, but may appear mid-string for some notations).
            i += 1;
        } else {
            return Err(FormulaParseError::UnknownElement(c.to_string()));
        }
    }
    Ok(i)
}

/// Parse an optional integer count starting at `pos`.
/// Returns `(count, new_pos)`. If no digit, returns `(1, pos)`.
fn parse_count(chars: &[char], pos: usize) -> Result<(u32, usize), FormulaParseError> {
    let mut i = pos;
    let mut n_str = String::new();
    while i < chars.len() && chars[i].is_ascii_digit() {
        n_str.push(chars[i]);
        i += 1;
    }
    if n_str.is_empty() {
        Ok((1, i))
    } else {
        n_str
            .parse()
            .map(|count| (count, i))
            .map_err(|_| FormulaParseError::InvalidCount(n_str))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counts(formula: &str) -> HashMap<String, u32> {
        parse_formula(formula).unwrap_or_else(|e| panic!("parse {formula:?}: {e}"))
    }

    #[test]
    fn test_empty_formula_error() {
        assert!(matches!(
            parse_formula(""),
            Err(FormulaParseError::EmptyFormula)
        ));
        assert!(matches!(
            parse_formula("   "),
            Err(FormulaParseError::EmptyFormula)
        ));
    }

    #[test]
    fn test_count_overflow_is_rejected_without_panicking() {
        assert!(matches!(
            parse_formula("C999999999999999999999999"),
            Err(FormulaParseError::InvalidCount(_))
        ));
    }

    #[test]
    fn test_water() {
        let m = counts("H2O");
        assert_eq!(m["H"], 2);
        assert_eq!(m["O"], 1);
    }

    #[test]
    fn test_glucose() {
        let m = counts("C6H12O6");
        assert_eq!(m["C"], 6);
        assert_eq!(m["H"], 12);
        assert_eq!(m["O"], 6);
    }

    #[test]
    fn test_calcium_hydroxide() {
        let m = counts("Ca(OH)2");
        assert_eq!(m["Ca"], 1);
        assert_eq!(m["O"], 2);
        assert_eq!(m["H"], 2);
    }

    #[test]
    fn test_ammonium_ion_bracket() {
        // SMILES-style bracket atom with charge
        let m = counts("[NH4]+");
        assert_eq!(m["N"], 1);
        assert_eq!(m["H"], 4);
        // Charge should be ignored
        assert!(!m.contains_key("+"));
    }

    #[test]
    fn test_sulfuric_acid() {
        let m = counts("H2SO4");
        assert_eq!(m["H"], 2);
        assert_eq!(m["S"], 1);
        assert_eq!(m["O"], 4);
    }

    #[test]
    fn test_ethanol() {
        let m = counts("C2H5OH");
        assert_eq!(m["C"], 2);
        assert_eq!(m["H"], 6); // 5 + 1 (OH)
        assert_eq!(m["O"], 1);
    }

    #[test]
    fn test_iron_iii_chloride() {
        let m = counts("FeCl3");
        assert_eq!(m["Fe"], 1);
        assert_eq!(m["Cl"], 3);
    }

    #[test]
    fn test_nested_parens() {
        // Al2(SO4)3
        let m = counts("Al2(SO4)3");
        assert_eq!(m["Al"], 2);
        assert_eq!(m["S"], 3);
        assert_eq!(m["O"], 12);
    }

    #[test]
    fn test_trailing_charge_sign_stripped() {
        // Simple sign-only charge (charge magnitude in sign position, not as a digit).
        // "NH4+" → strip "+" → "NH4" → {N:1, H:4}
        let m = counts("NH4+");
        assert_eq!(m["N"], 1);
        assert_eq!(m["H"], 4);

        // Anion notation with just "-"
        let m2 = counts("Cl-");
        assert_eq!(m2["Cl"], 1);
    }

    #[test]
    fn test_single_atom() {
        let m = counts("Fe");
        assert_eq!(m["Fe"], 1);
        assert_eq!(m.len(), 1);
    }

    #[test]
    fn test_methane() {
        let m = counts("CH4");
        assert_eq!(m["C"], 1);
        assert_eq!(m["H"], 4);
    }

    #[test]
    fn test_benzene() {
        let m = counts("C6H6");
        assert_eq!(m["C"], 6);
        assert_eq!(m["H"], 6);
    }

    #[test]
    fn test_formula_roundtrip_consistency() {
        // Verify parse_formula is consistent with calc_mol_formula output.
        use chematic_smiles::parse;
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap(); // aspirin
        let formula = crate::calc_mol_formula(&mol);
        let parsed = counts(&formula);
        assert_eq!(parsed["C"], 9);
        assert_eq!(parsed["H"], 8);
        assert_eq!(parsed["O"], 4);
    }
}
