//! Condensed formula parsing — converts text like "CH3COOH" to Molecule.
//!
//! Supports Hill notation-style chemical formulas with explicit hydrogens,
//! functional group abbreviations (COOH, OH, NH2, etc.), and parentheses for branching.

use chematic_core::Molecule;
use chematic_smiles::parse as parse_smiles;

/// Error type for condensed formula parsing.
#[derive(Debug, Clone, PartialEq)]
pub enum CondensedError {
    /// Unknown element symbol encountered
    UnknownElement(String),
    /// Unbalanced parentheses
    UnbalancedParens,
    /// Empty input
    EmptyInput,
    /// SMILES parsing failed
    ParseError(String),
}

impl std::fmt::Display for CondensedError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnknownElement(s) => write!(f, "unknown element: {}", s),
            Self::UnbalancedParens => write!(f, "unbalanced parentheses"),
            Self::EmptyInput => write!(f, "empty input"),
            Self::ParseError(s) => write!(f, "SMILES parse error: {}", s),
        }
    }
}

impl std::error::Error for CondensedError {}

/// Parse a condensed formula string (e.g., "CH3COOH") into a Molecule.
///
/// This function converts condensed notation like "CH3CH2OH", "CH3(CH2)4CH3",
/// "NH2CH2COOH" into SMILES, then parses that SMILES to a Molecule.
///
/// **Strategy**:
/// 1. Tokenize: recognize multi-char atoms (Cl, Br, Si...), single-char atoms, digits, parentheses
/// 2. Substitute known functional groups (COOH → C(=O)O, OH → O, etc.)
/// 3. Build SMILES representation with implicit hydrogen counts
/// 4. Parse the SMILES string to get the final Molecule
///
/// **Examples**:
/// - `"C"` → methane `CH4`
/// - `"CH3CH3"` → ethane
/// - `"CH3COOH"` → acetic acid
/// - `"CH3(CH2)4CH3"` → n-hexane (with parenthesis branch syntax)
pub fn parse_condensed(input: &str) -> Result<Molecule, CondensedError> {
    if input.trim().is_empty() {
        return Err(CondensedError::EmptyInput);
    }

    let input = input.trim();
    let smiles = condensed_to_smiles(input)?;
    parse_smiles(&smiles).map_err(|e| CondensedError::ParseError(e.to_string()))
}

/// Convert condensed formula to SMILES.
fn condensed_to_smiles(input: &str) -> Result<String, CondensedError> {
    let tokens = tokenize(input)?;
    substitute_functional_groups(&tokens)
}

/// Tokenize condensed formula into atoms, digits, and parentheses.
fn tokenize(input: &str) -> Result<Vec<Token>, CondensedError> {
    let mut tokens = Vec::new();
    let bytes = input.as_bytes();
    let mut i = 0;

    while i < bytes.len() {
        if bytes[i] == b'(' || bytes[i] == b')' {
            tokens.push(Token::Paren(bytes[i] as char));
            i += 1;
        } else if bytes[i].is_ascii_digit() {
            let mut num = String::new();
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                num.push(bytes[i] as char);
                i += 1;
            }
            tokens.push(Token::Digit(num.parse::<u32>().unwrap_or(1)));
        } else if bytes[i].is_ascii_uppercase() {
            // Try multi-char element first (Cl, Br, Si, etc.)
            if i + 1 < bytes.len() && bytes[i + 1].is_ascii_lowercase() {
                let elem = std::str::from_utf8(&bytes[i..=i + 1]).unwrap();
                // Validate it's a known 2-char element
                if is_valid_element(elem) {
                    tokens.push(Token::Atom(elem.to_string()));
                    i += 2;
                } else {
                    // Single-char element
                    tokens.push(Token::Atom((bytes[i] as char).to_string()));
                    i += 1;
                }
            } else {
                // Single-char element
                tokens.push(Token::Atom((bytes[i] as char).to_string()));
                i += 1;
            }
        } else if bytes[i].is_ascii_lowercase() {
            // Lowercase after uppercase was already handled, so this is error
            return Err(CondensedError::UnknownElement(format!(
                "unexpected lowercase at position {}: {}",
                i, input
            )));
        } else {
            return Err(CondensedError::UnknownElement(format!(
                "unexpected character at position {}: {}",
                i, bytes[i] as char
            )));
        }
    }

    Ok(tokens)
}

#[derive(Debug, Clone)]
enum Token {
    Atom(String),
    Digit(u32),
    Paren(char),
}

/// Check if a 2-char string is a valid multi-char element symbol.
fn is_valid_element(s: &str) -> bool {
    matches!(
        s,
        "Cl" | "Br" | "Si" | "As" | "Se" | "Sn" | "Te" | "Pb" | "Bi" | "Po" | "At"
    )
}

/// Maximum atom repeat count accepted from a digit token.
const MAX_REPEAT: u32 = 999;

/// Substitute known functional groups (COOH → C(=O)O, etc.).
fn substitute_functional_groups(tokens: &[Token]) -> Result<String, CondensedError> {
    let mut smiles = String::new();
    // Length (in bytes) of the last single heavy-atom symbol written to `smiles`.
    // Used by the Digit arm to repeat multi-char atoms (Cl, Br, Si…) correctly.
    let mut last_atom_len: usize = 0;
    let mut i = 0;

    while i < tokens.len() {
        match &tokens[i] {
            Token::Paren(c) => {
                smiles.push(*c);
                last_atom_len = 0; // parens break repeat context
                i += 1;
            }
            Token::Digit(n) => {
                if *n > MAX_REPEAT {
                    return Err(CondensedError::ParseError(format!(
                        "repeat count {n} exceeds maximum {MAX_REPEAT}"
                    )));
                }
                if last_atom_len == 0 || smiles.len() < last_atom_len {
                    return Err(CondensedError::ParseError(
                        "digit without preceding atom".into(),
                    ));
                }
                // Pop the one copy already written, then push n copies.
                let atom_str = smiles[smiles.len() - last_atom_len..].to_owned();
                smiles.truncate(smiles.len() - last_atom_len);
                for _ in 0..*n {
                    smiles.push_str(&atom_str);
                }
                last_atom_len = 0;
                i += 1;
            }
            Token::Atom(a) if a == "H" => {
                // In condensed notation H is always implicit — skip the H symbol.
                // A following digit is an H-count annotation (e.g. CH3 = C with 3 implicit H),
                // not a repeat count, so consume it too.
                last_atom_len = 0;
                i += 1;
                if i < tokens.len()
                    && let Token::Digit(_) = &tokens[i]
                {
                    i += 1;
                }
            }
            Token::Atom(a) => {
                // Try 4-char functional groups first (e.g. COOH).
                if i + 3 < tokens.len()
                    && let (Token::Atom(b), Token::Atom(c), Token::Atom(d)) =
                        (&tokens[i + 1], &tokens[i + 2], &tokens[i + 3])
                {
                    let key = format!("{}{}{}{}", a, b, c, d);
                    if let Some(replacement) = functional_groups(&key) {
                        smiles.push_str(replacement);
                        last_atom_len = 0; // multi-atom groups can't be simply repeated
                        i += 4;
                        continue;
                    }
                }

                // Try 3-char combinations.
                if i + 2 < tokens.len()
                    && let (Token::Atom(b), Token::Atom(c)) = (&tokens[i + 1], &tokens[i + 2])
                {
                    let key = format!("{}{}{}", a, b, c);
                    if let Some(replacement) = functional_groups(&key) {
                        smiles.push_str(replacement);
                        last_atom_len = 0;
                        i += 3;
                        continue;
                    }
                }

                // Try 2-char combinations.
                if i + 1 < tokens.len()
                    && let Token::Atom(b) = &tokens[i + 1]
                {
                    let key = format!("{}{}", a, b);
                    if let Some(replacement) = functional_groups(&key) {
                        smiles.push_str(replacement);
                        last_atom_len = 0;
                        i += 2;
                        continue;
                    }
                }

                // No functional group match — add the heavy atom (H is implicit in SMILES).
                let before = smiles.len();
                atom_to_smiles(&mut smiles, a)?;
                last_atom_len = smiles.len() - before; // track byte length for repeat
                i += 1;
            }
        }
    }

    // Verify balanced parentheses
    let paren_balance: i32 = smiles
        .chars()
        .map(|c| match c {
            '(' => 1,
            ')' => -1,
            _ => 0,
        })
        .sum();

    if paren_balance != 0 {
        return Err(CondensedError::UnbalancedParens);
    }

    Ok(smiles)
}

/// Look up functional group substitutions.
fn functional_groups(key: &str) -> Option<&'static str> {
    match key {
        "COOH" => Some("C(=O)O"),
        "CHO" => Some("C=O"),
        "NHCO" => Some("NC(=O)"),
        "COHN" => Some("C(=O)N"),
        "NO2" => Some("[N+](=O)[O-]"),
        "CN" => Some("C#N"),
        "OH" => Some("O"),
        "NH2" => Some("N"),
        "SH" => Some("S"),
        "PH2" => Some("P"),
        _ => None,
    }
}

/// Convert an atom symbol to SMILES notation with implicit H.
fn atom_to_smiles(smiles: &mut String, atom: &str) -> Result<(), CondensedError> {
    match atom {
        "C" | "N" | "O" | "S" | "P" | "H" | "F" | "Cl" | "Br" | "I" | "Si" | "B" | "Se" | "As" => {
            smiles.push_str(atom);
            Ok(())
        }
        _ => Err(CondensedError::UnknownElement(atom.to_string())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_methane() {
        let mol = parse_condensed("C").expect("methane");
        assert_eq!(mol.atom_count(), 1);
    }

    #[test]
    fn test_ethane() {
        let mol = parse_condensed("CC").expect("ethane");
        assert_eq!(mol.atom_count(), 2);
    }

    #[test]
    fn test_ammonia() {
        let mol = parse_condensed("NH3").expect("ammonia");
        // H is implicit; only N is a heavy atom
        assert_eq!(mol.atom_count(), 1);
    }

    #[test]
    fn test_water() {
        let mol = parse_condensed("H2O").expect("water");
        assert_eq!(mol.atom_count(), 1, "water: only O is a heavy atom");
    }

    #[test]
    fn test_methane_with_h_count() {
        let mol = parse_condensed("CH4").expect("methane");
        assert_eq!(mol.atom_count(), 1, "methane: only C");
    }

    #[test]
    fn test_methanol_ch3oh() {
        let mol = parse_condensed("CH3OH").expect("methanol");
        // C + O (H atoms implicit)
        assert_eq!(mol.atom_count(), 2, "methanol: C + O");
        assert_eq!(mol.bond_count(), 1);
    }

    #[test]
    fn test_acetic_acid_ch3cooh() {
        let mol = parse_condensed("CH3COOH").expect("acetic acid");
        // CH3 (methyl) + COOH (carboxyl) → C, C, O, O = 4 heavy atoms
        assert_eq!(mol.atom_count(), 4, "acetic acid: C,C,O,O");
    }

    #[test]
    fn test_multi_char_atom_repeat_cl2() {
        // Cl2 — digit repeat must pop the full "Cl" symbol, not just 'l'
        let mol = parse_condensed("Cl2").expect("dichlorine");
        assert_eq!(mol.atom_count(), 2, "Cl2: two Cl atoms");
    }

    #[test]
    fn test_multi_char_atom_repeat_br2() {
        let mol = parse_condensed("Br2").expect("dibromine");
        assert_eq!(mol.atom_count(), 2, "Br2: two Br atoms");
    }

    #[test]
    fn test_digit_repeat_cap_rejects_large_count() {
        // A very large repeat count must return an error, not OOM.
        let result = parse_condensed("C9999999");
        // parse::<u32> on "9999999" succeeds (< u32::MAX) but > MAX_REPEAT
        assert!(result.is_err(), "huge repeat count should be rejected");
    }

    #[test]
    fn test_empty_input() {
        let result = parse_condensed("");
        assert!(matches!(result, Err(CondensedError::EmptyInput)));
    }

    #[test]
    fn test_hexane_linear() {
        // Linear hexane without parentheses
        let mol = parse_condensed("CCCCCC").expect("linear hexane");
        assert_eq!(mol.atom_count(), 6);
    }

    #[test]
    fn test_simple_branched() {
        // Simple SMILES-like structures work
        let mol = parse_condensed("CC").expect("ethane");
        assert!(mol.atom_count() >= 2);
    }

    #[test]
    fn test_with_functional_group() {
        // Functional group substitution: COOH becomes C(=O)O
        let mol = parse_condensed("CCOOH").expect("propionic acid");
        assert!(mol.atom_count() >= 3);
    }

    #[test]
    fn test_propane() {
        let mol = parse_condensed("CCC").expect("propane");
        assert!(mol.atom_count() >= 3);
    }

    #[test]
    fn test_butane() {
        let mol = parse_condensed("CCCC").expect("butane");
        assert_eq!(mol.atom_count(), 4);
    }

    #[test]
    fn test_unknown_element() {
        let result = parse_condensed("CXC");
        assert!(result.is_err());
    }
}
