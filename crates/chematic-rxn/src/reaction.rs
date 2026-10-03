use std::collections::HashMap;

use chematic_core::{AtomIdx, Molecule};
use chematic_smiles::{parse as parse_smiles, write as write_smiles};

/// A chemical reaction with reactants, agents, and products.
pub struct Reaction {
    pub reactants: Vec<Molecule>,
    pub agents: Vec<Molecule>,
    pub products: Vec<Molecule>,
}

/// Resource limits for parsing untrusted reaction SMILES.
#[derive(Debug, Clone, Copy)]
pub struct ReactionParseLimits {
    /// Maximum UTF-8 input size in bytes.
    pub max_input_bytes: usize,
    /// Maximum number of non-empty dot-separated components per side.
    pub max_components_per_side: usize,
    /// Maximum atoms in one component.
    pub max_atoms_per_molecule: usize,
    /// Maximum bonds in one component.
    pub max_bonds_per_molecule: usize,
}

impl Default for ReactionParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 16 * 1024 * 1024,
            max_components_per_side: 10_000,
            max_atoms_per_molecule: 100_000,
            max_bonds_per_molecule: 200_000,
        }
    }
}

/// Error type for reaction SMILES parsing.
#[derive(Debug)]
pub enum RxnError {
    /// The string does not contain two `>` delimiters (reaction arrow).
    MissingArrow,
    /// The reaction exceeded a configured resource limit.
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
    /// A SMILES component failed to parse.
    SmilesParse { part: String, source: String },
    /// An atomic-number SMARTS primitive could not be expanded safely.
    UnsupportedAtomicNumberPrimitive { primitive: String },
    /// Expanding aromatic/aliphatic alternatives would exceed the safety cap.
    AtomicNumberExpansionLimit { actual: usize, limit: usize },
    /// A reactant template that is not SMILES-compatible also failed to parse
    /// as SMARTS.
    SmartsParse { part: String, source: String },
    /// A SMARTS-only reactant template uses a feature the reaction engine
    /// does not interpret (tetrahedral or double-bond stereo).
    UnsupportedReactantTemplate { part: String, reason: &'static str },
    /// A product template atom uses a SMARTS feature that cannot be applied
    /// as a product specification.
    UnsupportedProductPrimitive { primitive: String },
}

impl core::fmt::Display for RxnError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::MissingArrow => write!(f, "reaction SMILES must contain '>>'"),
            Self::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(
                f,
                "reaction {resource} exceeds limit {limit} (got {actual})"
            ),
            Self::SmilesParse { part, source } => {
                write!(f, "failed to parse SMILES '{part}': {source}")
            }
            Self::UnsupportedAtomicNumberPrimitive { primitive } => write!(
                f,
                "unsupported atomic-number SMARTS primitive '{primitive}'; expected [#N] or [#N:map]"
            ),
            Self::AtomicNumberExpansionLimit { actual, limit } => write!(
                f,
                "atomic-number SMARTS expansion exceeds limit {limit} (got {actual})"
            ),
            Self::SmartsParse { part, source } => {
                write!(f, "failed to parse reactant SMARTS '{part}': {source}")
            }
            Self::UnsupportedReactantTemplate { part, reason } => {
                write!(f, "unsupported reactant template '{part}': {reason}")
            }
            Self::UnsupportedProductPrimitive { primitive } => write!(
                f,
                "unsupported product template atom '{primitive}': only element, isotope, \
                 charge, chirality, H count (`H<n>` or `;H<n>`) and an atom map can be applied"
            ),
        }
    }
}

/// Expand simple atomic-number SMARTS atoms into SMILES-compatible variants.
///
/// Reaction application uses the SMILES parser for template construction, while
/// the SMARTS parser already accepts query atoms such as `[#7:2]`. This helper
/// is the explicit compatibility boundary between those representations. The
/// atom map is retained verbatim and the returned order is deterministic.
///
/// The supported forms are `[#number]`, `[#number:map]`, and the common
/// explicit-hydrogen forms `[#number;H<n>]` / `[#number;H<n>:map]` (one digit).
/// Any other reactant-side primitive (`[#6;X4:1]`) is left in place for the
/// SMARTS template parser; on the product side it is an explicit error rather
/// than a silent change of semantics.
///
/// Reactant-side atoms (every part before the last `>`) are queries: each
/// primitive expands to its aliphatic form and, where chemically valid, its
/// aromatic form (`[#7:2]` becomes `[N:2]` and `[n:2]`), and the variants
/// enumerate every combination. A string without `>` is treated as a
/// reactant-side pattern.
///
/// Product-side atoms (the last part) are literals, not queries (issue #679):
///
/// * an unmapped primitive becomes one aliphatic atom with implicit
///   hydrogens — organic-subset symbols are written bare (`[#6](=[#8])[#6]`
///   becomes `C(=O)C`), others in brackets (`[#14]` becomes `[Si]`), and
///   `;H1` keeps its explicit hydrogen (`[#6;H1]` becomes `[CH]`);
/// * a mapped primitive keeps the aromaticity of the matched reactant atom:
///   when the reactant side spells the same map as an atomic-number primitive
///   the two share one choice per variant (`[#7:1]>>[#7:1]C` becomes
///   `[N:1]>>[N:1]C` and `[n:1]>>[n:1]C`); otherwise it follows the case of
///   the reactant-side atom with that map, or is aliphatic when there is none.
pub fn expand_atomic_number_primitives(s: &str) -> Result<Vec<String>, RxnError> {
    let s = &normalize_product_query_atoms(s)?;
    expand_atomic_number_primitives_normalized(s)
}

/// Product-template bracket atoms are specifications, not queries. Following
/// RDKit's reaction semantics (issue #734), a product bracket atom that is not
/// valid SMILES is read as SMARTS and rewritten to the parts a product can
/// apply: element (and aromaticity when spelled), isotope, chirality, H count,
/// charge and atom map. Query-only features — `X`, `D`, `R`, `r`, `x`, `v`,
/// `h`, `^`, recursive `$()`, and `,`/`!` alternatives — constrain nothing in
/// a product and are dropped (`[OX2H1:2]` becomes `[OH1:2]`, `[C;H3:1]`
/// becomes `[CH3:1]`, `[#6;X4:1]` becomes `[#6:1]`). A mapped atom whose
/// query names no single element keeps the matched reactant atom's element
/// (`[c,n:1]` becomes `[*:1]`); an unmapped one cannot be built and is an
/// [`RxnError::UnsupportedProductPrimitive`] error rather than a dummy atom.
/// Reactant and agent parts, and every SMILES-valid product atom, are
/// returned unchanged.
/// An atom-map label's number (`"01"` and `"1"` are the same map); a label
/// that is not a number keeps a distinct sentinel.
fn map_number(label: &str) -> u64 {
    label.parse().unwrap_or(u64::MAX)
}

/// Byte offsets of the reaction-arrow `>` separators in `s`: top-level `>`
/// characters outside bracket atoms that are not the head of a dative bond
/// `->` (a `->` directly followed by `>` is a bond then the arrow).
pub(crate) fn reaction_separators(s: &str) -> Vec<usize> {
    let b = s.as_bytes();
    let mut depth = 0usize;
    let mut out: Vec<usize> = Vec::new();
    for (i, &c) in b.iter().enumerate() {
        match c {
            b'[' => depth += 1,
            b']' => depth = depth.saturating_sub(1),
            b'>' if depth == 0 => {
                let dative = i > 0
                    && b[i - 1] == b'-'
                    && b.get(i + 1) != Some(&b'>')
                    && out.last() != Some(&(i - 1));
                if !dative {
                    out.push(i);
                }
            }
            _ => {}
        }
    }
    out
}

/// `reactants>agents>products` split at the first two arrow separators
/// (see [`reaction_separators`]); `None` without two separators. The
/// product part keeps any further `>` for the caller to reject.
pub(crate) fn split_reaction_parts(s: &str) -> Option<[&str; 3]> {
    let seps = reaction_separators(s);
    let (&a, &b) = (seps.first()?, seps.get(1)?);
    Some([&s[..a], &s[a + 1..b], &s[b + 1..]])
}

pub fn normalize_product_query_atoms(s: &str) -> Result<String, RxnError> {
    // Product side = text after the last top-level `>` (none: nothing to do).
    let Some(product_start) = reaction_separators(s).last().map(|i| i + 1) else {
        return Ok(s.to_string());
    };
    let mut out = String::with_capacity(s.len());
    out.push_str(&s[..product_start]);
    let mut rest = &s[product_start..];
    while let Some(open) = rest.find('[') {
        out.push_str(&rest[..open]);
        // The matching `]`, skipping brackets nested in a recursive `$()`.
        let mut depth = 0usize;
        let mut close = None;
        for (i, b) in rest[open..].bytes().enumerate() {
            match b {
                b'[' => depth += 1,
                b']' => {
                    depth -= 1;
                    if depth == 0 {
                        close = Some(open + i);
                        break;
                    }
                }
                _ => {}
            }
        }
        let Some(close) = close else {
            out.push_str(&rest[open..]);
            return Ok(out);
        };
        let atom = &rest[open..=close];
        rest = &rest[close + 1..];
        out.push_str(&product_atom_spec(atom)?);
    }
    out.push_str(rest);
    Ok(out)
}

/// One product bracket atom rewritten per [`normalize_product_query_atoms`].
fn product_atom_spec(atom: &str) -> Result<String, RxnError> {
    use chematic_smarts::{AtomPrimitive, AtomQuery};

    fn conjuncts<'q>(q: &'q AtomQuery, out: &mut Vec<&'q AtomQuery>) {
        match q {
            AtomQuery::And(a, b) => {
                conjuncts(a, out);
                conjuncts(b, out);
            }
            other => out.push(other),
        }
    }
    fn put<T: PartialEq>(slot: &mut Option<T>, value: T, conflict: &mut bool) {
        if slot.as_ref().is_some_and(|v| *v != value) {
            *conflict = true;
        }
        *slot = Some(value);
    }

    // SMILES-valid atoms are already specifications.
    if parse_smiles(atom).is_ok() {
        return Ok(atom.to_string());
    }
    let unsupported = || RxnError::UnsupportedProductPrimitive {
        primitive: atom.to_string(),
    };
    let query = chematic_smarts::parse_smarts(atom).map_err(|_| unsupported())?;
    if query.atoms.len() != 1 {
        return Err(unsupported());
    }
    let map = query.atoms[0].atom_map;

    let mut parts = Vec::new();
    conjuncts(&query.atoms[0].query, &mut parts);
    let mut symbol: Option<String> = None;
    let mut atomic_number: Option<u8> = None;
    let mut aromatic: Option<bool> = None;
    let mut isotope: Option<u16> = None;
    let mut chirality: Option<u8> = None;
    let mut hcount: Option<u8> = None;
    let mut charge: Option<i8> = None;
    let mut conflict = false;
    for part in parts {
        // `,` / `!` alternatives are query-only.
        let AtomQuery::Primitive(p) = part else {
            continue;
        };
        match p {
            AtomPrimitive::Symbol(sym) => put(&mut symbol, sym.clone(), &mut conflict),
            AtomPrimitive::AtomicNum(n) => put(&mut atomic_number, *n, &mut conflict),
            AtomPrimitive::Aromatic(a) => put(&mut aromatic, *a, &mut conflict),
            AtomPrimitive::Isotope(m) => put(&mut isotope, *m, &mut conflict),
            AtomPrimitive::Chirality(c) => put(&mut chirality, *c, &mut conflict),
            AtomPrimitive::HCount(h) => put(&mut hcount, *h, &mut conflict),
            AtomPrimitive::Charge(c) => put(&mut charge, *c, &mut conflict),
            // X, D, R, r, x, v, h, ^, $(), *: query-only.
            _ => {}
        }
    }
    if conflict {
        return Err(unsupported());
    }
    let symbol_absent = symbol.is_none();
    let element = match (symbol, atomic_number) {
        (Some(sym), _) => {
            let e = chematic_core::Element::from_symbol(&sym).ok_or_else(unsupported)?;
            if atomic_number.is_some_and(|n| n != e.atomic_number()) {
                return Err(unsupported());
            }
            Some(e)
        }
        (None, Some(n)) => {
            Some(chematic_core::Element::from_atomic_number(n).ok_or_else(unsupported)?)
        }
        (None, None) => None,
    };

    // An atomic-number atom with unspelled aromaticity and nothing but an H
    // count, a charge and a map stays `[#n(;Hh)(;±c)(:m)]`, so the
    // atomic-number expansion picks its spelling from the reactant side (#679).
    if let Some(e) = element
        && symbol_absent
        && aromatic.is_none()
        && isotope.is_none()
        && chirality.is_none()
        && hcount.is_none_or(|h| h <= 9)
        && charge.is_none_or(|c| (-9..=9).contains(&c))
    {
        let mut spec = format!("[#{}", e.atomic_number());
        if let Some(h) = hcount {
            spec.push_str(&format!(";H{h}"));
        }
        match charge {
            Some(c) if c >= 0 => spec.push_str(&format!(";+{c}")),
            Some(c) => spec.push_str(&format!(";-{}", -i16::from(c))),
            None => {}
        }
        if let Some(m) = map {
            spec.push_str(&format!(":{m}"));
        }
        spec.push(']');
        return Ok(spec);
    }

    let mut spec = String::from("[");
    if let Some(m) = isotope {
        spec.push_str(&m.to_string());
    }
    match element {
        Some(e) if aromatic == Some(true) => spec.push_str(&e.symbol().to_ascii_lowercase()),
        Some(e) => spec.push_str(e.symbol()),
        // Mapped: the reactant atom's element is kept (RDKit semantics).
        None if map.is_some()
            && isotope.is_none()
            && chirality.is_none()
            && hcount.is_none()
            && charge.is_none() =>
        {
            spec.push('*')
        }
        None => return Err(unsupported()),
    }
    match chirality {
        Some(1) => spec.push('@'),
        Some(2) => spec.push_str("@@"),
        _ => {}
    }
    if let Some(h) = hcount {
        spec.push_str(&format!("H{h}"));
    }
    match charge {
        Some(0) => spec.push_str("+0"),
        Some(c) if c > 0 => spec.push_str(&format!("+{c}")),
        Some(c) => spec.push_str(&format!("-{}", -i16::from(c))),
        None => {}
    }
    if let Some(m) = map {
        spec.push_str(&format!(":{m}"));
    }
    spec.push(']');
    Ok(spec)
}

fn expand_atomic_number_primitives_normalized(s: &str) -> Result<Vec<String>, RxnError> {
    const MAX_VARIANTS: usize = 256;
    const ORGANIC_SUBSET: [&str; 10] = ["B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"];

    // Split into `>`-separated parts outside brackets; the last part is the
    // product side when there is at least one `>`.
    let bytes = s.as_bytes();
    let mut part_starts = vec![0usize];
    part_starts.extend(reaction_separators(s).into_iter().map(|i| i + 1));
    let has_products = part_starts.len() > 1;
    let product_start = *part_starts.last().unwrap_or(&0);
    let reactant_text = if has_products {
        &s[..product_start]
    } else {
        ""
    };

    // Aromaticity spelled for a mapped reactant-side bracket atom (`[n:1]` →
    // aromatic, `[N:1]` → aliphatic); `None` when the map is absent or the
    // atom is itself an atomic-number primitive.
    let reactant_case = |map: &str| -> Option<bool> {
        let mut rest = reactant_text;
        while let Some(open) = rest.find('[') {
            let Some(close_rel) = rest[open + 1..].find(']') else {
                break;
            };
            let inner = &rest[open + 1..open + 1 + close_rel];
            rest = &rest[open + 1 + close_rel + 1..];
            // Atom maps compare as numbers (`:01` is map 1).
            if inner.rsplit_once(':').map(|(_, m)| map_number(m)) != Some(map_number(map)) {
                continue;
            }
            let body = inner.trim_start_matches(|c: char| c.is_ascii_digit());
            match body.bytes().next() {
                Some(b'#') => return None,
                Some(c) if c.is_ascii_lowercase() => return Some(true),
                Some(c) if c.is_ascii_uppercase() => return Some(false),
                _ => return None,
            }
        }
        None
    };

    // Text segments and expansion slots, in order. A slot holds its
    // alternatives and the index of the choice group it belongs to (mapped
    // primitives with the same map share a group).
    enum Piece {
        Text(String),
        Slot { options: Vec<String>, group: usize },
    }
    let mut pieces: Vec<Piece> = Vec::new();
    let mut groups: Vec<usize> = Vec::new(); // option count per group
    let mut group_of_map: std::collections::HashMap<String, usize> =
        std::collections::HashMap::new();
    let mut text = String::new();

    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] != b'[' {
            let next = s[i..].find('[').map_or(s.len(), |n| i + n);
            text.push_str(&s[i..next]);
            i = next;
            continue;
        }
        let start = i;
        let Some(end_rel) = s[i + 1..].find(']') else {
            text.push_str(&s[i..]);
            break;
        };
        let end = i + 1 + end_rel;
        let primitive = &s[start..=end];
        let inner = &s[start + 1..end];
        if !inner.starts_with('#') {
            text.push_str(primitive);
            i = end + 1;
            continue;
        }
        let number_end = inner[1..]
            .bytes()
            .position(|b| !b.is_ascii_digit())
            .map(|p| p + 1)
            .unwrap_or(inner.len());
        let number_text = &inner[1..number_end];
        let suffix = &inner[number_end..];
        // Suffix: optional H count (`H`, `H2`) and charge (`+`, `++`, `-2`),
        // each optionally `;`-separated, then an optional `:map`.
        let (body, map_suffix) = match suffix.rfind(':') {
            Some(p)
                if p + 1 < suffix.len() && suffix[p + 1..].bytes().all(|b| b.is_ascii_digit()) =>
            {
                (&suffix[..p], &suffix[p..])
            }
            _ => (suffix, ""),
        };
        let mut hydrogen_count: Option<u8> = None;
        let mut charge: Option<i8> = None;
        let mut simple_suffix = true;
        let mut rest = body.as_bytes();
        while !rest.is_empty() {
            if rest[0] == b';' {
                rest = &rest[1..];
                continue;
            }
            match rest[0] {
                b'H' if hydrogen_count.is_none() => {
                    let digit = rest.get(1).filter(|b| b.is_ascii_digit());
                    hydrogen_count = Some(digit.map_or(1, |d| d - b'0'));
                    rest = &rest[1 + usize::from(digit.is_some())..];
                }
                sign @ (b'+' | b'-') if charge.is_none() => {
                    let mut n: i8 = 1;
                    let mut used = 1;
                    if let Some(d) = rest.get(1).filter(|b| b.is_ascii_digit()) {
                        n = (d - b'0') as i8;
                        used = 2;
                    } else {
                        while rest.get(used) == Some(&sign) && n < 9 {
                            n += 1;
                            used += 1;
                        }
                    }
                    charge = Some(if sign == b'+' { n } else { -n });
                    rest = &rest[used..];
                }
                _ => {
                    simple_suffix = false;
                    break;
                }
            }
        }
        let on_product_side = has_products && start >= product_start;
        if number_text.is_empty() || !simple_suffix {
            if !on_product_side && !number_text.is_empty() {
                // A compound reactant-side primitive (`[#6;X4:1]`) is a query
                // the SMARTS parser handles; leave it for the template parser.
                text.push_str(primitive);
                i = end + 1;
                continue;
            }
            return Err(RxnError::UnsupportedAtomicNumberPrimitive {
                primitive: primitive.to_string(),
            });
        }
        let atomic_number = number_text.parse::<u8>().ok();
        let Some(atomic_number) = atomic_number.filter(|n| (1..=118).contains(n)) else {
            return Err(RxnError::UnsupportedAtomicNumberPrimitive {
                primitive: primitive.to_string(),
            });
        };
        let Some(element) = chematic_core::Element::from_atomic_number(atomic_number) else {
            return Err(RxnError::UnsupportedAtomicNumberPrimitive {
                primitive: primitive.to_string(),
            });
        };
        let symbol = element.symbol();
        let mut hydrogen = match hydrogen_count {
            None => String::new(),
            Some(1) => "H".to_string(),
            Some(n) => format!("H{n}"),
        };
        match charge {
            Some(0) => hydrogen.push_str("+0"),
            Some(1) => hydrogen.push('+'),
            Some(-1) => hydrogen.push('-'),
            Some(c) if c > 0 => hydrogen.push_str(&format!("+{c}")),
            Some(c) => hydrogen.push_str(&format!("-{}", -i16::from(c))),
            None => {}
        }
        let can_be_aromatic = matches!(atomic_number, 5 | 6 | 7 | 8 | 15 | 16);
        let aliphatic = format!("[{symbol}{hydrogen}{map_suffix}]");
        let aromatic = format!(
            "[{}{}{}]",
            symbol.to_ascii_lowercase(),
            hydrogen,
            map_suffix
        );
        let map_key = map_suffix
            .strip_prefix(':')
            .map(|m| map_number(m).to_string());
        let map = map_key.as_deref();

        let options: Vec<String> = if !on_product_side {
            if can_be_aromatic {
                vec![aliphatic, aromatic]
            } else {
                vec![aliphatic]
            }
        } else if let Some(map) = map {
            if group_of_map.contains_key(map) {
                // Shares the reactant-side choice (same option order).
                if can_be_aromatic {
                    vec![aliphatic, aromatic]
                } else {
                    vec![aliphatic]
                }
            } else {
                match reactant_case(map) {
                    Some(true) if can_be_aromatic => vec![aromatic],
                    _ => vec![aliphatic],
                }
            }
        } else if hydrogen_count.is_none() && charge.is_none() && ORGANIC_SUBSET.contains(&symbol) {
            vec![symbol.to_string()]
        } else {
            vec![aliphatic]
        };

        if options.len() == 1 {
            text.push_str(&options[0]);
            i = end + 1;
            continue;
        }
        let group = match map {
            Some(map) if !on_product_side || group_of_map.contains_key(map) => {
                *group_of_map.entry(map.to_string()).or_insert_with(|| {
                    groups.push(options.len());
                    groups.len() - 1
                })
            }
            _ => {
                groups.push(options.len());
                groups.len() - 1
            }
        };
        if !text.is_empty() {
            pieces.push(Piece::Text(std::mem::take(&mut text)));
        }
        pieces.push(Piece::Slot { options, group });
        i = end + 1;
    }
    if !text.is_empty() {
        pieces.push(Piece::Text(text));
    }

    let total: usize = groups
        .iter()
        .try_fold(1usize, |acc, &n| {
            let next = acc.saturating_mul(n);
            (next <= MAX_VARIANTS).then_some(next)
        })
        .unwrap_or(usize::MAX);
    if total > MAX_VARIANTS {
        return Err(RxnError::AtomicNumberExpansionLimit {
            actual: total,
            limit: MAX_VARIANTS,
        });
    }

    // Enumerate choice vectors in lexicographic order (first group slowest).
    let mut choice = vec![0usize; groups.len()];
    let mut variants = Vec::with_capacity(total);
    loop {
        let mut out = String::with_capacity(s.len());
        for piece in &pieces {
            match piece {
                Piece::Text(t) => out.push_str(t),
                Piece::Slot { options, group } => out.push_str(&options[choice[*group]]),
            }
        }
        variants.push(out);
        let mut k = groups.len();
        loop {
            if k == 0 {
                return Ok(variants);
            }
            k -= 1;
            choice[k] += 1;
            if choice[k] < groups[k] {
                break;
            }
            choice[k] = 0;
        }
    }
}

impl std::error::Error for RxnError {}

/// Parse a reaction SMILES string of the form `"reactants>agents>products"`.
///
/// - `>` splits the string into 3 parts (reactants, agents, products).
/// - Each part is a dot-separated list of SMILES; empty parts yield empty `Vec`s.
/// - `"R>>P"` is the standard form with empty agents section.
/// - Returns `Err(RxnError::MissingArrow)` if fewer than two `>` delimiters are found.
pub fn parse_reaction(s: &str) -> Result<Reaction, RxnError> {
    parse_reaction_with_limits(s, &ReactionParseLimits::default())
}

/// Parse a reaction SMILES while enforcing input and per-component limits.
pub fn parse_reaction_with_limits(
    s: &str,
    limits: &ReactionParseLimits,
) -> Result<Reaction, RxnError> {
    if s.len() > limits.max_input_bytes {
        return Err(RxnError::ResourceLimit {
            resource: "input bytes",
            actual: s.len(),
            limit: limits.max_input_bytes,
        });
    }
    // Three parts need two arrow separators (a dative `->` is not one).
    let Some(parts) = split_reaction_parts(s) else {
        return Err(RxnError::MissingArrow);
    };

    let parse_part = |side: &'static str, s: &str| -> Result<Vec<Molecule>, RxnError> {
        if s.is_empty() {
            return Ok(Vec::new());
        }
        let components: Vec<&str> = s.split('.').filter(|p| !p.is_empty()).collect();
        if components.len() > limits.max_components_per_side {
            return Err(RxnError::ResourceLimit {
                resource: side,
                actual: components.len(),
                limit: limits.max_components_per_side,
            });
        }
        components
            .into_iter()
            .filter(|p| !p.is_empty())
            .map(|p| {
                let molecule = parse_smiles(p).map_err(|e| RxnError::SmilesParse {
                    part: p.to_string(),
                    source: e.to_string(),
                })?;
                if molecule.atom_count() > limits.max_atoms_per_molecule {
                    return Err(RxnError::ResourceLimit {
                        resource: "atoms per molecule",
                        actual: molecule.atom_count(),
                        limit: limits.max_atoms_per_molecule,
                    });
                }
                if molecule.bond_count() > limits.max_bonds_per_molecule {
                    return Err(RxnError::ResourceLimit {
                        resource: "bonds per molecule",
                        actual: molecule.bond_count(),
                        limit: limits.max_bonds_per_molecule,
                    });
                }
                Ok(molecule)
            })
            .collect()
    };

    Ok(Reaction {
        reactants: parse_part("reactants", parts[0])?,
        agents: parse_part("agents", parts[1])?,
        products: parse_part("products", parts[2])?,
    })
}

/// Serialize a `Reaction` back to a reaction SMILES string.
pub fn write_reaction(rxn: &Reaction) -> String {
    let join = |mols: &[Molecule]| -> String {
        mols.iter().map(write_smiles).collect::<Vec<_>>().join(".")
    };
    format!(
        "{}>{}>{}",
        join(&rxn.reactants),
        join(&rxn.agents),
        join(&rxn.products),
    )
}

/// The center of a chemical reaction: changed atoms and bonds.
#[derive(Debug, Clone)]
pub struct ReactionCenter {
    /// Bonds present in reactants but not in products (broken bonds).
    pub broken_bonds: Vec<(AtomIdx, AtomIdx)>,
    /// Bonds present in products but not in reactants (formed bonds).
    pub formed_bonds: Vec<(AtomIdx, AtomIdx)>,
    /// Atoms whose element, charge, or aromaticity changed between reactants and products.
    /// Uses reactant-side indexing.
    pub changed_atoms: Vec<AtomIdx>,
}

/// Identify the reaction center: bonds broken/formed and atoms changed.
///
/// Uses atom_map numbers (if present) to match reactant atoms to product atoms.
/// For each mapped atom pair:
/// - Compares bond connectivity to identify broken/formed bonds.
/// - Compares element, charge, aromaticity to identify changed atoms.
///
/// Returns empty vecs if atoms lack atom_map annotations.
pub fn find_reaction_center(rxn: &Reaction) -> ReactionCenter {
    let mut broken_bonds = Vec::new();
    let mut formed_bonds = Vec::new();
    let mut changed_atoms = Vec::new();

    // Build atom_map -> (mol_idx, atom_idx) for reactants and products
    let mut reactant_map: HashMap<u16, (usize, AtomIdx)> = HashMap::new();
    for (mol_idx, mol) in rxn.reactants.iter().enumerate() {
        for (atom_idx, atom) in mol.atoms() {
            if let Some(map_num) = atom.atom_map {
                reactant_map.insert(map_num, (mol_idx, atom_idx));
            }
        }
    }

    let mut product_map: HashMap<u16, (usize, AtomIdx)> = HashMap::new();
    for (mol_idx, mol) in rxn.products.iter().enumerate() {
        for (atom_idx, atom) in mol.atoms() {
            if let Some(map_num) = atom.atom_map {
                product_map.insert(map_num, (mol_idx, atom_idx));
            }
        }
    }

    // If no atom_map, return empty
    if reactant_map.is_empty() {
        return ReactionCenter {
            broken_bonds,
            formed_bonds,
            changed_atoms,
        };
    }

    // Identify broken bonds (edges in reactants not in products)
    for map_num in reactant_map.keys() {
        let (r_mol_idx, r_atom_idx) = reactant_map[map_num];
        let r_mol = &rxn.reactants[r_mol_idx];
        for (r_neighbor, _) in r_mol.neighbors(r_atom_idx) {
            if let Some(neighbor_map) = r_mol.atom(r_neighbor).atom_map
                && neighbor_map > *map_num
            {
                // Check if this bond exists in products
                if let Some((p_mol_idx, p_atom_idx)) = product_map.get(map_num) {
                    let p_mol = &rxn.products[*p_mol_idx];
                    if let Some((_, p_neighbor)) = product_map.get(&neighbor_map) {
                        let bond_exists = p_mol.bond_between(*p_atom_idx, *p_neighbor).is_some();
                        if !bond_exists {
                            broken_bonds.push((r_atom_idx, r_neighbor));
                        }
                    } else {
                        broken_bonds.push((r_atom_idx, r_neighbor));
                    }
                } else {
                    broken_bonds.push((r_atom_idx, r_neighbor));
                }
            }
        }
    }

    // Identify formed bonds (edges in products not in reactants)
    for map_num in product_map.keys() {
        let (p_mol_idx, p_atom_idx) = product_map[map_num];
        let p_mol = &rxn.products[p_mol_idx];
        for (p_neighbor, _) in p_mol.neighbors(p_atom_idx) {
            if let Some(neighbor_map) = p_mol.atom(p_neighbor).atom_map
                && neighbor_map > *map_num
            {
                // Check if this bond exists in reactants
                if let Some((r_mol_idx, r_atom_idx)) = reactant_map.get(map_num) {
                    let r_mol = &rxn.reactants[*r_mol_idx];
                    if let Some((_, r_neighbor)) = reactant_map.get(&neighbor_map) {
                        let bond_exists = r_mol.bond_between(*r_atom_idx, *r_neighbor).is_some();
                        if !bond_exists {
                            formed_bonds.push((p_atom_idx, p_neighbor));
                        }
                    } else {
                        formed_bonds.push((p_atom_idx, p_neighbor));
                    }
                } else {
                    formed_bonds.push((p_atom_idx, p_neighbor));
                }
            }
        }
    }

    // Identify changed atoms (element/charge/aromaticity changes)
    for map_num in reactant_map.keys() {
        let (r_mol_idx, r_atom_idx) = reactant_map[map_num];
        let r_atom = rxn.reactants[r_mol_idx].atom(r_atom_idx);

        if let Some((p_mol_idx, p_atom_idx)) = product_map.get(map_num) {
            let p_atom = rxn.products[*p_mol_idx].atom(*p_atom_idx);

            if r_atom.element != p_atom.element
                || r_atom.charge != p_atom.charge
                || r_atom.aromatic != p_atom.aromatic
            {
                changed_atoms.push(r_atom_idx);
            }
        }
    }

    ReactionCenter {
        broken_bonds,
        formed_bonds,
        changed_atoms,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::AtomIdx;

    #[test]
    fn dative_bonds_are_not_reaction_arrows(/* issue #734 bug check */) {
        assert_eq!(reaction_separators("[N:1]->[Pt:2]>>[N:1]"), vec![13, 14]);
        assert_eq!(reaction_separators("[Pt]<-N>[Pd]>N->[Pt]"), vec![7, 12]);
        // `->` directly followed by `>` is a bond end then the arrow.
        assert_eq!(reaction_separators("C->>C"), vec![2, 3]);
        assert_eq!(
            split_reaction_parts("N->[Pt]>>N.[Pt]"),
            Some(["N->[Pt]", "", "N.[Pt]"])
        );
        let rxn = parse_reaction("N->[Pt]>>N.[Pt]").unwrap();
        assert_eq!(rxn.reactants.len(), 1);
        assert_eq!(rxn.products.len(), 2);
        let query = crate::query::parse_reaction_smarts("[N:1]->[Pt:2]>>[N:1].[Pt:2]").unwrap();
        assert_eq!(query.reactant_patterns.len(), 1);
        assert_eq!(query.product_patterns.len(), 2);
    }

    #[test]
    fn test_simple_reaction() {
        let rxn = parse_reaction("[Na+].[Cl-]>>[Na+].[Cl-]").unwrap();
        assert_eq!(rxn.reactants.len(), 2);
        assert_eq!(rxn.agents.len(), 0);
        assert_eq!(rxn.products.len(), 2);
    }

    #[test]
    fn test_one_to_one() {
        let rxn = parse_reaction("CC>>CC").unwrap();
        assert_eq!(rxn.reactants.len(), 1);
        assert_eq!(rxn.agents.len(), 0);
        assert_eq!(rxn.products.len(), 1);
    }

    #[test]
    fn test_with_agent() {
        let rxn = parse_reaction("CC>O>CC=O").unwrap();
        assert_eq!(rxn.reactants.len(), 1);
        assert_eq!(rxn.agents.len(), 1);
        assert_eq!(rxn.products.len(), 1);
        // agent is water (O) — one heavy atom
        assert_eq!(rxn.agents[0].atom_count(), 1);
    }

    #[test]
    fn test_two_reactants() {
        let rxn = parse_reaction("C.C>>CC").unwrap();
        assert_eq!(rxn.reactants.len(), 2);
        assert_eq!(rxn.agents.len(), 0);
        assert_eq!(rxn.products.len(), 1);
    }

    #[test]
    fn test_empty_reaction() {
        let rxn = parse_reaction(">>").unwrap();
        assert_eq!(rxn.reactants.len(), 0);
        assert_eq!(rxn.agents.len(), 0);
        assert_eq!(rxn.products.len(), 0);
    }

    #[test]
    fn test_missing_arrow_error() {
        // "CC>CC" has only one ">" → splitn(3, '>') gives 2 parts → MissingArrow
        let err = parse_reaction("CC>CC");
        assert!(matches!(err, Err(RxnError::MissingArrow)));
    }

    #[test]
    fn test_invalid_smiles_error() {
        // "[X]" has an unrecognised element symbol → SmilesParse error
        let err = parse_reaction("[X]>>C");
        assert!(matches!(err, Err(RxnError::SmilesParse { .. })));
    }

    #[test]
    fn test_atom_map_preserved() {
        // [CH3:1] has atom_map = 1
        let rxn = parse_reaction("[CH3:1]>>[CH3:1]").unwrap();
        assert_eq!(rxn.reactants[0].atom(AtomIdx(0)).atom_map, Some(1));
        assert_eq!(rxn.products[0].atom(AtomIdx(0)).atom_map, Some(1));
    }

    #[test]
    fn test_product_atom_count() {
        let rxn = parse_reaction("CC>>CCO").unwrap();
        assert_eq!(rxn.products[0].atom_count(), 3); // C, C, O
    }

    #[test]
    fn test_empty_agents() {
        let rxn = parse_reaction("C>>C").unwrap();
        assert_eq!(rxn.agents.len(), 0);
    }

    #[test]
    fn test_complex_reaction() {
        // aspirin synthesis-like
        let rxn = parse_reaction("OC(=O)c1ccccc1.CC(=O)O>>CC(=O)Oc1ccccc1C(=O)O");
        assert!(rxn.is_ok());
        let rxn = rxn.unwrap();
        assert_eq!(rxn.reactants.len(), 2);
        assert_eq!(rxn.products.len(), 1);
    }

    #[test]
    fn test_reactant_only() {
        let rxn = parse_reaction("C>>").unwrap();
        assert_eq!(rxn.reactants.len(), 1);
        assert_eq!(rxn.products.len(), 0);
    }

    #[test]
    fn test_multi_product() {
        let rxn = parse_reaction(">>C.C").unwrap();
        assert_eq!(rxn.products.len(), 2);
    }

    #[test]
    fn test_roundtrip() {
        let s = "CC>>CC=O";
        let rxn = parse_reaction(s).unwrap();
        let written = write_reaction(&rxn);
        let rxn2 = parse_reaction(&written).unwrap();
        assert_eq!(rxn.reactants.len(), rxn2.reactants.len());
        assert_eq!(rxn.products.len(), rxn2.products.len());
        assert_eq!(
            rxn.reactants[0].atom_count(),
            rxn2.reactants[0].atom_count()
        );
    }

    #[test]
    fn test_write_reaction_format() {
        let rxn = parse_reaction("CC>O>CC=O").unwrap();
        let s = write_reaction(&rxn);
        // Standard reaction SMILES format: "reactants>agents>products"
        assert!(s.contains('>'), "written reaction should contain '>'");
        // Should have exactly 2 '>' characters
        assert_eq!(s.chars().filter(|&c| c == '>').count(), 2);
    }

    #[test]
    fn test_parse_reaction_limits_reject_input_and_components() {
        let limits = ReactionParseLimits {
            max_input_bytes: 3,
            ..ReactionParseLimits::default()
        };
        assert!(matches!(
            parse_reaction_with_limits("C>>C", &limits),
            Err(RxnError::ResourceLimit {
                resource: "input bytes",
                ..
            })
        ));

        let limits = ReactionParseLimits {
            max_components_per_side: 1,
            ..ReactionParseLimits::default()
        };
        assert!(matches!(
            parse_reaction_with_limits("C.C>>C", &limits),
            Err(RxnError::ResourceLimit {
                resource: "reactants",
                ..
            })
        ));
    }

    #[test]
    fn test_parse_reaction_limits_reject_component_size() {
        let limits = ReactionParseLimits {
            max_atoms_per_molecule: 1,
            ..ReactionParseLimits::default()
        };
        assert!(matches!(
            parse_reaction_with_limits("CC>>C", &limits),
            Err(RxnError::ResourceLimit {
                resource: "atoms per molecule",
                ..
            })
        ));
    }

    #[test]
    fn expands_atomic_number_primitives_and_preserves_maps() {
        let variants = expand_atomic_number_primitives("[#7:2][#6]").unwrap();
        assert_eq!(
            variants,
            vec![
                "[N:2][C]".to_string(),
                "[N:2][c]".to_string(),
                "[n:2][C]".to_string(),
                "[n:2][c]".to_string(),
            ]
        );
    }

    #[test]
    fn expands_multi_hydrogen_constraint() {
        // `;H<n>` with any single digit is a supported H count now (#734).
        assert_eq!(
            expand_atomic_number_primitives("[#7;H2:1]").unwrap(),
            vec!["[NH2:1]".to_string(), "[nH2:1]".to_string()],
        );
    }

    #[test]
    fn leaves_compound_reactant_primitive_for_smarts_parser() {
        // A non-H compound primitive on the reactant side is passed through
        // verbatim for the SMARTS template parser (#734)...
        assert_eq!(
            expand_atomic_number_primitives("[#7;X3:1]>>[N:1]").unwrap(),
            vec!["[#7;X3:1]>>[N:1]".to_string()],
        );
        // ...and on the product side its query-only part is dropped, as in
        // RDKit: a product atom is a specification, not a query.
        assert_eq!(
            expand_atomic_number_primitives("[#7:1]>>[#7;X3:1]").unwrap(),
            vec!["[N:1]>>[N:1]".to_string(), "[n:1]>>[n:1]".to_string()],
        );
        // An unmapped product atom naming no single element cannot be built.
        assert!(matches!(
            expand_atomic_number_primitives("[#7:1]>>[#7:1][#6,#7]"),
            Err(RxnError::UnsupportedProductPrimitive { .. })
        ));
    }

    #[test]
    fn product_query_atoms_reduce_to_specifications() {
        // Issue #734: RDKit-style product templates. Query-only features are
        // dropped; element, aromaticity, isotope, chirality, H count, charge
        // and map are kept.
        for (smirks, expected) in [
            (
                "[c:1][OX2:2][CH3:3]>>[c:1][OX2H1:2]",
                "[c:1][OX2:2][CH3:3]>>[c:1][OH1:2]",
            ),
            ("[C:1]>>[C;H3:1]", "[C:1]>>[CH3:1]"),
            ("[C:1]>>[CX4;!R:1]", "[C:1]>>[C:1]"),
            ("[N:1]>>[N+;H3:1]", "[N:1]>>[NH3+1:1]"),
            ("[C:1]>>[13C@;X4:1]", "[C:1]>>[13C@:1]"),
            ("[O:1]>>[O:1][CX4]", "[O:1]>>[O:1][C]"),
            ("[#6;X4:1]>>[#6;X4:1]O", "[#6;X4:1]>>[#6:1]O"),
            // A mapped atom without a single element keeps the reactant's.
            ("[C:1]>>[C,N:1]", "[C:1]>>[*:1]"),
            ("[C:1]>>[$(C):1]", "[C:1]>>[*:1]"),
            // Reactant side and SMILES-valid product atoms are untouched.
            ("[CX4:1][OH]>>[C:1]=O", "[CX4:1][OH]>>[C:1]=O"),
        ] {
            assert_eq!(
                normalize_product_query_atoms(smirks).unwrap(),
                expected,
                "{smirks}"
            );
        }
        for smirks in ["[O:1]>>[O:1][C,N]", "[O:1]>>[O:1][$(C)]", "[C:1]>>[C;N:1]"] {
            assert!(
                matches!(
                    normalize_product_query_atoms(smirks),
                    Err(RxnError::UnsupportedProductPrimitive { .. })
                ),
                "{smirks}"
            );
        }
    }

    #[test]
    fn expands_single_hydrogen_constraint_and_preserves_maps() {
        let variants = expand_atomic_number_primitives("[#7;H1:2][#6;H1]").unwrap();
        assert_eq!(
            variants,
            vec![
                "[NH:2][CH]".to_string(),
                "[NH:2][cH]".to_string(),
                "[nH:2][CH]".to_string(),
                "[nH:2][cH]".to_string(),
            ]
        );
    }
}
