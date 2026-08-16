//! CIF explicit symmetry-operation parsing and unit-cell expansion.
//!
//! Everything in this module works from symmetry operations **literally
//! written in a CIF's own `_space_group_symop_operation_xyz` /
//! `_symmetry_equiv_pos_as_xyz` (and their modern dotted-tag aliases) --
//! never from a space-group name or International Tables number. See
//! `crates/chematic-mol/src/cif.rs`'s [`crate::cif::CifSymmetryStatus`] docs
//! for the public-facing distinction this buys.
//!
//! IUCr's symmetry-operation convention (`X' = W*X + w`, a 3x3 rotation
//! part `W` plus a translation part `w`) is public crystallographic
//! convention (International Tables for Crystallography, Vol. A) -- cited
//! as that, not copied from any specific tool's implementation.
//!
//! # Module shape
//!
//! - [`Rational`]: a tiny internal exact-fraction type (no external crate).
//! - [`CifSymmetryOperation`]: one parsed operation (rotation + translation).
//! - [`CifSymmetryError`]: every way parsing/expansion can fail, fail-closed.
//! - [`scan_operation_sources`]: tokenizer-level extraction of raw
//!   `(id, expression text)` rows from every place an operation-list tag
//!   alias can appear (loop or standalone data item).
//! - [`resolve_symmetry_operations`]: parses, cross-checks multiple tag
//!   aliases against each other, and validates (duplicates, identity).
//! - [`expand_sites`]: applies a resolved operation list to an asymmetric
//!   unit, deduplicating special positions via
//!   [`chematic_crystal::minimum_image`].

use crate::cif::CifPeriodicError;
use chematic_crystal::{FractionalCoord, Lattice, PeriodicSite, SiteSpecies, minimum_image};

// ---------------------------------------------------------------------------
// Tolerance -- the ONE canonical site-merge distance for this adapter.
// ---------------------------------------------------------------------------

/// Canonical distance (Cartesian Angstrom, via [`chematic_crystal::minimum_image`])
/// within which two sites are treated as "the same site" -- both for
/// `_atom_site_*` disorder-row grouping (see `cif.rs`'s
/// `parse_cif_periodic_structure_with_options`) and for symmetry-expansion
/// special-position dedup ([`expand_sites`]).
///
/// Previously the disorder-grouping check used a *fractional*-coordinate
/// tolerance (`1e-4`, dropped in favor of this constant) compared
/// component-wise -- not rotation/skew-invariant, and expressed in
/// different units than what expansion needs (expansion must compare
/// distances between points related by an arbitrary rotation matrix, where
/// a per-component fractional check is not meaningful). Unified here to one
/// Cartesian-Angstrom value instead of letting two similar-but-different
/// tolerances drift apart.
///
/// Value: `1e-4` fractional on a representative ~10 Angstrom cell edge is
/// about `1e-3` Angstrom, so this keeps the old disorder-grouping behavior
/// at typical cell sizes. It sits about three orders of magnitude below
/// real interatomic distances (~1 Angstrom), so it is loose enough to
/// absorb the coordinate-precision CIF authors actually write (4-5 decimal
/// fractional places) while remaining far too tight to merge two genuinely
/// distinct sites.
pub(crate) const SITE_MERGE_TOLERANCE_ANGSTROM: f64 = 1e-3;

/// `true` if `a` and `b` are the same site under [`SITE_MERGE_TOLERANCE_ANGSTROM`].
///
/// ponytail: O(1) per call, but every caller currently uses this inside an
/// O(n) or O(n^2) scan over already-small (tens-to-low-hundreds row) CIF
/// site lists -- fine at that scale, not built for bulk/corpus-scale input
/// (see [`MAX_EXPANSION_PRODUCT`] for where this is capped on the expansion
/// side).
pub(crate) fn sites_within_tolerance(
    lattice: &Lattice,
    a: FractionalCoord,
    b: FractionalCoord,
) -> bool {
    minimum_image(lattice, a, b).distance <= SITE_MERGE_TOLERANCE_ANGSTROM
}

/// Cap on `operation_count * asymmetric_site_count` before [`expand_sites`]
/// refuses with [`CifSymmetryError::ExpansionTooLarge`].
///
/// This bounds the O(n^2) dedup scan in [`expand_sites`] (every newly
/// expanded site is checked against every previously accepted output site),
/// not just allocation size -- at the cap, dedup is already up to
/// `MAX_EXPANSION_PRODUCT^2` = 16,000,000 [`chematic_crystal::minimum_image`]
/// calls in the worst case. `4000` comfortably covers realistic small-molecule
/// and inorganic asymmetric units (real CIFs rarely exceed a few dozen
/// distinct sites times the largest published operation counts, ~192 for
/// cubic groups), while still failing closed well before a pathological or
/// adversarial CIF could hang this function.
const MAX_EXPANSION_PRODUCT: usize = 4000;

// ---------------------------------------------------------------------------
// Rational
// ---------------------------------------------------------------------------

/// A small internal exact-fraction type. Always stored in lowest terms with
/// a strictly positive denominator, so `PartialEq` (derived) is exact value
/// equality, not representation equality.
///
/// No external crate: this workspace has no general-purpose rational-number
/// type (`chematic-cip::rational::RationalAtomicNumber` is CIP-specific and
/// not reusable here), and CIF symop translations only ever need small
/// exact fractions (halves, thirds, quarters, sixths).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Rational {
    num: i64,
    den: i64,
}

fn gcd_u64(a: u64, b: u64) -> u64 {
    let (mut a, mut b) = (a, b);
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

impl Rational {
    pub(crate) const ZERO: Rational = Rational { num: 0, den: 1 };

    /// Construct and reduce to lowest terms. Rejects `den == 0`.
    fn new(num: i64, den: i64, text: &str) -> Result<Self, CifSymmetryError> {
        if den == 0 {
            return Err(CifSymmetryError::ZeroDenominator {
                text: text.to_string(),
            });
        }
        // Normalize sign onto the numerator.
        let (num, den) = if den < 0 {
            let num = num
                .checked_neg()
                .ok_or_else(|| overflow(text, "translation numerator negation"))?;
            let den = den
                .checked_neg()
                .ok_or_else(|| overflow(text, "translation denominator negation"))?;
            (num, den)
        } else {
            (num, den)
        };
        let g = gcd_u64(num.unsigned_abs(), den.unsigned_abs());
        let g = if g == 0 { 1 } else { g as i64 };
        Ok(Rational {
            num: num / g,
            den: den / g,
        })
    }

    fn checked_add(self, other: Rational, text: &str) -> Result<Rational, CifSymmetryError> {
        let a = self
            .num
            .checked_mul(other.den)
            .ok_or_else(|| overflow(text, "translation sum"))?;
        let b = other
            .num
            .checked_mul(self.den)
            .ok_or_else(|| overflow(text, "translation sum"))?;
        let num = a
            .checked_add(b)
            .ok_or_else(|| overflow(text, "translation sum"))?;
        let den = self
            .den
            .checked_mul(other.den)
            .ok_or_else(|| overflow(text, "translation sum"))?;
        Rational::new(num, den, text)
    }

    fn checked_neg(self, text: &str) -> Result<Rational, CifSymmetryError> {
        Ok(Rational {
            num: self
                .num
                .checked_neg()
                .ok_or_else(|| overflow(text, "translation negation"))?,
            den: self.den,
        })
    }

    fn to_f64(self) -> f64 {
        self.num as f64 / self.den as f64
    }

    /// Reduce into `[0, 1)` -- used only for identity/duplicate/conflict
    /// comparisons (two operations whose translations differ by a whole
    /// lattice vector produce identical expanded positions after wrapping,
    /// so they must compare equal here even though their raw `Rational`
    /// values differ).
    fn mod1(self) -> Rational {
        Rational {
            num: self.num.rem_euclid(self.den),
            den: self.den,
        }
    }
}

fn overflow(text: &str, what: &str) -> CifSymmetryError {
    CifSymmetryError::MalformedSymmetryOperation {
        text: text.to_string(),
        reason: format!("arithmetic overflow while parsing {what}"),
    }
}

// ---------------------------------------------------------------------------
// CifSymmetryOperation
// ---------------------------------------------------------------------------

/// One parsed CIF symmetry operation: `X' = rotation * X + translation`.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CifSymmetryOperation {
    pub(crate) rotation: [[i8; 3]; 3],
    pub(crate) translation: [Rational; 3],
    /// The `_space_group_symop.id`/`_symmetry_equiv_pos_site_id` value (or
    /// alias), if the source loop had an id column. Carried through
    /// parsing so an id-column CIF round-trips its column shape correctly
    /// (see [`scan_operation_sources`]); not otherwise consulted by
    /// identity/duplicate/dedup logic, which all key on parsed content.
    pub(crate) source_id: Option<String>,
    pub(crate) source_text: String,
}

const IDENTITY_ROTATION: [[i8; 3]; 3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

fn is_identity(op: &CifSymmetryOperation) -> bool {
    op.rotation == IDENTITY_ROTATION && op.translation.iter().all(|t| t.mod1() == Rational::ZERO)
}

fn operations_equal_content(a: &CifSymmetryOperation, b: &CifSymmetryOperation) -> bool {
    a.rotation == b.rotation && (0..3).all(|i| a.translation[i].mod1() == b.translation[i].mod1())
}

fn operations_equivalent_as_sets(a: &[CifSymmetryOperation], b: &[CifSymmetryOperation]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    // Both sides are already known duplicate-free (checked before this is
    // called), so simple existential matching is enough -- no multiset
    // bookkeeping needed.
    a.iter()
        .all(|x| b.iter().any(|y| operations_equal_content(x, y)))
}

/// Determinant of a 3x3 matrix with entries in `{-1, 0, 1}` (guaranteed by
/// [`parse_coordinate_row`]'s per-variable-used-once rule). Plain `i64`
/// arithmetic, not checked: with every entry bounded to `{-1, 0, 1}`, every
/// intermediate product/sum is bounded by a small constant, so overflow is
/// not reachable regardless of input -- checked arithmetic would be
/// dead-code discipline theater here, unlike the [`Rational`] arithmetic
/// above (which parses attacker-controlled magnitudes from CIF text).
fn determinant3(m: &[[i8; 3]; 3]) -> i64 {
    let a = m.map(|row| row.map(i64::from));
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Every way parsing or expanding a CIF's explicit symmetry-operation list
/// can fail. Nested inside [`crate::cif::CifPeriodicError::Symmetry`].
///
/// Fail-closed throughout: none of these are ever downgraded to a silent
/// fallback (e.g. a malformed operation list never silently degrades to
/// "treat as P1" or "skip that operation") -- see this crate's `cif.rs`
/// module docs and the workspace-wide CLAUDE.md/PR-policy discipline this
/// mirrors.
#[derive(Debug, Clone, PartialEq)]
pub enum CifSymmetryError {
    /// The operation-expression text could not be tokenized at all (wrong
    /// comma count, dangling/doubled `+`/`-`, empty expression, integer
    /// overflow while parsing a translation constant, ...).
    MalformedSymmetryOperation { text: String, reason: String },
    /// The text tokenized but uses syntax this grammar deliberately does
    /// not support: an unknown variable (anything but x/y/z), a
    /// non-±1 coefficient (a variable used more than once in one
    /// coordinate expression), parentheses, a multiplication symbol, or a
    /// non-finite result.
    UnsupportedSymmetryExpression { text: String, reason: String },
    /// A translation fraction's denominator was written as `0`.
    ZeroDenominator { text: String },
    /// The parsed 3x3 rotation matrix's determinant is not exactly `+1` or
    /// `-1` (not a valid crystallographic rotation/rotoinversion).
    InvalidRotationMatrix { text: String, determinant: i64 },
    /// The same operation (rotation + translation, compared modulo integer
    /// translation) appears more than once in the resolved operation list.
    DuplicateSymmetryOperation { text: String, duplicate_of: String },
    /// The resolved operation list contains no identity operation. Per
    /// IUCr convention the identity must be explicitly listed; this
    /// implementation does not silently add one. A list with exactly one
    /// operation that *is* the identity is not an error (trivial,
    /// P1-equivalent expansion) -- this only fires when identity is truly
    /// absent.
    MissingIdentityOperation,
    /// Two different operation-list tag aliases (e.g.
    /// `_space_group_symop_operation_xyz` and `_symmetry_equiv_pos_as_xyz`)
    /// are both present in the same CIF and parse to genuinely different
    /// operation sets.
    ConflictingSymmetryOperationLists,
    /// Two different asymmetric-unit sites expanded (under some pair of
    /// operations) to the same position but with different species
    /// compositions -- not a legitimate special-position merge.
    SymmetrySiteCollision {
        label_a: Option<String>,
        label_b: Option<String>,
    },
    /// `operation_count * asymmetric_site_count` exceeds
    /// [`MAX_EXPANSION_PRODUCT`] -- refused before attempting the
    /// allocation/dedup rather than risking unbounded work on a
    /// pathological or adversarial CIF.
    ExpansionTooLarge {
        operation_count: usize,
        asymmetric_site_count: usize,
        limit: usize,
    },
}

impl core::fmt::Display for CifSymmetryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::MalformedSymmetryOperation { text, reason } => {
                write!(f, "malformed symmetry operation '{text}': {reason}")
            }
            Self::UnsupportedSymmetryExpression { text, reason } => {
                write!(f, "unsupported symmetry expression '{text}': {reason}")
            }
            Self::ZeroDenominator { text } => {
                write!(f, "zero denominator in symmetry operation '{text}'")
            }
            Self::InvalidRotationMatrix { text, determinant } => write!(
                f,
                "symmetry operation '{text}' has rotation-matrix determinant {determinant} \
                 (must be exactly +1 or -1)"
            ),
            Self::DuplicateSymmetryOperation { text, duplicate_of } => write!(
                f,
                "symmetry operation '{text}' duplicates '{duplicate_of}' (same rotation and \
                 translation modulo integer translation)"
            ),
            Self::MissingIdentityOperation => write!(
                f,
                "no identity operation ('x,y,z') found in the CIF's explicit symmetry-operation \
                 list -- IUCr convention requires it to be listed explicitly, so it is not added \
                 silently"
            ),
            Self::ConflictingSymmetryOperationLists => write!(
                f,
                "multiple symmetry-operation tag aliases are present in this CIF and parse to \
                 genuinely different operation sets"
            ),
            Self::SymmetrySiteCollision { label_a, label_b } => write!(
                f,
                "symmetry expansion produced a position shared by two sites with different \
                 species compositions (labels {label_a:?} and {label_b:?})"
            ),
            Self::ExpansionTooLarge {
                operation_count,
                asymmetric_site_count,
                limit,
            } => write!(
                f,
                "symmetry expansion would produce up to {operation_count} * \
                 {asymmetric_site_count} = {} candidate sites, exceeding the {limit} cap",
                operation_count.saturating_mul(*asymmetric_site_count)
            ),
        }
    }
}

impl std::error::Error for CifSymmetryError {}

// ---------------------------------------------------------------------------
// Tag scanning
// ---------------------------------------------------------------------------

const OP_TAGS: &[&str] = &[
    "_space_group_symop.operation_xyz",
    "_space_group_symop_operation_xyz",
    "_symmetry_equiv.pos_as_xyz",
    "_symmetry_equiv_pos_as_xyz",
];

const ID_TAGS: &[&str] = &[
    "_space_group_symop.id",
    "_space_group_symop_id",
    "_symmetry_equiv.pos_site_id",
    "_symmetry_equiv_pos_site_id",
];

/// One raw `(id, expression text)` row, before parsing.
type RawOp = (Option<String>, String);

/// Find every place in the token stream that declares an operation list --
/// a `loop_` block whose headers include one of [`OP_TAGS`], or a
/// standalone (non-loop) `_tag value` data item using one of [`OP_TAGS`].
/// Each match is one "source"; multiple sources exist when a CIF declares
/// the same conceptual list under more than one tag alias (see
/// [`resolve_symmetry_operations`] for how those are cross-checked).
///
/// Pure tokenizing -- does not parse/validate expression text, so this is
/// also reused (via `.first().map(Vec::len)`) to report an `operation_count`
/// when the caller opted out of expansion entirely
/// (`CifPeriodicParseOptions::expand_explicit_symmetry = false`).
pub(crate) fn scan_operation_sources(tokens: &[String]) -> Vec<Vec<RawOp>> {
    let mut sources: Vec<Vec<RawOp>> = Vec::new();
    let mut i = 0;
    while i < tokens.len() {
        if tokens[i] == "loop_" {
            let mut j = i + 1;
            let mut headers: Vec<String> = Vec::new();
            while j < tokens.len() && tokens[j].starts_with('_') {
                headers.push(tokens[j].to_ascii_lowercase());
                j += 1;
            }
            let ncols = headers.len();
            let op_col = headers.iter().position(|h| OP_TAGS.contains(&h.as_str()));
            if let (Some(op_col), true) = (op_col, ncols > 0) {
                let id_col = headers.iter().position(|h| ID_TAGS.contains(&h.as_str()));
                let mut k = j;
                let mut rows: Vec<RawOp> = Vec::new();
                while k + ncols <= tokens.len() {
                    let row = &tokens[k..k + ncols];
                    if row[0] == "loop_" || row[0].starts_with("data_") || row[0].starts_with('_') {
                        break;
                    }
                    rows.push((id_col.map(|c| row[c].clone()), row[op_col].clone()));
                    k += ncols;
                }
                sources.push(rows);
                i = k;
                continue;
            }
            i = j;
            continue;
        }
        if tokens[i].starts_with('_') {
            let tag = tokens[i].to_ascii_lowercase();
            if OP_TAGS.contains(&tag.as_str()) && i + 1 < tokens.len() {
                // Standalone single data item -- one operation, no id column
                // (see module docs: this is the "genuine P1 CIF declares its
                // one identity operation explicitly" shape).
                sources.push(vec![(None, tokens[i + 1].clone())]);
                i += 2;
                continue;
            }
        }
        i += 1;
    }
    sources
}

// ---------------------------------------------------------------------------
// Expression parsing
// ---------------------------------------------------------------------------

enum TermKind {
    Var(usize),
    Const(Rational),
}

/// Split `expr` (already lowercased and whitespace-stripped) into
/// `(sign, chunk_text)` pairs at each top-level `+`/`-`. Rejects a leading
/// double sign, a doubled/adjacent operator, and a dangling trailing
/// operator -- all of which would otherwise silently produce an empty
/// chunk.
fn split_signed_terms(
    expr: &str,
    source_text: &str,
) -> Result<Vec<(i8, String)>, CifSymmetryError> {
    let chars: Vec<char> = expr.chars().collect();
    if chars.is_empty() {
        return Err(CifSymmetryError::MalformedSymmetryOperation {
            text: source_text.to_string(),
            reason: "empty coordinate expression".to_string(),
        });
    }
    let mut terms = Vec::new();
    let mut sign: i8 = 1;
    let mut current = String::new();
    for idx in 0..chars.len() {
        let c = chars[idx];
        if c == '+' || c == '-' {
            if idx > 0 && (chars[idx - 1] == '+' || chars[idx - 1] == '-') {
                return Err(CifSymmetryError::MalformedSymmetryOperation {
                    text: source_text.to_string(),
                    reason: "two operators in a row".to_string(),
                });
            }
            if idx == chars.len() - 1 {
                return Err(CifSymmetryError::MalformedSymmetryOperation {
                    text: source_text.to_string(),
                    reason: "dangling operator at end of expression".to_string(),
                });
            }
            if idx > 0 {
                if current.is_empty() {
                    return Err(CifSymmetryError::MalformedSymmetryOperation {
                        text: source_text.to_string(),
                        reason: "empty term".to_string(),
                    });
                }
                terms.push((sign, current.clone()));
                current.clear();
            }
            sign = if c == '+' { 1 } else { -1 };
        } else {
            current.push(c);
        }
    }
    if current.is_empty() {
        return Err(CifSymmetryError::MalformedSymmetryOperation {
            text: source_text.to_string(),
            reason: "empty term".to_string(),
        });
    }
    terms.push((sign, current));
    Ok(terms)
}

/// Classify one sign-stripped term chunk: a bare `x`/`y`/`z` (variable,
/// coefficient always exactly ±1 via the term's sign) or an integer/rational
/// constant (`"2"`, `"1/2"`). Anything else -- an unknown letter,
/// parentheses, a multiplication symbol, more than one `/`, a decimal point
/// -- is rejected.
fn classify_term(text: &str, source_text: &str) -> Result<TermKind, CifSymmetryError> {
    if text.chars().count() == 1 {
        let c = text.chars().next().expect("length checked above");
        if c.is_ascii_alphabetic() {
            return match c {
                'x' => Ok(TermKind::Var(0)),
                'y' => Ok(TermKind::Var(1)),
                'z' => Ok(TermKind::Var(2)),
                other => Err(CifSymmetryError::UnsupportedSymmetryExpression {
                    text: source_text.to_string(),
                    reason: format!("unknown variable '{other}' (only x, y, z are supported)"),
                }),
            };
        }
    }
    if text.chars().all(|c| c.is_ascii_digit() || c == '/') {
        let parts: Vec<&str> = text.split('/').collect();
        return match parts.as_slice() {
            [n] if !n.is_empty() => {
                let num =
                    n.parse::<i64>()
                        .map_err(|_| CifSymmetryError::MalformedSymmetryOperation {
                            text: source_text.to_string(),
                            reason: format!("integer constant '{n}' out of range"),
                        })?;
                Ok(TermKind::Const(Rational::new(num, 1, source_text)?))
            }
            [n, d] if !n.is_empty() && !d.is_empty() => {
                let num =
                    n.parse::<i64>()
                        .map_err(|_| CifSymmetryError::MalformedSymmetryOperation {
                            text: source_text.to_string(),
                            reason: format!("fraction numerator '{n}' out of range"),
                        })?;
                let den =
                    d.parse::<i64>()
                        .map_err(|_| CifSymmetryError::MalformedSymmetryOperation {
                            text: source_text.to_string(),
                            reason: format!("fraction denominator '{d}' out of range"),
                        })?;
                Ok(TermKind::Const(Rational::new(num, den, source_text)?))
            }
            _ => Err(CifSymmetryError::MalformedSymmetryOperation {
                text: source_text.to_string(),
                reason: format!("'{text}' is not a valid integer or single fraction"),
            }),
        };
    }
    Err(CifSymmetryError::UnsupportedSymmetryExpression {
        text: source_text.to_string(),
        reason: format!(
            "'{text}' is neither a bare x/y/z variable nor an integer/fractional constant \
             (parentheses, coefficients other than ±1, and multiplication are not supported)"
        ),
    })
}

/// Parse one coordinate expression (e.g. `"x-y+1/2"`) into a rotation-matrix
/// row (coefficients for x, y, z) and a translation constant.
fn parse_coordinate_row(
    expr: &str,
    source_text: &str,
) -> Result<([i8; 3], Rational), CifSymmetryError> {
    let terms = split_signed_terms(expr, source_text)?;
    let mut row = [0i8; 3];
    let mut seen = [false; 3];
    let mut translation = Rational::ZERO;
    for (sign, text) in terms {
        match classify_term(&text, source_text)? {
            TermKind::Var(idx) => {
                if seen[idx] {
                    return Err(CifSymmetryError::UnsupportedSymmetryExpression {
                        text: source_text.to_string(),
                        reason: format!(
                            "variable '{}' appears more than once (only a ±1 coefficient is \
                             supported)",
                            ['x', 'y', 'z'][idx]
                        ),
                    });
                }
                seen[idx] = true;
                row[idx] = sign;
            }
            TermKind::Const(r) => {
                let signed = if sign < 0 {
                    r.checked_neg(source_text)?
                } else {
                    r
                };
                translation = translation.checked_add(signed, source_text)?;
            }
        }
    }
    Ok((row, translation))
}

/// Parse one full operation expression (`"x,y,z"`, `"-x,-y,-z"`, ...) into a
/// [`CifSymmetryOperation`]. Requires exactly 3 comma-separated coordinate
/// expressions and a rotation-matrix determinant of exactly `+1` or `-1`.
fn parse_operation(
    raw_text: &str,
    id: Option<String>,
) -> Result<CifSymmetryOperation, CifSymmetryError> {
    let normalized: String = raw_text
        .chars()
        .filter(|c| !c.is_whitespace())
        .flat_map(|c| c.to_lowercase())
        .collect();
    let parts: Vec<&str> = normalized.split(',').collect();
    if parts.len() != 3 {
        return Err(CifSymmetryError::MalformedSymmetryOperation {
            text: raw_text.to_string(),
            reason: format!(
                "expected exactly 3 comma-separated coordinate expressions, found {}",
                parts.len()
            ),
        });
    }

    let mut rotation = [[0i8; 3]; 3];
    let mut translation = [Rational::ZERO; 3];
    for (i, part) in parts.iter().enumerate() {
        let (row, t) = parse_coordinate_row(part, raw_text)?;
        rotation[i] = row;
        translation[i] = t;
    }

    let det = determinant3(&rotation);
    if det != 1 && det != -1 {
        return Err(CifSymmetryError::InvalidRotationMatrix {
            text: raw_text.to_string(),
            determinant: det,
        });
    }

    Ok(CifSymmetryOperation {
        rotation,
        translation,
        source_id: id,
        source_text: raw_text.to_string(),
    })
}

// ---------------------------------------------------------------------------
// Resolution: parse every source, cross-check aliases, validate the result
// ---------------------------------------------------------------------------

fn check_no_duplicates(ops: &[CifSymmetryOperation]) -> Result<(), CifSymmetryError> {
    for i in 0..ops.len() {
        for j in (i + 1)..ops.len() {
            if operations_equal_content(&ops[i], &ops[j]) {
                return Err(CifSymmetryError::DuplicateSymmetryOperation {
                    text: ops[j].source_text.clone(),
                    duplicate_of: ops[i].source_text.clone(),
                });
            }
        }
    }
    Ok(())
}

/// Scan, parse, cross-check, and validate a CIF's explicit symmetry
/// operation list.
///
/// Returns `Ok(None)` if no operation-list tag alias is present anywhere in
/// the token stream (nothing to expand from). Returns `Ok(Some(ops))` with
/// the file's original declared order otherwise (identity-first reordering
/// for expansion output happens later, in [`expand_sites`]).
///
/// If more than one tag alias declares an operation list, each is parsed
/// and duplicate-checked independently (order-independent: a duplicate
/// inside *either* list is caught regardless of which alias happens to
/// appear first in the file), then compared as order-independent sets --
/// genuinely different content is [`CifSymmetryError::ConflictingSymmetryOperationLists`],
/// equivalent content silently collapses to one (the first-encountered, by
/// file order) list.
pub(crate) fn resolve_symmetry_operations(
    tokens: &[String],
) -> Result<Option<Vec<CifSymmetryOperation>>, CifSymmetryError> {
    let raw_sources = scan_operation_sources(tokens);
    if raw_sources.is_empty() {
        return Ok(None);
    }

    let mut parsed_sources: Vec<Vec<CifSymmetryOperation>> = Vec::with_capacity(raw_sources.len());
    for source in &raw_sources {
        let mut ops = Vec::with_capacity(source.len());
        for (id, text) in source {
            ops.push(parse_operation(text, id.clone())?);
        }
        check_no_duplicates(&ops)?;
        parsed_sources.push(ops);
    }

    for other in &parsed_sources[1..] {
        if !operations_equivalent_as_sets(&parsed_sources[0], other) {
            return Err(CifSymmetryError::ConflictingSymmetryOperationLists);
        }
    }

    let canonical = parsed_sources
        .into_iter()
        .next()
        .expect("checked non-empty above");
    if !canonical.iter().any(is_identity) {
        return Err(CifSymmetryError::MissingIdentityOperation);
    }
    Ok(Some(canonical))
}

// ---------------------------------------------------------------------------
// Expansion
// ---------------------------------------------------------------------------

fn apply_operation(op: &CifSymmetryOperation, frac: FractionalCoord) -> FractionalCoord {
    let f = frac.0;
    let mut out = [0.0f64; 3];
    for (i, out_i) in out.iter_mut().enumerate() {
        let mut v = op.translation[i].to_f64();
        for (rot_ij, f_j) in op.rotation[i].iter().zip(f.iter()) {
            v += f64::from(*rot_ij) * f_j;
        }
        *out_i = v;
    }
    FractionalCoord::new(out).wrapped()
}

fn species_match(a: &[SiteSpecies], b: &[SiteSpecies]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    let key =
        |s: &SiteSpecies| -> (&'static str, f64) { (s.element.symbol(), s.occupancy.value()) };
    let mut a_sorted: Vec<(&str, f64)> = a.iter().map(key).collect();
    let mut b_sorted: Vec<(&str, f64)> = b.iter().map(key).collect();
    a_sorted.sort_by(|x, y| x.0.cmp(y.0).then(x.1.total_cmp(&y.1)));
    b_sorted.sort_by(|x, y| x.0.cmp(y.0).then(x.1.total_cmp(&y.1)));
    a_sorted
        .iter()
        .zip(b_sorted.iter())
        .all(|(x, y)| x.0 == y.0 && (x.1 - y.1).abs() < 1e-9)
}

/// Expand `sites` (the CIF's original asymmetric-unit order) under
/// `operations` (identity presence already validated by
/// [`resolve_symmetry_operations`]).
///
/// Output order: asymmetric-unit site order outermost, then within each
/// site the identity operation's image first (regardless of where identity
/// sits in `operations`' declared order), then the remaining operations in
/// their original declared order. A candidate that lands (within
/// [`SITE_MERGE_TOLERANCE_ANGSTROM`], via [`chematic_crystal::minimum_image`])
/// on an already-accepted output site is dropped if species match (special-
/// position merge -- occupancy is never summed/duplicated) or rejected with
/// [`CifSymmetryError::SymmetrySiteCollision`] if species differ.
///
/// A non-identity image's label is the source site's label (if any) with
/// `@sym{N}` appended, where `N` is the 1-based position of the producing
/// operation in `operations`' *original declared* order (stable regardless
/// of the identity-first output reordering). A site with no source label
/// gets no label on any of its images (never synthesized).
pub(crate) fn expand_sites(
    sites: &[PeriodicSite],
    operations: &[CifSymmetryOperation],
    lattice: &Lattice,
) -> Result<Vec<PeriodicSite>, CifPeriodicError> {
    let product = operations
        .len()
        .checked_mul(sites.len())
        .filter(|p| *p <= MAX_EXPANSION_PRODUCT);
    if product.is_none() {
        return Err(CifPeriodicError::Symmetry(
            CifSymmetryError::ExpansionTooLarge {
                operation_count: operations.len(),
                asymmetric_site_count: sites.len(),
                limit: MAX_EXPANSION_PRODUCT,
            },
        ));
    }

    let identity_idx = operations
        .iter()
        .position(is_identity)
        .expect("resolve_symmetry_operations already validated identity presence");

    let mut ordered_ops: Vec<(usize, &CifSymmetryOperation)> = Vec::with_capacity(operations.len());
    ordered_ops.push((identity_idx + 1, &operations[identity_idx]));
    for (idx, op) in operations.iter().enumerate() {
        if idx != identity_idx {
            ordered_ops.push((idx + 1, op));
        }
    }

    let mut output: Vec<PeriodicSite> = Vec::new();
    for site in sites {
        for (op_number, op) in &ordered_ops {
            let expanded_frac = apply_operation(op, site.fractional);
            let collision = output.iter().position(|existing| {
                sites_within_tolerance(lattice, existing.fractional, expanded_frac)
            });
            match collision {
                Some(existing_idx) => {
                    if !species_match(&output[existing_idx].species, &site.species) {
                        return Err(CifPeriodicError::Symmetry(
                            CifSymmetryError::SymmetrySiteCollision {
                                label_a: output[existing_idx].label.clone(),
                                label_b: site.label.clone(),
                            },
                        ));
                    }
                    // Legitimate special-position merge: keep the
                    // already-accepted site, do not duplicate/sum.
                }
                None => {
                    let label = if *op_number == identity_idx + 1 {
                        site.label.clone()
                    } else {
                        site.label.as_ref().map(|l| format!("{l}@sym{op_number}"))
                    };
                    output.push(PeriodicSite::new(
                        site.species.clone(),
                        expanded_frac,
                        label,
                    )?);
                }
            }
        }
    }
    Ok(output)
}

// ---------------------------------------------------------------------------
// Tests -- parser/expansion internals in isolation (white-box: this module's
// pub(crate) items are directly reachable). End-to-end CIF-text-level
// fixtures live in `cif.rs`'s `crystal_adapter::tests`.
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::Element;
    use chematic_crystal::{Occupancy, SiteSpecies};

    fn op(text: &str) -> CifSymmetryOperation {
        parse_operation(text, None).unwrap_or_else(|e| panic!("failed to parse '{text}': {e}"))
    }

    fn site(el: Element, frac: [f64; 3], label: Option<&str>) -> PeriodicSite {
        PeriodicSite::new(
            vec![SiteSpecies::full(el)],
            FractionalCoord::new(frac),
            label.map(str::to_string),
        )
        .unwrap()
    }

    // -- Rational -----------------------------------------------------

    #[test]
    fn rational_reduces_to_lowest_terms() {
        let r = Rational::new(2, 4, "2/4").unwrap();
        assert_eq!(r, Rational::new(1, 2, "1/2").unwrap());
    }

    #[test]
    fn rational_normalizes_sign_onto_numerator() {
        let r = Rational::new(1, -2, "1/-2").unwrap();
        assert_eq!(r, Rational::new(-1, 2, "-1/2").unwrap());
    }

    #[test]
    fn rational_zero_denominator_is_rejected() {
        assert_eq!(
            Rational::new(1, 0, "1/0"),
            Err(CifSymmetryError::ZeroDenominator {
                text: "1/0".to_string()
            })
        );
    }

    #[test]
    fn rational_mod1_wraps_negative_and_over_one() {
        let r = Rational::new(-1, 3, "-1/3").unwrap();
        assert_eq!(r.mod1(), Rational::new(2, 3, "2/3").unwrap());
        let r2 = Rational::new(4, 3, "4/3").unwrap();
        assert_eq!(r2.mod1(), Rational::new(1, 3, "1/3").unwrap());
    }

    // -- parse_operation: accepted grammar ----------------------------

    #[test]
    fn parses_identity() {
        let o = op("x,y,z");
        assert_eq!(o.rotation, IDENTITY_ROTATION);
        assert!(o.translation.iter().all(|t| *t == Rational::ZERO));
    }

    #[test]
    fn parses_inversion_case_and_whitespace_insensitive() {
        let a = op("-x,-y,-z");
        let b = op(" -X , -Y , -Z ");
        assert_eq!(a.rotation, b.rotation);
        assert_eq!(a.translation, b.translation);
        assert_eq!(a.rotation, [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]);
    }

    #[test]
    fn parses_half_translation_and_leading_sign() {
        let o = op("x+1/2,y,z");
        assert_eq!(o.rotation, IDENTITY_ROTATION);
        assert_eq!(o.translation[0], Rational::new(1, 2, "").unwrap());
        assert_eq!(o.translation[1], Rational::ZERO);
    }

    #[test]
    fn parses_leading_fraction_before_variable() {
        let o = op("1/2-x,-y,1/2+z");
        assert_eq!(o.rotation[0], [-1, 0, 0]);
        assert_eq!(o.translation[0], Rational::new(1, 2, "").unwrap());
    }

    #[test]
    fn parses_combined_variable_terms_and_rational_translation() {
        let o = op("-y+x,-y,1/3+z");
        assert_eq!(o.rotation, [[1, -1, 0], [0, -1, 0], [0, 0, 1]]);
        assert_eq!(o.translation[2], Rational::new(1, 3, "").unwrap());
        assert_eq!(determinant3(&o.rotation), -1);
    }

    #[test]
    fn parses_x_minus_y_two_thirds() {
        let o = op("x-y,x,z+2/3");
        assert_eq!(o.rotation, [[1, -1, 0], [1, 0, 0], [0, 0, 1]]);
        assert_eq!(o.translation[2], Rational::new(2, 3, "").unwrap());
    }

    // -- parse_operation: rejections ----------------------------------

    #[test]
    fn rejects_wrong_comma_count() {
        assert!(matches!(
            parse_operation("x,y", None),
            Err(CifSymmetryError::MalformedSymmetryOperation { .. })
        ));
        assert!(matches!(
            parse_operation("x,y,z,w", None),
            Err(CifSymmetryError::MalformedSymmetryOperation { .. })
        ));
    }

    #[test]
    fn rejects_unknown_variable() {
        assert!(matches!(
            parse_operation("x,y,w", None),
            Err(CifSymmetryError::UnsupportedSymmetryExpression { .. })
        ));
    }

    #[test]
    fn rejects_malformed_double_slash_fraction() {
        assert!(matches!(
            parse_operation("x,y,1/2/3", None),
            Err(CifSymmetryError::MalformedSymmetryOperation { .. })
        ));
    }

    #[test]
    fn rejects_zero_denominator() {
        assert_eq!(
            parse_operation("x,y,1/0", None),
            Err(CifSymmetryError::ZeroDenominator {
                text: "x,y,1/0".to_string()
            })
        );
    }

    #[test]
    fn rejects_repeated_variable_as_unsupported_coefficient() {
        assert!(matches!(
            parse_operation("x+x,y,z", None),
            Err(CifSymmetryError::UnsupportedSymmetryExpression { .. })
        ));
    }

    #[test]
    fn rejects_doubled_operator() {
        assert!(matches!(
            parse_operation("x+-y,y,z", None),
            Err(CifSymmetryError::MalformedSymmetryOperation { .. })
        ));
    }

    #[test]
    fn rejects_dangling_operator() {
        assert!(matches!(
            parse_operation("x+,y,z", None),
            Err(CifSymmetryError::MalformedSymmetryOperation { .. })
        ));
    }

    #[test]
    fn rejects_parentheses_and_multiplication() {
        assert!(matches!(
            parse_operation("(x),y,z", None),
            Err(CifSymmetryError::UnsupportedSymmetryExpression { .. })
        ));
        assert!(matches!(
            parse_operation("2*x,y,z", None),
            Err(CifSymmetryError::UnsupportedSymmetryExpression { .. })
        ));
    }

    #[test]
    fn rejects_singular_rotation_matrix() {
        // Two rows both encode "x" (row0 == row1): determinant 0.
        match parse_operation("x,x,z", None) {
            Err(CifSymmetryError::InvalidRotationMatrix { determinant, .. }) => {
                assert_eq!(determinant, 0);
            }
            other => panic!("expected InvalidRotationMatrix, got {other:?}"),
        }
    }

    // -- resolve_symmetry_operations ------------------------------------

    #[test]
    fn missing_identity_is_rejected() {
        let tokens: Vec<String> = "loop_ _symmetry_equiv_pos_as_xyz -x,-y,-z"
            .split_whitespace()
            .map(str::to_string)
            .collect();
        assert_eq!(
            resolve_symmetry_operations(&tokens),
            Err(CifSymmetryError::MissingIdentityOperation)
        );
    }

    #[test]
    fn single_identity_operation_resolves_successfully() {
        let tokens: Vec<String> = vec![
            "_symmetry_equiv_pos_as_xyz".to_string(),
            "x,y,z".to_string(),
        ];
        let ops = resolve_symmetry_operations(&tokens).unwrap().unwrap();
        assert_eq!(ops.len(), 1);
        assert!(is_identity(&ops[0]));
    }

    #[test]
    fn duplicate_operation_within_one_source_is_rejected() {
        let tokens: Vec<String> = vec![
            "loop_".to_string(),
            "_symmetry_equiv_pos_as_xyz".to_string(),
            "x,y,z".to_string(),
            "-x,-y,-z".to_string(),
            "-x,-y,-z".to_string(),
        ];
        assert!(matches!(
            resolve_symmetry_operations(&tokens),
            Err(CifSymmetryError::DuplicateSymmetryOperation { .. })
        ));
    }

    #[test]
    fn no_operation_tag_present_returns_none() {
        let tokens: Vec<String> = vec!["_cell_length_a".to_string(), "5.0".to_string()];
        assert_eq!(resolve_symmetry_operations(&tokens), Ok(None));
    }

    #[test]
    fn conflicting_alias_lists_are_rejected() {
        let tokens: Vec<String> = vec![
            "_symmetry_equiv_pos_as_xyz".to_string(),
            "x,y,z".to_string(),
            "_space_group_symop_operation_xyz".to_string(),
            "-x,-y,-z".to_string(),
        ];
        // Two standalone single-operation sources: [x,y,z] vs [-x,-y,-z] --
        // both individually valid (each contains no identity issue on its
        // own duplicate check), but genuinely different content.
        assert_eq!(
            resolve_symmetry_operations(&tokens),
            Err(CifSymmetryError::ConflictingSymmetryOperationLists)
        );
    }

    #[test]
    fn equivalent_alias_lists_collapse_to_one_without_error() {
        let tokens: Vec<String> = vec![
            "loop_".to_string(),
            "_symmetry_equiv_pos_as_xyz".to_string(),
            "x,y,z".to_string(),
            "-x,-y,-z".to_string(),
            "loop_".to_string(),
            "_space_group_symop_operation_xyz".to_string(),
            "-x,-y,-z".to_string(),
            "x,y,z".to_string(),
        ];
        let ops = resolve_symmetry_operations(&tokens).unwrap().unwrap();
        assert_eq!(ops.len(), 2);
    }

    #[test]
    fn multi_column_loop_with_id_column_is_parsed() {
        let tokens: Vec<String> = vec![
            "loop_".to_string(),
            "_space_group_symop_id".to_string(),
            "_space_group_symop_operation_xyz".to_string(),
            "1".to_string(),
            "x,y,z".to_string(),
            "2".to_string(),
            "-x,-y,-z".to_string(),
        ];
        let sources = scan_operation_sources(&tokens);
        assert_eq!(sources.len(), 1);
        assert_eq!(sources[0].len(), 2);
        assert_eq!(sources[0][0], (Some("1".to_string()), "x,y,z".to_string()));
        assert_eq!(
            sources[0][1],
            (Some("2".to_string()), "-x,-y,-z".to_string())
        );

        // The id survives all the way into the parsed CifSymmetryOperation,
        // not just the raw scan.
        let ops = resolve_symmetry_operations(&tokens).unwrap().unwrap();
        assert_eq!(ops[0].source_id.as_deref(), Some("1"));
        assert_eq!(ops[1].source_id.as_deref(), Some("2"));
    }

    // -- expand_sites -----------------------------------------------------

    fn cubic() -> Lattice {
        Lattice::cubic(10.0).unwrap()
    }

    #[test]
    fn expand_special_position_at_origin_is_not_duplicated() {
        let ops = vec![op("x,y,z"), op("-x,-y,-z")];
        let sites = vec![site(Element::C, [0.0, 0.0, 0.0], Some("C1"))];
        let out = expand_sites(&sites, &ops, &cubic()).unwrap();
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].label.as_deref(), Some("C1"));
    }

    #[test]
    fn expand_general_position_produces_one_site_per_operation() {
        let ops = vec![op("x,y,z"), op("-x,-y,-z")];
        let sites = vec![site(Element::C, [0.1, 0.2, 0.3], Some("C1"))];
        let out = expand_sites(&sites, &ops, &cubic()).unwrap();
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].label.as_deref(), Some("C1"));
        assert_eq!(out[1].label.as_deref(), Some("C1@sym2"));
        let f = out[1].fractional.0;
        assert!((f[0] - 0.9).abs() < 1e-9);
        assert!((f[1] - 0.8).abs() < 1e-9);
        assert!((f[2] - 0.7).abs() < 1e-9);
    }

    #[test]
    fn identity_not_first_in_declared_order_still_produces_identity_image_first() {
        // Inversion declared BEFORE identity in the source list.
        let ops = vec![op("-x,-y,-z"), op("x,y,z")];
        let sites = vec![site(Element::C, [0.1, 0.2, 0.3], Some("C1"))];
        let out = expand_sites(&sites, &ops, &cubic()).unwrap();
        assert_eq!(out.len(), 2);
        // Identity's image is first in output and keeps the label
        // unchanged; the inversion op (declared list position 1) produced
        // the second output site, suffixed by its ORIGINAL declared index.
        assert_eq!(out[0].label.as_deref(), Some("C1"));
        let f0 = out[0].fractional.0;
        assert!(
            (f0[0] - 0.1).abs() < 1e-9 && (f0[1] - 0.2).abs() < 1e-9 && (f0[2] - 0.3).abs() < 1e-9
        );
        assert_eq!(out[1].label.as_deref(), Some("C1@sym1"));
    }

    #[test]
    fn expand_preserves_unlabeled_site_as_unlabeled_throughout() {
        let ops = vec![op("x,y,z"), op("-x,-y,-z")];
        let sites = vec![site(Element::C, [0.1, 0.2, 0.3], None)];
        let out = expand_sites(&sites, &ops, &cubic()).unwrap();
        assert_eq!(out.len(), 2);
        assert!(out.iter().all(|s| s.label.is_none()));
    }

    #[test]
    fn expand_wraps_translated_coordinate_into_unit_cell() {
        let ops = vec![op("x,y,z"), op("x+1/2,y,z")];
        let sites = vec![site(Element::C, [0.8, 0.1, 0.1], Some("C1"))];
        let out = expand_sites(&sites, &ops, &cubic()).unwrap();
        assert_eq!(out.len(), 2);
        let f = out[1].fractional.0;
        assert!(
            (f[0] - 0.3).abs() < 1e-9,
            "expected wrap to 0.3, got {}",
            f[0]
        );
    }

    #[test]
    fn expand_preserves_disorder_composition_across_every_image() {
        let ops = vec![op("x,y,z"), op("-x,-y,-z")];
        let disordered = PeriodicSite::new(
            vec![
                SiteSpecies {
                    element: Element::FE,
                    occupancy: Occupancy::new(0.6).unwrap(),
                },
                SiteSpecies {
                    element: Element::NI,
                    occupancy: Occupancy::new(0.4).unwrap(),
                },
            ],
            FractionalCoord::new([0.1, 0.2, 0.3]),
            Some("Fe1".to_string()),
        )
        .unwrap();
        let out = expand_sites(&[disordered], &ops, &cubic()).unwrap();
        assert_eq!(out.len(), 2);
        for s in &out {
            assert_eq!(s.species.len(), 2);
            let fe = s
                .species
                .iter()
                .find(|sp| sp.element == Element::FE)
                .unwrap();
            let ni = s
                .species
                .iter()
                .find(|sp| sp.element == Element::NI)
                .unwrap();
            assert!((fe.occupancy.value() - 0.6).abs() < 1e-9);
            assert!((ni.occupancy.value() - 0.4).abs() < 1e-9);
        }
    }

    #[test]
    fn expand_different_species_colliding_at_same_position_is_a_typed_error() {
        let ops = vec![op("x,y,z"), op("-x,-y,-z")];
        // Both sites sit exactly at the origin, an inversion-invariant
        // special position, but hold different elements -- must not
        // silently pick one.
        let sites = vec![
            site(Element::FE, [0.0, 0.0, 0.0], Some("Fe1")),
            site(Element::NI, [0.0, 0.0, 0.0], Some("Ni1")),
        ];
        match expand_sites(&sites, &ops, &cubic()) {
            Err(CifPeriodicError::Symmetry(CifSymmetryError::SymmetrySiteCollision { .. })) => {}
            other => panic!("expected SymmetrySiteCollision, got {other:?}"),
        }
    }

    #[test]
    fn expansion_too_large_is_rejected_before_doing_the_work() {
        // operation_count (1) * asymmetric_site_count (MAX_EXPANSION_PRODUCT+1)
        // exceeds the cap -- the guard fires before the O(n^2) dedup scan
        // runs, so constructing the (cheap) oversized site list is fine;
        // this proves the actual `expand_sites` entry point rejects it, not
        // just the arithmetic in isolation.
        let ops: Vec<CifSymmetryOperation> = vec![op("x,y,z")];
        let too_many_sites: Vec<PeriodicSite> = (0..(MAX_EXPANSION_PRODUCT + 1))
            .map(|i| site(Element::C, [i as f64 / 1e6, 0.0, 0.0], None))
            .collect();
        match expand_sites(&too_many_sites, &ops, &cubic()) {
            Err(CifPeriodicError::Symmetry(CifSymmetryError::ExpansionTooLarge {
                operation_count,
                asymmetric_site_count,
                limit,
            })) => {
                assert_eq!(operation_count, 1);
                assert_eq!(asymmetric_site_count, MAX_EXPANSION_PRODUCT + 1);
                assert_eq!(limit, MAX_EXPANSION_PRODUCT);
            }
            other => panic!("expected ExpansionTooLarge, got {other:?}"),
        }

        // Sanity: the ordinary tiny case still succeeds (guard isn't
        // over-tight).
        let tiny_sites = vec![site(Element::C, [0.1, 0.2, 0.3], None)];
        assert!(expand_sites(&tiny_sites, &ops, &cubic()).is_ok());
    }
}
