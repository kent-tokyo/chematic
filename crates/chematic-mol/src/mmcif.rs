//! PDBx/mmCIF (macromolecular Crystallographic Information File) reader and
//! writer.
//!
//! ## Provenance
//!
//! Implemented independently from the public wwPDB PDBx/mmCIF dictionary
//! (`_atom_site`, `_cell`, `_symmetry` category definitions --
//! <https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Index/>) and the
//! general CIF/STAR syntax rules (IUCr CIF1.1 working spec). No source code,
//! comments, or tables were copied from Open Babel or any other tool; those
//! were consulted only informally as a behavioral sanity check, never as a
//! source of implementation text.
//!
//! ## Relationship to [`crate::cif`]
//!
//! `chematic-mol` already has a CIF reader/writer for **small-molecule /
//! crystallographic** CIF (`crate::cif`), which uses the older
//! `_atom_site_fract_x`-style (DDL1, underscore-joined) tag convention. This
//! module targets **macromolecular mmCIF** (DDL2, dot-joined
//! `_atom_site.Cartn_x`-style tags -- the convention every PDB-deposited
//! `.cif` file actually uses) and the residue/chain/model bookkeeping
//! (`_atom_site.label_asym_id`, `.auth_seq_id`, `.pdbx_PDB_model_num`, ...)
//! that small-molecule CIF has no equivalent of.
//!
//! Both formats share the same underlying STAR tokenizing rules (quoting,
//! comments, `loop_` blocks, semicolon text fields), so this module reuses
//! `crate::cif`'s existing tokenizer (`tokenize_cif`, `strip_cif_comment`,
//! `strip_esd`/`parse_esd`, `resolve_element`) rather than re-implementing
//! it, and reuses [`crate::cif::UnitCell`] for `_cell.*` bookkeeping. The
//! category-tag scanning (dot convention, `atom_site`/`cell`/`symmetry`
//! category handling) is new, mmCIF-specific work.
//!
//! ## Scope
//!
//! - Reads/writes the `_atom_site` loop: `group_PDB`, `id`, `type_symbol`,
//!   `label_atom_id`/`auth_atom_id`, `label_alt_id`, `label_comp_id`/
//!   `auth_comp_id`, `label_asym_id`/`auth_asym_id`, `label_entity_id`,
//!   `label_seq_id`/`auth_seq_id`, `pdbx_PDB_ins_code`, `Cartn_x/y/z`,
//!   `occupancy`, `B_iso_or_equiv`, `pdbx_formal_charge`,
//!   `pdbx_PDB_model_num`.
//! - Reads/writes `_cell.length_a/b/c`, `_cell.angle_alpha/beta/gamma`
//!   (reusing [`crate::cif::UnitCell`]) and `_symmetry.space_group_name_H-M`,
//!   when present. These are *not* used to transform `_atom_site.Cartn_x/y/z`
//!   -- per the mmCIF dictionary, `_atom_site.Cartn_*` is always given in
//!   orthogonal Ångströms already (unlike small-molecule CIF's
//!   `_atom_site_fract_x`, which this module's sibling does convert).
//! - Only the first `data_` block is read (matches `crate::cif`'s existing
//!   scope note).
//! - **Judgment call:** `label_*` and `auth_*` variants of `atom_id`/
//!   `comp_id`/`asym_id` are collapsed to one canonical value per atom
//!   (`auth_*` preferred, falling back to `label_*` -- both are treated as
//!   "no value" if empty or a `.`/`?` placeholder, so the placeholder marker
//!   itself never leaks into the stored value) rather than tracked as two
//!   independently-round-tripping values. This matches what "chain id"/
//!   "residue name" mean to most consumers (the human-facing
//!   author-assigned identifiers) and is what the vast majority of files
//!   have `label_*` == `auth_*` anyway; [`write_mmcif`] writes the single
//!   stored value into both the `label_*` and `auth_*` columns, so read ->
//!   write -> read is still an exact fixed point. A file whose
//!   `label_*`/`auth_*` genuinely disagree loses that distinction on
//!   round-trip -- this is a deliberate simplification, not a silent bug.
//!   `seq_id` is the one exception: [`MmcifAtomRecord::label_seq_id`] is
//!   kept as its own field rather than folded into this collapse, because
//!   `label_seq_id` is legitimately `.` for essentially every non-polymer
//!   atom (waters, ligands, ions) even when `auth_seq_id` is a normal
//!   number, and fabricating an integer there on write would be a real
//!   semantic change, not just a lost distinction.
//! - **Known inherited limitation:** [`crate::cif`]'s `strip_cif_comment`
//!   strips `#`-to-end-of-line per physical line, which does not know about
//!   multi-line semicolon (`;...;`) text fields -- a `#` inside one would be
//!   (incorrectly) treated as a comment delimiter. Harmless for every
//!   category this module reads (none of `_atom_site`/`_cell`/`_symmetry`'s
//!   values are typically given as semicolon text blocks), but semicolon
//!   blocks are common elsewhere in real mmCIF (titles, sequences,
//!   citations) -- this module doesn't read those categories at all, so the
//!   gap is currently unreachable here, not fixed.
//! - No bond table: mmCIF's `_atom_site` category carries no explicit
//!   connectivity (unlike a MOL/SDF Ctab). This module never infers bonds,
//!   not even opt-in -- callers needing bonds should run their own
//!   perception step (e.g. distance-based) against the returned coordinates.
//! - Unknown/unmapped `_atom_site` loop columns are never silently dropped:
//!   their tag names are collected into [`MmcifResult::unhandled_columns`]
//!   so a caller can tell the file carried more than this reader modeled.
//!   Categories entirely outside this module's scope (`_entity.*` beyond a
//!   per-atom id, `_struct_conn`, `_pdbx_struct_assembly*`, ...) are simply
//!   not read at all -- this is a scope boundary stated up front, not a
//!   silent drop of something this module claimed to understand.

use chematic_core::{Atom, Element, Molecule, MoleculeBuilder};

use crate::cif::{
    UnitCell, parse_esd, resolve_element, strip_cif_comment, strip_esd, tokenize_cif,
};

// ---------------------------------------------------------------------------
// Limits
// ---------------------------------------------------------------------------

/// Security/robustness limits enforced before any chemistry is interpreted.
/// Mirrors the shape of [`crate::mrv::MrvParseLimits`].
#[derive(Debug, Clone, Copy)]
pub struct MmcifParseLimits {
    /// Maximum input size in bytes.
    pub max_input_bytes: usize,
    /// Maximum number of `_atom_site` rows (across all models).
    pub max_atoms: usize,
    /// Maximum length of any single line, in bytes.
    pub max_line_len: usize,
}

impl Default for MmcifParseLimits {
    fn default() -> Self {
        Self {
            // Macromolecular mmCIF files (large assemblies, cryo-EM models)
            // routinely run tens of MB; small-molecule CIF's 32 MiB default
            // (see MrvParseLimits) would be too tight here.
            max_input_bytes: 128 << 20, // 128 MiB
            max_atoms: 2_000_000,
            max_line_len: 8192,
        }
    }
}

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Errors that can occur when parsing mmCIF files.
#[derive(Debug, Clone, PartialEq)]
pub enum MmcifError {
    /// Input exceeded [`MmcifParseLimits::max_input_bytes`].
    InputTooLarge { limit: usize },
    /// A single line exceeded [`MmcifParseLimits::max_line_len`].
    LineTooLong { line: usize, limit: usize },
    /// `_atom_site` row count exceeded [`MmcifParseLimits::max_atoms`].
    TooManyAtoms { limit: usize },
    /// No `loop_` block with `_atom_site.*` tags was found.
    NoAtomSiteLoop,
    /// The `_atom_site` loop lacked a `type_symbol` column (element is
    /// mandatory in the mmCIF dictionary).
    MissingTypeSymbolColumn,
    /// The `_atom_site` loop lacked one or more of `Cartn_x`/`Cartn_y`/`Cartn_z`.
    MissingCoordinateColumns,
    /// An element symbol could not be resolved.
    UnknownElement(String),
    /// A coordinate/occupancy/B-factor value could not be parsed as a
    /// finite float.
    InvalidNumber { column: &'static str, raw: String },
    /// An integer-valued column (`id`, `auth_seq_id`, `pdbx_PDB_model_num`,
    /// `pdbx_formal_charge`) could not be parsed as an integer.
    InvalidInteger { column: &'static str, raw: String },
}

impl core::fmt::Display for MmcifError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InputTooLarge { limit } => write!(f, "mmCIF input exceeds {limit}-byte limit"),
            Self::LineTooLong { line, limit } => {
                write!(f, "mmCIF line {line} exceeds {limit}-byte limit")
            }
            Self::TooManyAtoms { limit } => {
                write!(f, "mmCIF atom_site row count exceeds {limit}-atom limit")
            }
            Self::NoAtomSiteLoop => write!(f, "no _atom_site.* loop_ block found in mmCIF"),
            Self::MissingTypeSymbolColumn => {
                write!(
                    f,
                    "_atom_site loop is missing the required type_symbol column"
                )
            }
            Self::MissingCoordinateColumns => {
                write!(
                    f,
                    "_atom_site loop is missing Cartn_x/Cartn_y/Cartn_z columns"
                )
            }
            Self::UnknownElement(s) => write!(f, "unknown element '{s}' in mmCIF atom_site"),
            Self::InvalidNumber { column, raw } => {
                write!(f, "invalid {column} value '{raw}' in mmCIF atom_site")
            }
            Self::InvalidInteger { column, raw } => {
                write!(
                    f,
                    "invalid integer {column} value '{raw}' in mmCIF atom_site"
                )
            }
        }
    }
}

impl std::error::Error for MmcifError {}

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// One `_atom_site` row, in file order.
#[derive(Debug, Clone, PartialEq)]
pub struct MmcifAtomRecord {
    /// `_atom_site.group_PDB` -- `"ATOM"` or `"HETATM"` (or, rarely, some
    /// other file-specific value; kept verbatim rather than restricted to
    /// an enum so an unusual value round-trips instead of erroring).
    pub group_pdb: String,
    /// `_atom_site.id` -- the atom serial number.
    pub serial: i64,
    /// Resolved from `_atom_site.type_symbol`.
    pub element: Element,
    /// `_atom_site.auth_atom_id`, falling back to `.label_atom_id`.
    pub atom_name: String,
    /// `_atom_site.label_alt_id` -- `None` when absent or `"."`/`"?"`.
    pub alt_loc: Option<char>,
    /// `_atom_site.auth_comp_id`, falling back to `.label_comp_id`.
    pub res_name: String,
    /// `_atom_site.auth_asym_id`, falling back to `.label_asym_id`.
    pub chain_id: String,
    /// `_atom_site.auth_seq_id`, falling back to `.label_seq_id`, defaulting
    /// to `0` if neither is present/non-placeholder (rare -- `label_seq_id`
    /// is dictionary-mandated, so this only happens on a file this module
    /// doesn't fully understand).
    pub res_seq: i64,
    /// `_atom_site.label_seq_id` specifically, kept separate from
    /// [`Self::res_seq`] rather than collapsed into it: `label_seq_id` is
    /// `.` (`None` here) for essentially every non-polymer atom (waters,
    /// ligands, ions) in real mmCIF even when `auth_seq_id` (-> `res_seq`)
    /// is a normal number -- collapsing the two would silently turn every
    /// such atom's "not applicable" into a fabricated `label_seq_id`
    /// integer on write. `None` when the column is absent, `.`, or `?`.
    pub label_seq_id: Option<i64>,
    /// `_atom_site.pdbx_PDB_ins_code` -- `None` when absent or `"."`/`"?"`.
    pub icode: Option<char>,
    /// `_atom_site.Cartn_x` (Å).
    pub x: f64,
    /// `_atom_site.Cartn_y` (Å).
    pub y: f64,
    /// `_atom_site.Cartn_z` (Å).
    pub z: f64,
    /// `_atom_site.occupancy`, defaulting to `1.0` when the column is
    /// absent (a `.`/`?` placeholder is also treated as `1.0`).
    pub occupancy: f64,
    /// `_atom_site.B_iso_or_equiv`, defaulting to `0.0` when the column is
    /// absent (a `.`/`?` placeholder is also treated as `0.0`).
    pub b_iso: f64,
    /// `_atom_site.pdbx_formal_charge`, if present and not `.`/`?`.
    pub formal_charge: Option<i32>,
    /// `_atom_site.label_entity_id`, if present.
    pub entity_id: Option<String>,
    /// `_atom_site.pdbx_PDB_model_num`, defaulting to `1` when the column
    /// is absent (a file with no model column is, by definition, single-
    /// model).
    pub model_num: i32,
}

/// Result of parsing an mmCIF file.
#[derive(Debug, Clone, PartialEq)]
pub struct MmcifResult {
    /// `_atom_site` rows, in file order (spanning every model, if the file
    /// has more than one -- filter on [`MmcifAtomRecord::model_num`]).
    pub atoms: Vec<MmcifAtomRecord>,
    /// Unit cell parameters from `_cell.length_*`/`_cell.angle_*`, if all
    /// six were present.
    pub cell: Option<UnitCell>,
    /// `_symmetry.space_group_name_H-M`, if present (CIF quoting already
    /// stripped by the tokenizer).
    pub space_group: Option<String>,
    /// `_atom_site` loop column tags this reader saw but does not model
    /// (e.g. `_atom_site.pdbx_PDB_strand_id`-adjacent extensions, ADP
    /// columns). Never silently dropped -- surfaced here instead of a
    /// per-atom raw pass-through field (see module docs).
    pub unhandled_columns: Vec<String>,
}

impl MmcifResult {
    /// Build a plain [`Molecule`] + per-atom Cartesian coordinates from
    /// [`Self::atoms`], in the same file order. No bonds are added (see
    /// module docs). Includes every model's atoms; filter `self.atoms` by
    /// `model_num` first if only one model is wanted.
    pub fn to_molecule(&self) -> (Molecule, Vec<(f64, f64, f64)>) {
        let mut builder = MoleculeBuilder::new();
        let mut coords = Vec::with_capacity(self.atoms.len());
        for a in &self.atoms {
            builder.add_atom(Atom::new(a.element));
            coords.push((a.x, a.y, a.z));
        }
        (builder.build(), coords)
    }
}

// ---------------------------------------------------------------------------
// Tag helpers (DDL2 `_category.field` convention)
// ---------------------------------------------------------------------------

/// Split a `_category.field` tag into `(category, field)`, both borrowed.
/// Returns `None` for tags without a `.` (e.g. DDL1-style small-molecule
/// tags, or `loop_`/`data_` keywords) -- those simply aren't mmCIF category
/// tags as far as this module is concerned.
fn split_tag(tag: &str) -> Option<(&str, &str)> {
    let rest = tag.strip_prefix('_')?;
    let dot = rest.find('.')?;
    Some((&rest[..dot], &rest[dot + 1..]))
}

/// Find the first `loop_` block whose header tags belong to `category`
/// (case-insensitive, e.g. `"atom_site"`). Returns `(lowercased field
/// names with the category prefix stripped, token index where data rows
/// begin)`. mmCIF loops are always single-category by dictionary
/// convention, so only the first header tag needs to be checked.
fn find_loop_by_category(tokens: &[String], category: &str) -> Option<(Vec<String>, usize)> {
    let mut i = 0;
    while i < tokens.len() {
        if tokens[i] == "loop_" {
            let mut j = i + 1;
            let mut headers: Vec<String> = Vec::new();
            while j < tokens.len() && tokens[j].starts_with('_') {
                headers.push(tokens[j].clone());
                j += 1;
            }
            if let Some((cat, _)) = headers.first().and_then(|h| split_tag(h))
                && cat.eq_ignore_ascii_case(category)
            {
                let fields = headers
                    .iter()
                    .map(|h| {
                        split_tag(h)
                            .map(|(_, f)| f.to_ascii_lowercase())
                            .unwrap_or_default()
                    })
                    .collect();
                return Some((fields, j));
            }
            i = j;
        } else {
            i += 1;
        }
    }
    None
}

/// Scan `tokens` for bare (non-loop) `_cell.*`/`_symmetry.*` singleton
/// tag/value pairs. Mirrors [`crate::cif::scan_cell`]'s shape for the
/// dot-tag convention.
fn scan_cell_and_symmetry(tokens: &[String]) -> (Option<UnitCell>, Option<String>) {
    let mut cell = UnitCell::default();
    let (mut a, mut b, mut c, mut alpha, mut beta, mut gamma) =
        (false, false, false, false, false, false);
    let mut space_group: Option<String> = None;

    let mut i = 0;
    while i + 1 < tokens.len() {
        if let Some((cat, field)) = split_tag(&tokens[i]) {
            let field_lc = field.to_ascii_lowercase();
            if cat.eq_ignore_ascii_case("cell") {
                match field_lc.as_str() {
                    "length_a" => {
                        if let Some(v) = parse_esd(&tokens[i + 1]) {
                            cell.a = v;
                            a = true;
                        }
                    }
                    "length_b" => {
                        if let Some(v) = parse_esd(&tokens[i + 1]) {
                            cell.b = v;
                            b = true;
                        }
                    }
                    "length_c" => {
                        if let Some(v) = parse_esd(&tokens[i + 1]) {
                            cell.c = v;
                            c = true;
                        }
                    }
                    "angle_alpha" => {
                        if let Some(v) = parse_esd(&tokens[i + 1]) {
                            cell.alpha = v;
                            alpha = true;
                        }
                    }
                    "angle_beta" => {
                        if let Some(v) = parse_esd(&tokens[i + 1]) {
                            cell.beta = v;
                            beta = true;
                        }
                    }
                    "angle_gamma" => {
                        if let Some(v) = parse_esd(&tokens[i + 1]) {
                            cell.gamma = v;
                            gamma = true;
                        }
                    }
                    _ => {}
                }
            } else if cat.eq_ignore_ascii_case("symmetry") && field_lc == "space_group_name_h-m" {
                space_group = Some(tokens[i + 1].clone());
            }
        }
        i += 1;
    }

    let cell_opt = if a && b && c && alpha && beta && gamma {
        Some(cell)
    } else {
        None
    };
    (cell_opt, space_group)
}

/// `.`/`?` are CIF's "inapplicable"/"unknown" placeholders; both are
/// treated as "value not given" here.
fn is_cif_placeholder(s: &str) -> bool {
    s == "." || s == "?"
}

fn opt_str(s: &str) -> Option<String> {
    if is_cif_placeholder(s) {
        None
    } else {
        Some(s.to_string())
    }
}

fn opt_char(s: &str) -> Option<char> {
    opt_str(s).and_then(|v| v.chars().next())
}

/// Pick the `auth_*` value, falling back to `label_*`, treating an absent
/// column, empty string, and CIF's `.`/`?` placeholders all as "no value"
/// on *both* sides -- not just on `auth`. Returns `""` if neither source
/// has a real value.
///
/// The naive version of this fallback (checking only `auth` for
/// emptiness/placeholder and using `label` unconditionally) leaks CIF's
/// `.` placeholder into the stored data whenever `label_*` itself happens
/// to be a placeholder too (e.g. a file with `label_asym_id` absent and
/// `auth_asym_id` present would be fine, but the reverse -- `auth_asym_id`
/// literally `.` and `label_asym_id` also `.` -- would previously store
/// the string `"."` as a chain id). Returning `""` in that case instead
/// keeps [`write_mmcif`]'s round trip a fixed point: `""` writes back out
/// as `.` (`quote_cif_value("")`) and re-parses as `""` again.
fn pick(auth: Option<&str>, label: Option<&str>) -> String {
    for candidate in [auth, label].into_iter().flatten() {
        if !candidate.is_empty() && !is_cif_placeholder(candidate) {
            return candidate.to_string();
        }
    }
    String::new()
}

/// `_atom_site.type_symbol` is dictionary-mandated to carry the IUPAC
/// title-case symbol (`"Zn"`), but real deposited/tool-written files
/// sometimes use all-caps (`"ZN"`) instead -- `resolve_element` (shared
/// with `crate::cif`) is strict, so retry with a title-cased fallback
/// rather than rejecting an otherwise-unambiguous, very common case. This
/// wrapper is local to mmCIF reading; the shared `resolve_element` itself
/// is left as-is for `crate::cif`'s small-molecule callers.
fn resolve_atom_site_element(raw: &str) -> Result<Element, MmcifError> {
    resolve_element(raw)
        .or_else(|_| {
            let mut chars = raw.chars();
            let title_cased = match chars.next() {
                Some(first) => {
                    first.to_ascii_uppercase().to_string() + &chars.as_str().to_ascii_lowercase()
                }
                None => String::new(),
            };
            resolve_element(&title_cased)
        })
        .map_err(|e| MmcifError::UnknownElement(e.to_string()))
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Parse an mmCIF file with default limits ([`MmcifParseLimits::default`]).
pub fn parse_mmcif(input: &str) -> Result<MmcifResult, MmcifError> {
    parse_mmcif_with_limits(input, &MmcifParseLimits::default())
}

/// Parse an mmCIF file, enforcing `limits`. See module docs for scope.
pub fn parse_mmcif_with_limits(
    input: &str,
    limits: &MmcifParseLimits,
) -> Result<MmcifResult, MmcifError> {
    if input.len() > limits.max_input_bytes {
        return Err(MmcifError::InputTooLarge {
            limit: limits.max_input_bytes,
        });
    }
    for (lineno, line) in input.lines().enumerate() {
        if line.len() > limits.max_line_len {
            return Err(MmcifError::LineTooLong {
                line: lineno + 1,
                limit: limits.max_line_len,
            });
        }
    }

    let clean: String = input
        .lines()
        .map(strip_cif_comment)
        .collect::<Vec<_>>()
        .join("\n");
    let tokens = tokenize_cif(&clean);

    let (cell, space_group) = scan_cell_and_symmetry(&tokens);

    let (headers, data_start) =
        find_loop_by_category(&tokens, "atom_site").ok_or(MmcifError::NoAtomSiteLoop)?;
    let ncols = headers.len();
    let col = |name: &str| -> Option<usize> { headers.iter().position(|h| h == name) };

    let col_group_pdb = col("group_pdb");
    let col_id = col("id");
    let col_type_symbol = col("type_symbol").ok_or(MmcifError::MissingTypeSymbolColumn)?;
    let col_label_atom_id = col("label_atom_id");
    let col_auth_atom_id = col("auth_atom_id");
    let col_alt_id = col("label_alt_id");
    let col_label_comp_id = col("label_comp_id");
    let col_auth_comp_id = col("auth_comp_id");
    let col_label_asym_id = col("label_asym_id");
    let col_auth_asym_id = col("auth_asym_id");
    let col_entity_id = col("label_entity_id");
    let col_label_seq_id = col("label_seq_id");
    let col_auth_seq_id = col("auth_seq_id");
    let col_ins_code = col("pdbx_pdb_ins_code");
    let col_x = col("cartn_x");
    let col_y = col("cartn_y");
    let col_z = col("cartn_z");
    let col_occ = col("occupancy");
    let col_biso = col("b_iso_or_equiv");
    let col_charge = col("pdbx_formal_charge");
    let col_model = col("pdbx_pdb_model_num");

    let (col_x, col_y, col_z) = match (col_x, col_y, col_z) {
        (Some(x), Some(y), Some(z)) => (x, y, z),
        _ => return Err(MmcifError::MissingCoordinateColumns),
    };

    // Known columns (by loop-header position) -- everything else is
    // reported via `unhandled_columns` rather than silently ignored.
    let known: Vec<usize> = [
        col_group_pdb,
        col_id,
        Some(col_type_symbol),
        col_label_atom_id,
        col_auth_atom_id,
        col_alt_id,
        col_label_comp_id,
        col_auth_comp_id,
        col_label_asym_id,
        col_auth_asym_id,
        col_entity_id,
        col_label_seq_id,
        col_auth_seq_id,
        col_ins_code,
        Some(col_x),
        Some(col_y),
        Some(col_z),
        col_occ,
        col_biso,
        col_charge,
        col_model,
    ]
    .into_iter()
    .flatten()
    .collect();
    let unhandled_columns: Vec<String> = headers
        .iter()
        .enumerate()
        .filter(|(i, _)| !known.contains(i))
        .map(|(_, h)| format!("_atom_site.{h}"))
        .collect();

    let data_tokens = &tokens[data_start..];
    let mut atoms: Vec<MmcifAtomRecord> = Vec::new();
    let mut row = 0usize;
    while (row + 1) * ncols <= data_tokens.len() {
        let base = row * ncols;
        let tok = &data_tokens[base..base + ncols];

        // Stop at the next loop_/data_ block, or a bare singleton tag
        // (`_category.field value`) -- stricter than crate::cif's small-
        // molecule loop scanner, since real mmCIF files commonly place a
        // bare tag directly after the last loop with no separating
        // `loop_`/`data_` token.
        if tok[0] == "loop_" || tok[0].starts_with("data_") || tok[0].starts_with('_') {
            break;
        }
        if row >= limits.max_atoms {
            return Err(MmcifError::TooManyAtoms {
                limit: limits.max_atoms,
            });
        }

        let get = |c: usize| -> &str { tok[c].as_str() };

        let group_pdb = col_group_pdb
            .map(|c| get(c).to_string())
            .unwrap_or_else(|| "ATOM".to_string());
        let serial = match col_id {
            Some(c) => parse_int(get(c), "id")?,
            None => row as i64 + 1,
        };
        let element = resolve_atom_site_element(get(col_type_symbol))?;

        let atom_name = pick(col_auth_atom_id.map(get), col_label_atom_id.map(get));
        let alt_loc = col_alt_id.map(get).and_then(opt_char);
        let res_name = pick(col_auth_comp_id.map(get), col_label_comp_id.map(get));
        let chain_id = pick(col_auth_asym_id.map(get), col_label_asym_id.map(get));
        let entity_id = col_entity_id.map(get).and_then(opt_str);

        let label_seq_id = col_label_seq_id
            .map(get)
            .filter(|v| !is_cif_placeholder(v))
            .map(|v| parse_int(v, "label_seq_id"))
            .transpose()?;
        let res_seq = match col_auth_seq_id.map(get) {
            Some(v) if !is_cif_placeholder(v) => parse_int(v, "auth_seq_id")?,
            _ => label_seq_id.unwrap_or(0),
        };

        let icode = col_ins_code.map(get).and_then(opt_char);

        let x = parse_finite(get(col_x), "Cartn_x")?;
        let y = parse_finite(get(col_y), "Cartn_y")?;
        let z = parse_finite(get(col_z), "Cartn_z")?;

        let occupancy = match col_occ.map(get) {
            Some(v) if !is_cif_placeholder(v) => parse_finite(v, "occupancy")?,
            _ => 1.0,
        };
        let b_iso = match col_biso.map(get) {
            Some(v) if !is_cif_placeholder(v) => parse_finite(v, "B_iso_or_equiv")?,
            _ => 0.0,
        };
        let formal_charge = match col_charge.map(get) {
            Some(v) if !is_cif_placeholder(v) => Some(parse_int(v, "pdbx_formal_charge")? as i32),
            _ => None,
        };
        let model_num = match col_model.map(get) {
            Some(v) if !is_cif_placeholder(v) => parse_int(v, "pdbx_PDB_model_num")? as i32,
            _ => 1,
        };

        atoms.push(MmcifAtomRecord {
            group_pdb,
            serial,
            element,
            atom_name,
            alt_loc,
            res_name,
            chain_id,
            res_seq,
            label_seq_id,
            icode,
            x,
            y,
            z,
            occupancy,
            b_iso,
            formal_charge,
            entity_id,
            model_num,
        });
        row += 1;
    }

    if atoms.is_empty() {
        return Err(MmcifError::NoAtomSiteLoop);
    }

    Ok(MmcifResult {
        atoms,
        cell,
        space_group,
        unhandled_columns,
    })
}

fn parse_finite(s: &str, column: &'static str) -> Result<f64, MmcifError> {
    let v = strip_esd(s)
        .parse::<f64>()
        .map_err(|_| MmcifError::InvalidNumber {
            column,
            raw: s.to_string(),
        })?;
    if !v.is_finite() {
        return Err(MmcifError::InvalidNumber {
            column,
            raw: s.to_string(),
        });
    }
    Ok(v)
}

/// Parse an mmCIF integer field. The `pdbx_formal_charge` column is
/// dictionary-mandated to be a plain signed integer (`"2"`, `"-1"`), but
/// files converted from legacy PDB carry the old postfix-sign notation
/// from PDB columns 79-80 instead (`"2+"`, `"1-"`) -- accepted here too
/// (postfix sign only; a value can't have both a prefix and postfix sign).
/// This widens what can be read; it never affects what [`write_mmcif`]
/// produces (always plain signed-integer form).
fn parse_int(s: &str, column: &'static str) -> Result<i64, MmcifError> {
    let s = strip_esd(s).trim();
    if let Some(digits) = s.strip_suffix('+') {
        return digits
            .parse::<i64>()
            .map_err(|_| MmcifError::InvalidInteger {
                column,
                raw: s.to_string(),
            });
    }
    if let Some(digits) = s.strip_suffix('-') {
        return digits
            .parse::<i64>()
            .map(|v| -v)
            .map_err(|_| MmcifError::InvalidInteger {
                column,
                raw: s.to_string(),
            });
    }
    s.parse::<i64>().map_err(|_| MmcifError::InvalidInteger {
        column,
        raw: s.to_string(),
    })
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Quote a CIF value if required (contains whitespace, is empty, or starts
/// with a character that is otherwise syntactically significant). Prefers
/// an unquoted token; falls back to single quotes, then double quotes if
/// the value itself contains a single quote.
fn quote_cif_value(s: &str) -> String {
    if s.is_empty() {
        return ".".to_string();
    }
    let needs_quoting = s.chars().any(char::is_whitespace)
        || matches!(s, "." | "?")
        || s.starts_with(['_', '#', '$', '\'', '"', ';', '[', ']']);
    if !needs_quoting {
        return s.to_string();
    }
    if !s.contains('\'') {
        format!("'{s}'")
    } else if !s.contains('"') {
        format!("\"{s}\"")
    } else {
        // Contains both quote characters -- fall back to a semicolon text
        // field is the fully-correct CIF answer, but no field this module
        // writes can realistically contain both quote characters (atom
        // names, residue names, chain ids); documented rather than
        // implemented.
        format!("'{s}'")
    }
}

/// Write an mmCIF file from parsed atom records plus optional cell/symmetry.
///
/// Writes both `label_*` and `auth_*` columns with the same value (see the
/// "Judgment call" note in the module docs). `data_block_name` becomes the
/// `data_<name>` header; non-alphanumeric characters are not sanitized
/// (callers should pass a valid mmCIF data-block name).
pub fn write_mmcif(
    atoms: &[MmcifAtomRecord],
    cell: Option<&UnitCell>,
    space_group: Option<&str>,
    data_block_name: &str,
) -> String {
    let mut out = format!("data_{data_block_name}\n#\n");

    if let Some(c) = cell {
        out.push_str(&format!(
            "_cell.length_a      {:.4}\n\
             _cell.length_b      {:.4}\n\
             _cell.length_c      {:.4}\n\
             _cell.angle_alpha   {:.3}\n\
             _cell.angle_beta    {:.3}\n\
             _cell.angle_gamma   {:.3}\n#\n",
            c.a, c.b, c.c, c.alpha, c.beta, c.gamma
        ));
    }
    if let Some(sg) = space_group {
        out.push_str(&format!(
            "_symmetry.space_group_name_H-M   {}\n#\n",
            quote_cif_value(sg)
        ));
    }

    out.push_str(
        "loop_\n\
         _atom_site.group_PDB\n\
         _atom_site.id\n\
         _atom_site.type_symbol\n\
         _atom_site.label_atom_id\n\
         _atom_site.label_alt_id\n\
         _atom_site.label_comp_id\n\
         _atom_site.label_asym_id\n\
         _atom_site.label_entity_id\n\
         _atom_site.label_seq_id\n\
         _atom_site.pdbx_PDB_ins_code\n\
         _atom_site.Cartn_x\n\
         _atom_site.Cartn_y\n\
         _atom_site.Cartn_z\n\
         _atom_site.occupancy\n\
         _atom_site.B_iso_or_equiv\n\
         _atom_site.pdbx_formal_charge\n\
         _atom_site.auth_seq_id\n\
         _atom_site.auth_comp_id\n\
         _atom_site.auth_asym_id\n\
         _atom_site.auth_atom_id\n\
         _atom_site.pdbx_PDB_model_num\n",
    );

    for a in atoms {
        let alt = a
            .alt_loc
            .map(|c| c.to_string())
            .unwrap_or_else(|| ".".to_string());
        let icode = a
            .icode
            .map(|c| c.to_string())
            .unwrap_or_else(|| ".".to_string());
        let entity = a
            .entity_id
            .as_deref()
            .map(quote_cif_value)
            .unwrap_or_else(|| ".".to_string());
        let charge = a
            .formal_charge
            .map(|c| c.to_string())
            .unwrap_or_else(|| ".".to_string());
        // label_seq_id is written separately from auth_seq_id (res_seq) --
        // see MmcifAtomRecord::label_seq_id's doc comment: collapsing the
        // two would fabricate a label_seq_id for every non-polymer atom.
        let label_seq = a
            .label_seq_id
            .map(|v| v.to_string())
            .unwrap_or_else(|| ".".to_string());

        out.push_str(&format!(
            "{group} {serial} {sym} {name} {alt} {comp} {asym} {entity} {label_seq} {icode} \
             {x:.3} {y:.3} {z:.3} {occ:.2} {biso:.2} {charge} \
             {auth_seq} {comp} {asym} {name} {model}\n",
            group = quote_cif_value(&a.group_pdb),
            serial = a.serial,
            sym = a.element.symbol(),
            name = quote_cif_value(&a.atom_name),
            alt = alt,
            comp = quote_cif_value(&a.res_name),
            asym = quote_cif_value(&a.chain_id),
            entity = entity,
            label_seq = label_seq,
            icode = icode,
            x = a.x,
            y = a.y,
            z = a.z,
            occ = a.occupancy,
            biso = a.b_iso,
            charge = charge,
            auth_seq = a.res_seq,
            model = a.model_num,
        ));
    }

    out
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Hand-authored per the wwPDB mmCIF dictionary's documented
    /// `_atom_site` column layout (not copied from any real deposited PDB
    /// entry or tool output) -- a minimal two-residue, two-chain, one
    /// altloc-bearing fragment.
    const FIXTURE: &str = "data_TEST\n\
        #\n\
        _cell.length_a      61.200\n\
        _cell.length_b      61.200\n\
        _cell.length_c      82.500\n\
        _cell.angle_alpha   90.000\n\
        _cell.angle_beta    90.000\n\
        _cell.angle_gamma   120.000\n\
        #\n\
        _symmetry.space_group_name_H-M   'P 32 2 1'\n\
        #\n\
        loop_\n\
        _atom_site.group_PDB\n\
        _atom_site.id\n\
        _atom_site.type_symbol\n\
        _atom_site.label_atom_id\n\
        _atom_site.label_alt_id\n\
        _atom_site.label_comp_id\n\
        _atom_site.label_asym_id\n\
        _atom_site.label_entity_id\n\
        _atom_site.label_seq_id\n\
        _atom_site.pdbx_PDB_ins_code\n\
        _atom_site.Cartn_x\n\
        _atom_site.Cartn_y\n\
        _atom_site.Cartn_z\n\
        _atom_site.occupancy\n\
        _atom_site.B_iso_or_equiv\n\
        _atom_site.pdbx_formal_charge\n\
        _atom_site.auth_seq_id\n\
        _atom_site.auth_comp_id\n\
        _atom_site.auth_asym_id\n\
        _atom_site.auth_atom_id\n\
        _atom_site.pdbx_PDB_model_num\n\
        ATOM   1  N  N   . SER A 1 1  ? 10.123 20.456 5.789  1.00 25.30 ? 1 SER A N  1\n\
        ATOM   2  C  CA  A SER A 1 1  ? 11.234 21.567 6.890  0.60 24.10 ? 1 SER A CA 1\n\
        ATOM   3  C  CA  B SER A 1 1  ? 11.334 21.667 6.990  0.40 24.50 ? 1 SER A CA 1\n\
        HETATM 4  ZN ZN  . ZN  B 2 .  ? 30.000 30.000 30.000 1.00 15.00 2 101 ZN  B ZN 1\n\
        ATOM   5  N  N   . SER A 1 1  ? 10.123 20.456 15.789 1.00 25.30 ? 1 SER A N  2\n\
        ATOM   6  C  CA  . SER A 1 1  ? 11.234 21.567 16.890 1.00 24.10 ? 1 SER A CA 2\n";

    #[test]
    fn parses_atom_count_and_models() {
        let r = parse_mmcif(FIXTURE).unwrap();
        assert_eq!(r.atoms.len(), 6);
        let models: Vec<i32> = r.atoms.iter().map(|a| a.model_num).collect();
        assert_eq!(models, vec![1, 1, 1, 1, 2, 2]);
    }

    #[test]
    fn parses_cell_and_symmetry() {
        let r = parse_mmcif(FIXTURE).unwrap();
        let cell = r.cell.unwrap();
        assert!((cell.a - 61.2).abs() < 1e-6);
        assert!((cell.gamma - 120.0).abs() < 1e-6);
        assert_eq!(r.space_group.as_deref(), Some("P 32 2 1"));
    }

    #[test]
    fn parses_chain_residue_altloc_and_charge() {
        let r = parse_mmcif(FIXTURE).unwrap();
        assert_eq!(r.atoms[0].chain_id, "A");
        assert_eq!(r.atoms[0].res_name, "SER");
        assert_eq!(r.atoms[0].res_seq, 1);
        assert_eq!(r.atoms[0].alt_loc, None);
        assert_eq!(r.atoms[1].alt_loc, Some('A'));
        assert_eq!(r.atoms[2].alt_loc, Some('B'));
        assert!((r.atoms[1].occupancy - 0.60).abs() < 1e-9);
        assert!((r.atoms[2].occupancy - 0.40).abs() < 1e-9);
        assert_eq!(r.atoms[3].chain_id, "B");
        assert_eq!(r.atoms[3].group_pdb, "HETATM");
        assert_eq!(r.atoms[3].formal_charge, Some(2));
        assert_eq!(r.atoms[3].element, Element::ZN);
    }

    #[test]
    fn insertion_code_placeholder_is_none() {
        let r = parse_mmcif(FIXTURE).unwrap();
        assert_eq!(r.atoms[0].icode, None);
    }

    #[test]
    fn unhandled_columns_is_empty_for_fully_modeled_fixture() {
        let r = parse_mmcif(FIXTURE).unwrap();
        assert!(r.unhandled_columns.is_empty());
    }

    #[test]
    fn unhandled_columns_reports_unmapped_tag() {
        let cif = "data_x\n\
            loop_\n\
            _atom_site.group_PDB\n\
            _atom_site.id\n\
            _atom_site.type_symbol\n\
            _atom_site.Cartn_x\n\
            _atom_site.Cartn_y\n\
            _atom_site.Cartn_z\n\
            _atom_site.pdbx_PDB_strand_id\n\
            ATOM 1 C 0.0 0.0 0.0 X\n";
        let r = parse_mmcif(cif).unwrap();
        assert_eq!(r.unhandled_columns, vec!["_atom_site.pdbx_pdb_strand_id"]);
    }

    #[test]
    fn to_molecule_builds_atoms_with_no_bonds() {
        let r = parse_mmcif(FIXTURE).unwrap();
        let (mol, coords) = r.to_molecule();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(coords.len(), 6);
        assert_eq!(mol.atom(chematic_core::AtomIdx(0)).element, Element::N);
        assert!((coords[0].0 - 10.123).abs() < 1e-6);
    }

    #[test]
    fn round_trip_preserves_metadata() {
        let r = parse_mmcif(FIXTURE).unwrap();
        let written = write_mmcif(&r.atoms, r.cell.as_ref(), r.space_group.as_deref(), "TEST");
        let r2 = parse_mmcif(&written).unwrap();
        assert_eq!(r2.atoms.len(), r.atoms.len());
        for (a, b) in r.atoms.iter().zip(r2.atoms.iter()) {
            assert_eq!(a.group_pdb, b.group_pdb);
            assert_eq!(a.serial, b.serial);
            assert_eq!(a.element, b.element);
            assert_eq!(a.atom_name, b.atom_name);
            assert_eq!(a.alt_loc, b.alt_loc);
            assert_eq!(a.res_name, b.res_name);
            assert_eq!(a.chain_id, b.chain_id);
            assert_eq!(a.res_seq, b.res_seq);
            assert_eq!(a.label_seq_id, b.label_seq_id);
            assert_eq!(a.icode, b.icode);
            assert!((a.x - b.x).abs() < 1e-3);
            assert!((a.y - b.y).abs() < 1e-3);
            assert!((a.z - b.z).abs() < 1e-3);
            assert!((a.occupancy - b.occupancy).abs() < 1e-3);
            assert!((a.b_iso - b.b_iso).abs() < 1e-3);
            assert_eq!(a.formal_charge, b.formal_charge);
            assert_eq!(a.model_num, b.model_num);
        }
        let cell1 = r.cell.unwrap();
        let cell2 = r2.cell.unwrap();
        assert!((cell1.a - cell2.a).abs() < 1e-3);
        assert!((cell1.gamma - cell2.gamma).abs() < 1e-3);
        assert_eq!(r.space_group, r2.space_group);
    }

    #[test]
    fn entity_id_round_trips() {
        let cif = "data_x\n\
            loop_\n\
            _atom_site.group_PDB\n\
            _atom_site.id\n\
            _atom_site.type_symbol\n\
            _atom_site.label_entity_id\n\
            _atom_site.Cartn_x\n\
            _atom_site.Cartn_y\n\
            _atom_site.Cartn_z\n\
            ATOM 1 C 1 0.0 0.0 0.0\n\
            HETATM 2 ZN 2 5.0 5.0 5.0\n";
        let r = parse_mmcif(cif).unwrap();
        assert_eq!(r.atoms[0].entity_id.as_deref(), Some("1"));
        assert_eq!(r.atoms[1].entity_id.as_deref(), Some("2"));
        let written = write_mmcif(&r.atoms, None, None, "x");
        let r2 = parse_mmcif(&written).unwrap();
        assert_eq!(r2.atoms[0].entity_id.as_deref(), Some("1"));
        assert_eq!(r2.atoms[1].entity_id.as_deref(), Some("2"));
    }

    #[test]
    fn all_caps_type_symbol_still_resolves() {
        // Real-world files vary between IUPAC title case ("Zn") and
        // all-caps ("ZN") for two-letter type_symbol values.
        let cif = "data_x\n\
            loop_\n_atom_site.type_symbol\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
            ZN 0.0 0.0 0.0\n";
        let r = parse_mmcif(cif).unwrap();
        assert_eq!(r.atoms[0].element, Element::ZN);
    }

    #[test]
    fn missing_atom_site_loop_is_typed_error_not_panic() {
        let err = parse_mmcif("data_empty\n#\nno atom site here\n").unwrap_err();
        assert_eq!(err, MmcifError::NoAtomSiteLoop);
    }

    #[test]
    fn missing_coordinate_columns_is_typed_error() {
        let cif = "data_x\nloop_\n_atom_site.group_PDB\n_atom_site.type_symbol\nATOM C\n";
        let err = parse_mmcif(cif).unwrap_err();
        assert_eq!(err, MmcifError::MissingCoordinateColumns);
    }

    #[test]
    fn unknown_element_is_typed_error_not_panic() {
        let cif = "data_x\n\
            loop_\n_atom_site.type_symbol\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
            Qq 0.0 0.0 0.0\n";
        let err = parse_mmcif(cif).unwrap_err();
        assert!(matches!(err, MmcifError::UnknownElement(_)));
    }

    #[test]
    fn nan_coordinate_is_rejected() {
        let cif = "data_x\n\
            loop_\n_atom_site.type_symbol\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
            C NaN 0.0 0.0\n";
        let err = parse_mmcif(cif).unwrap_err();
        assert!(matches!(
            err,
            MmcifError::InvalidNumber {
                column: "Cartn_x",
                ..
            }
        ));
    }

    #[test]
    fn infinity_coordinate_is_rejected() {
        let cif = "data_x\n\
            loop_\n_atom_site.type_symbol\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
            C inf 0.0 0.0\n";
        let err = parse_mmcif(cif).unwrap_err();
        assert!(matches!(
            err,
            MmcifError::InvalidNumber {
                column: "Cartn_x",
                ..
            }
        ));
    }

    #[test]
    fn malformed_input_never_panics() {
        // A grab-bag of adversarial inputs that must all fail closed
        // (Err), never panic.
        let inputs = [
            "",
            "data_x",
            "loop_\n_atom_site.Cartn_x\n",
            "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nC 1.0\n",
            "\u{0}\u{0}\u{0}",
            "data_x\nloop_\n_atom_site.type_symbol\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nC ( ) [",
        ];
        for input in inputs {
            let _ = parse_mmcif(input); // must not panic regardless of Ok/Err
        }
    }

    #[test]
    fn input_too_large_is_rejected() {
        let limits = MmcifParseLimits {
            max_input_bytes: 16,
            ..MmcifParseLimits::default()
        };
        let err = parse_mmcif_with_limits(FIXTURE, &limits).unwrap_err();
        assert_eq!(err, MmcifError::InputTooLarge { limit: 16 });
    }

    #[test]
    fn too_many_atoms_is_rejected() {
        let limits = MmcifParseLimits {
            max_atoms: 3,
            ..MmcifParseLimits::default()
        };
        let err = parse_mmcif_with_limits(FIXTURE, &limits).unwrap_err();
        assert_eq!(err, MmcifError::TooManyAtoms { limit: 3 });
    }

    #[test]
    fn line_too_long_is_rejected() {
        let long_line = "x".repeat(20_000);
        let cif = format!("data_x\n# {long_line}\n");
        let limits = MmcifParseLimits {
            max_line_len: 1000,
            ..MmcifParseLimits::default()
        };
        let err = parse_mmcif_with_limits(&cif, &limits).unwrap_err();
        assert!(matches!(err, MmcifError::LineTooLong { .. }));
    }

    #[test]
    fn quote_cif_value_quotes_values_with_whitespace() {
        assert_eq!(quote_cif_value("P 32 2 1"), "'P 32 2 1'");
        assert_eq!(quote_cif_value("A"), "A");
        assert_eq!(quote_cif_value(""), ".");
    }

    #[test]
    fn atom_name_with_quote_character_round_trips() {
        // Nucleic-acid atom names like O5' are common and contain a
        // syntactically significant quote character.
        let cif = "data_x\n\
            loop_\n\
            _atom_site.type_symbol\n_atom_site.label_atom_id\n\
            _atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
            O \"O5'\" 1.0 2.0 3.0\n";
        let r = parse_mmcif(cif).unwrap();
        assert_eq!(r.atoms[0].atom_name, "O5'");
        let written = write_mmcif(&r.atoms, None, None, "x");
        let r2 = parse_mmcif(&written).unwrap();
        assert_eq!(r2.atoms[0].atom_name, "O5'");
    }

    #[test]
    fn legacy_postfix_sign_formal_charge_is_accepted() {
        // Files converted from legacy PDB carry the old postfix-sign
        // notation from PDB columns 79-80 ("2+", "1-") instead of the
        // mmCIF-mandated plain signed integer.
        let cif = "data_x\n\
            loop_\n_atom_site.type_symbol\n_atom_site.pdbx_formal_charge\n\
            _atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
            Ca 2+ 0.0 0.0 0.0\n\
            Cl 1- 1.0 0.0 0.0\n";
        let r = parse_mmcif(cif).unwrap();
        assert_eq!(r.atoms[0].formal_charge, Some(2));
        assert_eq!(r.atoms[1].formal_charge, Some(-1));
    }

    #[test]
    fn label_seq_id_is_not_fabricated_for_hetatm_on_write() {
        // label_seq_id is "." (None) for the ZN HETATM in FIXTURE even
        // though auth_seq_id (-> res_seq) is a real number (101) -- write
        // must not turn that "not applicable" into a fabricated integer.
        let r = parse_mmcif(FIXTURE).unwrap();
        let zn = r.atoms.iter().find(|a| a.group_pdb == "HETATM").unwrap();
        assert_eq!(zn.label_seq_id, None);
        assert_eq!(zn.res_seq, 101);

        let written = write_mmcif(&r.atoms, r.cell.as_ref(), r.space_group.as_deref(), "TEST");
        let r2 = parse_mmcif(&written).unwrap();
        let zn2 = r2.atoms.iter().find(|a| a.group_pdb == "HETATM").unwrap();
        assert_eq!(zn2.label_seq_id, None);
        assert_eq!(zn2.res_seq, 101);
    }

    #[test]
    fn both_auth_and_label_as_placeholder_does_not_leak_the_dot_into_stored_data() {
        // Previously, when auth_asym_id was itself "." (not just absent)
        // and label_asym_id was also ".", the fallback used label_asym_id
        // unconditionally and stored the literal string "." as chain_id.
        let cif = "data_x\n\
            loop_\n\
            _atom_site.type_symbol\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n\
            _atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
            C . . 0.0 0.0 0.0\n";
        let r = parse_mmcif(cif).unwrap();
        assert_eq!(
            r.atoms[0].chain_id, "",
            "must not store the literal '.' placeholder"
        );
    }

    #[test]
    fn oracle_two_atom_water_like_fragment() {
        // Values independently chosen from the wwPDB dictionary's example
        // atom_site rows (not copied from a specific real PDB entry) --
        // exact known values, asserted precisely rather than just "parses".
        let cif = "data_water\n\
            loop_\n\
            _atom_site.group_PDB\n\
            _atom_site.id\n\
            _atom_site.type_symbol\n\
            _atom_site.label_atom_id\n\
            _atom_site.label_alt_id\n\
            _atom_site.label_comp_id\n\
            _atom_site.label_asym_id\n\
            _atom_site.label_seq_id\n\
            _atom_site.Cartn_x\n\
            _atom_site.Cartn_y\n\
            _atom_site.Cartn_z\n\
            _atom_site.occupancy\n\
            _atom_site.B_iso_or_equiv\n\
            HETATM 1 O O . HOH A . 0.000 0.000 0.000 1.00 30.00\n";
        let r = parse_mmcif(cif).unwrap();
        assert_eq!(r.atoms.len(), 1);
        let a = &r.atoms[0];
        assert_eq!(a.group_pdb, "HETATM");
        assert_eq!(a.element, Element::O);
        assert_eq!(a.res_name, "HOH");
        assert_eq!(a.chain_id, "A");
        assert_eq!(a.res_seq, 0); // no auth_seq_id/label_seq_id given -> 0
        assert_eq!((a.x, a.y, a.z), (0.0, 0.0, 0.0));
        assert!((a.occupancy - 1.0).abs() < 1e-9);
        assert!((a.b_iso - 30.0).abs() < 1e-9);
    }
}
