//! CIF (Crystallographic Information File) parser and writer.
//!
//! ## Scope
//!
//! Extracts **atomic positions** and converts fractional coordinates to
//! orthogonal Ångström coordinates.  Symmetry expansion (space group /
//! symmetry operations) is **not** performed — only atoms listed in the
//! `_atom_site_*` loop are returned (effectively P1 treatment).
//!
//! ## CIF structure
//!
//! ```text
//! data_NaCl
//! _cell_length_a   5.6402
//! _cell_length_b   5.6402
//! _cell_length_c   5.6402
//! _cell_angle_alpha  90.000
//! _cell_angle_beta   90.000
//! _cell_angle_gamma  90.000
//!
//! loop_
//! _atom_site_label
//! _atom_site_type_symbol
//! _atom_site_fract_x
//! _atom_site_fract_y
//! _atom_site_fract_z
//! Na1  Na  0.00000  0.00000  0.00000
//! Cl1  Cl  0.50000  0.50000  0.50000
//! ```
//!
//! ## Fractional → Orthogonal conversion (IUCr convention)
//!
//! ```text
//! X = a·fx + b·cos(γ)·fy + c·cos(β)·fz
//! Y =        b·sin(γ)·fy + c·(cos(α)−cos(β)cos(γ))/sin(γ)·fz
//! Z =                      (V/(a·b·sin(γ)))·fz
//! ```

use chematic_core::{Atom, AtomIdx, Element, Molecule, MoleculeBuilder};

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Unit cell parameters (lengths in Å, angles in degrees).
#[derive(Debug, Clone, PartialEq)]
pub struct UnitCell {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub alpha: f64,
    pub beta: f64,
    pub gamma: f64,
}

impl UnitCell {
    /// Volume in Å³.
    pub fn volume(&self) -> f64 {
        let (ca, cb, cg) = (
            self.alpha.to_radians().cos(),
            self.beta.to_radians().cos(),
            self.gamma.to_radians().cos(),
        );
        self.a * self.b * self.c * (1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg).sqrt()
    }

    /// Convert fractional coordinates to orthogonal Å.
    pub fn frac_to_cart(&self, fx: f64, fy: f64, fz: f64) -> (f64, f64, f64) {
        let ca = self.alpha.to_radians().cos();
        let cb = self.beta.to_radians().cos();
        let sg = self.gamma.to_radians().sin();
        let cg = self.gamma.to_radians().cos();
        let x = self.a * fx + self.b * cg * fy + self.c * cb * fz;
        let y = self.b * sg * fy + self.c * ((ca - cb * cg) / sg) * fz;
        let z = (self.volume() / (self.a * self.b * sg)) * fz;
        (x, y, z)
    }
}

impl Default for UnitCell {
    fn default() -> Self {
        Self {
            a: 1.0,
            b: 1.0,
            c: 1.0,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
        }
    }
}

/// Result of parsing a CIF file.
pub struct CifResult {
    /// Molecular topology (atoms; no bonds).
    pub mol: Molecule,
    /// Atomic coordinates in orthogonal Ångströms as `(x, y, z)` tuples.
    pub coords: Vec<(f64, f64, f64)>,
    /// Unit cell parameters, if present in the file.
    pub cell: Option<UnitCell>,
}

// ---------------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------------

/// Errors that can occur when parsing CIF files.
#[derive(Debug, Clone, PartialEq)]
pub enum CifError {
    /// No `_atom_site_*` loop with coordinates was found.
    NoAtomSiteLoop,
    /// The atom site loop lacked at least one coordinate column.
    MissingCoordinateColumns,
    /// An element symbol or label could not be resolved.
    UnknownElement(String),
    /// A coordinate value could not be parsed as a float.
    InvalidCoordinate(String),
    /// The unit cell has degenerate angles (e.g. γ = 0° or 180°) that make
    /// the fractional → Cartesian transformation undefined.
    InvalidCellParameters(String),
    /// The atom site loop uses fractional coordinates but the CIF contains no
    /// `_cell_length_*` / `_cell_angle_*` parameters to convert them.
    MissingCellParameters,
}

impl core::fmt::Display for CifError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NoAtomSiteLoop => write!(f, "no _atom_site_* loop found in CIF"),
            Self::MissingCoordinateColumns => {
                write!(f, "atom_site loop missing fract_x/y/z columns")
            }
            Self::UnknownElement(s) => write!(f, "unknown element '{s}' in CIF"),
            Self::InvalidCoordinate(s) => write!(f, "invalid coordinate '{s}' in CIF"),
            Self::InvalidCellParameters(s) => write!(f, "invalid cell parameters: {s}"),
            Self::MissingCellParameters => write!(
                f,
                "fractional coordinates present but no _cell_length_*/_cell_angle_* parameters found"
            ),
        }
    }
}

impl std::error::Error for CifError {}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Parse a CIF file and return atoms with orthogonal coordinates.
///
/// Only the first `data_` block is parsed.  Symmetry expansion is **not**
/// performed; returned atoms are exactly those listed in the `_atom_site_*`
/// loop.
/// Strip a CIF comment from one line, respecting single- and double-quoted
/// strings (a `#` inside quotes is not a comment delimiter).
fn strip_cif_comment(line: &str) -> &str {
    let mut in_single = false;
    let mut in_double = false;
    for (i, c) in line.char_indices() {
        match c {
            '\'' if !in_double => in_single = !in_single,
            '"' if !in_single => in_double = !in_double,
            '#' if !in_single && !in_double => return &line[..i],
            _ => {}
        }
    }
    line
}

pub fn parse_cif(input: &str) -> Result<CifResult, CifError> {
    // Strip CIF comments (# to end of line), but not '#' inside quoted strings.
    let clean: String = input
        .lines()
        .map(strip_cif_comment)
        .collect::<Vec<_>>()
        .join("\n");

    let tokens = tokenize_cif(&clean);

    // --- Extract cell parameters ---
    let mut cell = UnitCell::default();
    let mut has_cell = false;
    let mut i = 0;
    while i + 1 < tokens.len() {
        match tokens[i].to_ascii_lowercase().as_str() {
            "_cell_length_a" => {
                cell.a = parse_esd(&tokens[i + 1]).unwrap_or(cell.a);
                has_cell = true;
            }
            "_cell_length_b" => {
                cell.b = parse_esd(&tokens[i + 1]).unwrap_or(cell.b);
            }
            "_cell_length_c" => {
                cell.c = parse_esd(&tokens[i + 1]).unwrap_or(cell.c);
            }
            "_cell_angle_alpha" => {
                cell.alpha = parse_esd(&tokens[i + 1]).unwrap_or(cell.alpha);
            }
            "_cell_angle_beta" => {
                cell.beta = parse_esd(&tokens[i + 1]).unwrap_or(cell.beta);
            }
            "_cell_angle_gamma" => {
                cell.gamma = parse_esd(&tokens[i + 1]).unwrap_or(cell.gamma);
            }
            _ => {}
        }
        i += 1;
    }

    // --- Find _atom_site loop ---
    i = 0;
    let mut atom_site_data_start: Option<usize> = None;
    let mut col_headers: Vec<String> = Vec::new();

    while i < tokens.len() {
        if tokens[i].as_str() == "loop_" {
            let mut j = i + 1;
            let mut headers: Vec<String> = Vec::new();
            while j < tokens.len() && tokens[j].starts_with('_') {
                headers.push(tokens[j].to_ascii_lowercase());
                j += 1;
            }
            if headers.iter().any(|h| h.starts_with("_atom_site_")) {
                col_headers = headers;
                atom_site_data_start = Some(j);
                break;
            }
            i = j;
        } else {
            i += 1;
        }
    }

    let data_start = atom_site_data_start.ok_or(CifError::NoAtomSiteLoop)?;
    let ncols = col_headers.len();
    if ncols == 0 {
        return Err(CifError::NoAtomSiteLoop);
    }

    let col = |name: &str| -> Option<usize> { col_headers.iter().position(|h| h.as_str() == name) };

    let col_type = col("_atom_site_type_symbol");
    let col_label = col("_atom_site_label");
    // Prefer fractional; fall back to Cartesian.
    let use_cartesian = col("_atom_site_fract_x").is_none();
    let col_x = col("_atom_site_fract_x").or_else(|| col("_atom_site_cartn_x"));
    let col_y = col("_atom_site_fract_y").or_else(|| col("_atom_site_cartn_y"));
    let col_z = col("_atom_site_fract_z").or_else(|| col("_atom_site_cartn_z"));

    let (col_x, col_y, col_z) = match (col_x, col_y, col_z) {
        (Some(x), Some(y), Some(z)) => (x, y, z),
        _ => return Err(CifError::MissingCoordinateColumns),
    };

    // Validate that fractional→Cartesian conversion is well-defined.
    if !use_cartesian {
        if !has_cell {
            return Err(CifError::MissingCellParameters);
        }
        let sg = cell.gamma.to_radians().sin();
        if sg.abs() < 1e-10 {
            return Err(CifError::InvalidCellParameters(format!(
                "_cell_angle_gamma = {} makes sin(γ) ≈ 0, transformation undefined",
                cell.gamma
            )));
        }
        // Also guard against a non-physical cell volume (≤ 0 or NaN).
        let vol = cell.volume();
        if !vol.is_finite() || vol <= 0.0 {
            return Err(CifError::InvalidCellParameters(format!(
                "unit cell volume is non-positive ({vol:.6} Å³); check _cell_angle_* values"
            )));
        }
    }

    let mut builder = MoleculeBuilder::new();
    let mut coords: Vec<(f64, f64, f64)> = Vec::new();
    let data_tokens = &tokens[data_start..];

    let mut row = 0;
    while (row + 1) * ncols <= data_tokens.len() {
        let base = row * ncols;
        let tok = &data_tokens[base..base + ncols];

        // Stop at next loop_ or data_ block.
        if tok[0].as_str() == "loop_" || tok[0].starts_with("data_") {
            break;
        }

        // Resolve element.
        let elem_raw: &str = col_type
            .and_then(|c| tok.get(c))
            .map(|s| s.as_str())
            .or_else(|| col_label.and_then(|c| tok.get(c)).map(|s| s.as_str()))
            .unwrap_or("X");
        // Strip trailing digits and oxidation-state signs from labels like
        // "Na1", "Cu2+", "Fe3+", "O2-".
        let elem_str =
            elem_raw.trim_end_matches(|c: char| c.is_ascii_digit() || c == '+' || c == '-');
        let elem = Element::from_symbol(elem_str)
            .ok_or_else(|| CifError::UnknownElement(elem_str.to_string()))?;

        let parse_coord = |s: &str| -> Result<f64, CifError> {
            strip_esd(s)
                .parse::<f64>()
                .map_err(|_| CifError::InvalidCoordinate(s.to_string()))
        };
        let fx = parse_coord(&tok[col_x])?;
        let fy = parse_coord(&tok[col_y])?;
        let fz = parse_coord(&tok[col_z])?;

        let (x, y, z) = if use_cartesian {
            (fx, fy, fz)
        } else {
            cell.frac_to_cart(fx, fy, fz)
        };

        builder.add_atom(Atom::new(elem));
        coords.push((x, y, z));
        row += 1;
    }

    if coords.is_empty() {
        return Err(CifError::NoAtomSiteLoop);
    }

    Ok(CifResult {
        mol: builder.build(),
        coords,
        cell: if has_cell { Some(cell) } else { None },
    })
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// Write a minimal CIF file (P1, no symmetry operators).
///
/// When `cell` is `None`, fractional coordinates are written as-is
/// (the coordinates are treated as already orthogonal with a=b=c=1).
pub fn write_cif(mol: &Molecule, coords: &[(f64, f64, f64)], cell: Option<&UnitCell>) -> String {
    let mut out = String::from("data_chematic\n\n");

    if let Some(c) = cell {
        out.push_str(&format!(
            "_cell_length_a   {:.4}\n\
             _cell_length_b   {:.4}\n\
             _cell_length_c   {:.4}\n\
             _cell_angle_alpha  {:.3}\n\
             _cell_angle_beta   {:.3}\n\
             _cell_angle_gamma  {:.3}\n\n\
             _symmetry_space_group_name_H-M  'P 1'\n\n",
            c.a, c.b, c.c, c.alpha, c.beta, c.gamma
        ));
    }

    out.push_str(
        "loop_\n\
         _atom_site_label\n\
         _atom_site_type_symbol\n\
         _atom_site_fract_x\n\
         _atom_site_fract_y\n\
         _atom_site_fract_z\n",
    );

    let fallback_cell = UnitCell::default();
    let cell_ref = cell.unwrap_or(&fallback_cell);

    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let sym = mol.atom(idx).element.symbol();
        let (x, y, z) = coords.get(i).copied().unwrap_or((0.0, 0.0, 0.0));
        let (fx, fy, fz) = cart_to_frac(cell_ref, x, y, z);
        out.push_str(&format!(
            "{sym}{} {sym}  {fx:.5}  {fy:.5}  {fz:.5}\n",
            i + 1
        ));
    }
    out
}

fn cart_to_frac(cell: &UnitCell, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    // Exact inverse for orthogonal cells; approximate for triclinic.
    let sg = cell.gamma.to_radians().sin();
    let cg = cell.gamma.to_radians().cos();
    let cb = cell.beta.to_radians().cos();
    let ca = cell.alpha.to_radians().cos();
    let fz = z * cell.a * cell.b * sg / cell.volume();
    let fy = (y - cell.c * ((ca - cb * cg) / sg) * fz) / (cell.b * sg);
    let fx = (x - cell.b * cg * fy - cell.c * cb * fz) / cell.a;
    (fx, fy, fz)
}

// ---------------------------------------------------------------------------
// CIF tokenizer
// ---------------------------------------------------------------------------

fn tokenize_cif(input: &str) -> Vec<String> {
    let mut tokens: Vec<String> = Vec::new();
    let chars: Vec<char> = input.chars().collect();
    let len = chars.len();
    let mut i = 0;

    while i < len {
        let ch = chars[i];
        // Skip whitespace.
        if ch.is_whitespace() {
            i += 1;
            continue;
        }

        // Semicolon text block (starts at column 0 on a new line).
        if ch == ';' && (i == 0 || chars[i - 1] == '\n') {
            i += 1;
            let start = i;
            while i < len {
                if chars[i] == '\n' && i + 1 < len && chars[i + 1] == ';' {
                    break;
                }
                i += 1;
            }
            // Skip the text content (we don't use it for atom sites).
            tokens.push(chars[start..i].iter().collect());
            i += 2; // skip '\n;'
            continue;
        }

        // Single-quoted string.
        if ch == '\'' {
            i += 1;
            let start = i;
            while i < len && chars[i] != '\'' {
                i += 1;
            }
            tokens.push(chars[start..i].iter().collect());
            if i < len {
                i += 1;
            }
            continue;
        }

        // Double-quoted string.
        if ch == '"' {
            i += 1;
            let start = i;
            while i < len && chars[i] != '"' {
                i += 1;
            }
            tokens.push(chars[start..i].iter().collect());
            if i < len {
                i += 1;
            }
            continue;
        }

        // Regular token.
        let start = i;
        while i < len && !chars[i].is_whitespace() {
            i += 1;
        }
        tokens.push(chars[start..i].iter().collect());
    }
    tokens
}

fn strip_esd(s: &str) -> &str {
    s.find('(').map_or(s, |pos| &s[..pos])
}

fn parse_esd(s: &str) -> Option<f64> {
    strip_esd(s).parse::<f64>().ok()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const NACL_CIF: &str = r#"data_NaCl
_cell_length_a   5.6402
_cell_length_b   5.6402
_cell_length_c   5.6402
_cell_angle_alpha  90.000
_cell_angle_beta   90.000
_cell_angle_gamma  90.000

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Na1  Na  0.00000  0.00000  0.00000
Cl1  Cl  0.50000  0.50000  0.50000
"#;

    #[test]
    fn parse_nacl_atom_count() {
        let r = parse_cif(NACL_CIF).unwrap();
        assert_eq!(r.mol.atom_count(), 2);
        assert_eq!(r.coords.len(), 2);
    }

    #[test]
    fn parse_nacl_cell_params() {
        let r = parse_cif(NACL_CIF).unwrap();
        let cell = r.cell.as_ref().unwrap();
        assert!((cell.a - 5.6402).abs() < 1e-4);
        assert!((cell.alpha - 90.0).abs() < 1e-4);
    }

    #[test]
    fn parse_nacl_na_at_origin() {
        let r = parse_cif(NACL_CIF).unwrap();
        let (x, y, z) = r.coords[0];
        assert!(
            x.abs() < 1e-6 && y.abs() < 1e-6 && z.abs() < 1e-6,
            "Na at origin: got ({x}, {y}, {z})"
        );
    }

    #[test]
    fn parse_nacl_cl_at_body_center() {
        let r = parse_cif(NACL_CIF).unwrap();
        let (x, y, z) = r.coords[1];
        let half = 5.6402 / 2.0;
        assert!((x - half).abs() < 1e-3, "Cl x: got {x}, expected ~{half}");
        assert!((y - half).abs() < 1e-3);
        assert!((z - half).abs() < 1e-3);
    }

    #[test]
    fn unit_cell_volume_cubic() {
        let cell = UnitCell {
            a: 5.0,
            b: 5.0,
            c: 5.0,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
        };
        assert!((cell.volume() - 125.0).abs() < 1e-6);
    }

    #[test]
    fn write_cif_roundtrip() {
        let r = parse_cif(NACL_CIF).unwrap();
        let out = write_cif(&r.mol, &r.coords, r.cell.as_ref());
        assert!(out.contains("data_chematic"));
        let r2 = parse_cif(&out).unwrap();
        assert_eq!(r2.mol.atom_count(), 2);
    }

    #[test]
    fn parse_esd_stripped() {
        let cif = "data_x\n_cell_length_a 5.64(2)\n_cell_length_b 5.64\n_cell_length_c 5.64\n\
                   _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                   loop_\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                   Na 0.0 0.0 0.0\n";
        let r = parse_cif(cif).unwrap();
        assert!((r.cell.unwrap().a - 5.64).abs() < 1e-3);
    }

    #[test]
    fn frac_to_cart_cubic_roundtrip() {
        let cell = UnitCell {
            a: 5.0,
            b: 5.0,
            c: 5.0,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
        };
        let (x, y, z) = cell.frac_to_cart(0.3, 0.4, 0.5);
        let (fx, fy, fz) = cart_to_frac(&cell, x, y, z);
        assert!((fx - 0.3).abs() < 1e-10);
        assert!((fy - 0.4).abs() < 1e-10);
        assert!((fz - 0.5).abs() < 1e-10);
    }
}
