//! Plain CIF atomic-position writers.
use super::{CifError, UnitCell};
use chematic_core::{AtomIdx, Molecule};

/// Write a minimal CIF file (P1, no symmetry operators).
///
/// Input coordinates are Cartesian Å. With a cell, they are converted to
/// fractional coordinates. Without a cell, they are emitted unchanged under
/// `_atom_site_Cartn_x/y/z`, without inventing a unit cell (Open Babel #2995).
///
/// This legacy writer pads missing coordinates with zero and does not validate
/// the cell or numbers. Prefer [`write_cif_checked`] for caller-supplied input.
pub fn write_cif(mol: &Molecule, coords: &[(f64, f64, f64)], cell: Option<&UnitCell>) -> String {
    write_cif_impl(mol, coords, cell, false)
}

/// Write a minimal P1 CIF after validating the cell and coordinates.
///
/// Requires exactly one finite Cartesian coordinate per atom. Rejects invalid
/// or numerically near-singular cells (see [`UnitCell::validate`]) and non-finite
/// results of the Cartesian-to-fractional conversion (Open Babel #2818).
/// Uses round-trip float formatting so a small valid cell is not rounded to
/// zero by the legacy writer's fixed decimal precision. No bonds, symmetry
/// expansion, occupancies or molecular metadata are serialized.
pub fn write_cif_checked(
    mol: &Molecule,
    coords: &[(f64, f64, f64)],
    cell: Option<&UnitCell>,
) -> Result<String, CifError> {
    if coords.len() != mol.atom_count() {
        return Err(CifError::CoordinateCountMismatch {
            expected: mol.atom_count(),
            actual: coords.len(),
        });
    }
    if let Some(cell) = cell {
        cell.validate()?;
    }
    for (i, &(x, y, z)) in coords.iter().enumerate() {
        if !x.is_finite() || !y.is_finite() || !z.is_finite() {
            return Err(CifError::InvalidCoordinate(format!(
                "atom {} has a non-finite Cartesian coordinate",
                i + 1
            )));
        }
        let (fx, fy, fz) = cell.map_or((x, y, z), |cell| cart_to_frac(cell, x, y, z));
        if !fx.is_finite() || !fy.is_finite() || !fz.is_finite() {
            return Err(CifError::InvalidCoordinate(format!(
                "atom {} overflows during Cartesian-to-fractional conversion",
                i + 1
            )));
        }
    }
    Ok(write_cif_impl(mol, coords, cell, true))
}

fn write_cif_impl(
    mol: &Molecule,
    coords: &[(f64, f64, f64)],
    cell: Option<&UnitCell>,
    round_trip: bool,
) -> String {
    let mut out = String::from("data_chematic\n\n");

    if let Some(c) = cell {
        let cell_text = if round_trip {
            format!(
                "_cell_length_a   {}\n_cell_length_b   {}\n_cell_length_c   {}\n\
                 _cell_angle_alpha  {}\n_cell_angle_beta   {}\n_cell_angle_gamma  {}\n\n\
                 _symmetry_space_group_name_H-M  'P 1'\n\n",
                c.a, c.b, c.c, c.alpha, c.beta, c.gamma
            )
        } else {
            format!(
                "_cell_length_a   {:.4}\n\
             _cell_length_b   {:.4}\n\
             _cell_length_c   {:.4}\n\
             _cell_angle_alpha  {:.3}\n\
             _cell_angle_beta   {:.3}\n\
             _cell_angle_gamma  {:.3}\n\n\
             _symmetry_space_group_name_H-M  'P 1'\n\n",
                c.a, c.b, c.c, c.alpha, c.beta, c.gamma
            )
        };
        out.push_str(&cell_text);
    }

    out.push_str(
        "loop_\n\
         _atom_site_label\n\
         _atom_site_type_symbol\n",
    );
    let family = if cell.is_some() { "fract" } else { "Cartn" };
    out.push_str(&format!(
        "_atom_site_{family}_x\n_atom_site_{family}_y\n_atom_site_{family}_z\n"
    ));

    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let sym = mol.atom(idx).element.symbol();
        let (x, y, z) = coords.get(i).copied().unwrap_or((0.0, 0.0, 0.0));
        let (fx, fy, fz) = cell.map_or((x, y, z), |cell| cart_to_frac(cell, x, y, z));
        if round_trip {
            out.push_str(&format!("{sym}{} {sym}  {fx}  {fy}  {fz}\n", i + 1));
        } else {
            out.push_str(&format!(
                "{sym}{} {sym}  {fx:.5}  {fy:.5}  {fz:.5}\n",
                i + 1
            ));
        }
    }
    out
}

pub(super) fn cart_to_frac(cell: &UnitCell, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    // Inverse of UnitCell::frac_to_cart, including triclinic cells.
    let sg = cell.gamma.to_radians().sin();
    let cg = cell.gamma.to_radians().cos();
    let cb = cell.beta.to_radians().cos();
    let ca = cell.alpha.to_radians().cos();
    let fz = z * cell.a * cell.b * sg / cell.volume();
    let fy = (y - cell.c * ((ca - cb * cg) / sg) * fz) / (cell.b * sg);
    let fx = (x - cell.b * cg * fy - cell.c * cb * fz) / cell.a;
    (fx, fy, fz)
}
