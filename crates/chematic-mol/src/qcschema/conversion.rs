use super::*;

// ─── Conversion to/from chematic_core::Molecule ─────────────────────────────

/// CODATA 2018 Bohr radius in Angstrom (1 a0 = 0.529177210903 Å), matching
/// the constant set `qcelemental` uses as of this writing. QCSchema
/// `geometry` is always in Bohr; `chematic_core::Coords3D` is always in
/// Angstrom (see `chematic_core::coords3d` module docs) -- this constant is
/// the single conversion point between the two.
pub const BOHR_TO_ANGSTROM: f64 = 0.529177210903;

/// Errors converting between [`QcMolecule`] and `chematic_core`'s own
/// molecule/coordinate types.
#[derive(Debug, Clone, PartialEq)]
pub enum QcConvertError {
    UnknownElement(String),
    InvalidBondIndex {
        index: usize,
        atom_count: usize,
    },
    /// `coords.atom_count() != mol.atom_count()` on the chematic -> QCSchema
    /// direction.
    AtomCountMismatch {
        molecule: usize,
        coords: usize,
    },
}

impl std::fmt::Display for QcConvertError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnknownElement(s) => {
                write!(f, "QCSchema conversion: unknown element symbol '{s}'")
            }
            Self::InvalidBondIndex { index, atom_count } => write!(
                f,
                "QCSchema conversion: connectivity references atom index {index} but there are only {atom_count} atoms"
            ),
            Self::AtomCountMismatch { molecule, coords } => write!(
                f,
                "QCSchema conversion: molecule has {molecule} atoms but coords has {coords}"
            ),
        }
    }
}

impl std::error::Error for QcConvertError {}

/// Bundle returned by [`qc_molecule_to_chematic`]: `chematic_core::Molecule`
/// has no molecule-level charge/multiplicity field and no built-in
/// coordinate storage of its own (coordinates live in the separate
/// [`Coords3D`] type, same convention as this crate's XYZ reader), so both
/// are threaded through here rather than silently dropped.
// `chematic_core::Molecule` implements neither `Debug` nor `PartialEq`, so
// this bundle can't derive them either.
#[derive(Clone)]
pub struct ChematicMoleculeView {
    pub molecule: Molecule,
    /// Angstrom (converted from the source's Bohr geometry).
    pub coords: Coords3D,
    pub molecular_charge: f64,
    pub molecular_multiplicity: i64,
}

fn qc_bond_order_to_chematic(v: f64) -> BondOrder {
    if (v - 1.5).abs() < 1e-6 {
        return BondOrder::Aromatic;
    }
    match v.round() as i64 {
        0 => BondOrder::Zero,
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        4 => BondOrder::Quadruple,
        _ => BondOrder::Single,
    }
}

fn chematic_bond_order_to_qc(o: BondOrder) -> f64 {
    match o {
        BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => 1.0,
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Quadruple => 4.0,
        BondOrder::Aromatic => 1.5,
        BondOrder::Zero => 0.0,
        BondOrder::QueryAny
        | BondOrder::QuerySingleOrDouble
        | BondOrder::QuerySingleOrAromatic
        | BondOrder::QueryDoubleOrAromatic => 1.0,
    }
}

/// Convert a QCSchema [`QcMolecule`] into chematic's own graph + coordinate
/// representation.
///
/// **Gained:** nothing beyond what QCSchema already carries.
/// **Lost:** `fix_com`/`fix_orientation`/`fragments`/`fragment_charges`/
/// `fragment_multiplicities`/`real` (ghost-atom flag)/`masses` (explicit
/// isotope masses)/`atom_labels`/`name`/`comment`/`fix_symmetry`/
/// `provenance` -- `chematic_core::Molecule` and [`Coords3D`] have no
/// equivalent slots for any of these; call sites that need them should keep
/// the original [`QcMolecule`] around. `connectivity` bond orders (plain
/// floats, e.g. `1.0`, `2.0`, `1.5`) are mapped onto the nearest
/// [`BondOrder`]; QCSchema has no stereo/aromaticity-flag concept distinct
/// from that numeric order. `mass_numbers` (when present and not `-1`) is
/// mapped onto `Atom::isotope`.
pub fn qc_molecule_to_chematic(qc: &QcMolecule) -> Result<ChematicMoleculeView, QcConvertError> {
    let mut builder = MoleculeBuilder::new();
    let mut coords = Coords3D::new_zeroed(qc.symbols.len());

    for (i, sym) in qc.symbols.iter().enumerate() {
        let element =
            Element::from_symbol(sym).ok_or_else(|| QcConvertError::UnknownElement(sym.clone()))?;
        let mut atom = Atom::new(element);
        if let Some(mn) = qc.mass_numbers.as_ref().and_then(|v| v.get(i))
            && *mn >= 0
        {
            atom.isotope = Some(*mn as u16);
        }
        let idx = builder.add_atom(atom);
        let p = Point3::new(
            qc.geometry[i * 3] * BOHR_TO_ANGSTROM,
            qc.geometry[i * 3 + 1] * BOHR_TO_ANGSTROM,
            qc.geometry[i * 3 + 2] * BOHR_TO_ANGSTROM,
        );
        coords.set(idx, p);
    }

    if let Some(conn) = &qc.connectivity {
        for (a, b, order) in conn {
            if *a >= qc.symbols.len() || *b >= qc.symbols.len() {
                return Err(QcConvertError::InvalidBondIndex {
                    index: (*a).max(*b),
                    atom_count: qc.symbols.len(),
                });
            }
            let bo = qc_bond_order_to_chematic(*order);
            // Length/duplicate/self-bond invariants were already enforced by
            // `parse_qc_molecule_value`; `add_bond` cannot fail here.
            let _ = builder.add_bond(AtomIdx(*a as u32), AtomIdx(*b as u32), bo);
        }
    }

    Ok(ChematicMoleculeView {
        molecule: builder.build(),
        coords,
        molecular_charge: qc.molecular_charge,
        molecular_multiplicity: qc.molecular_multiplicity,
    })
}

/// Convert a chematic `Molecule` + its `Coords3D` (Angstrom) + molecule-level
/// charge/multiplicity into a QCSchema [`QcMolecule`] (Bohr).
///
/// **Gained:** nothing. **Lost:** per-atom formal charge, isotope-vs-mass-
/// number distinction (only the isotope's mass number survives, as
/// `mass_numbers`; `masses` -- fractional atomic weight -- is left `None`,
/// letting a QCSchema consumer fall back to canonical weights), aromaticity,
/// stereochemistry (`Chirality`/CIP code), and wildcard/atom-map annotations
/// -- QCSchema's `Molecule` is a geometry+connectivity schema, not a full
/// cheminformatics model, and has no field for any of these.
pub fn chematic_to_qc_molecule(
    mol: &Molecule,
    coords: &Coords3D,
    molecular_charge: f64,
    molecular_multiplicity: i64,
) -> Result<QcMolecule, QcConvertError> {
    if coords.atom_count() != mol.atom_count() {
        return Err(QcConvertError::AtomCountMismatch {
            molecule: mol.atom_count(),
            coords: coords.atom_count(),
        });
    }

    let mut symbols = Vec::with_capacity(mol.atom_count());
    let mut geometry = Vec::with_capacity(mol.atom_count() * 3);
    let mut mass_numbers = Vec::with_capacity(mol.atom_count());
    let mut any_isotope = false;

    for (idx, atom) in mol.atoms() {
        symbols.push(atom.element.symbol().to_string());
        let p = coords.get(idx);
        geometry.push(p.x / BOHR_TO_ANGSTROM);
        geometry.push(p.y / BOHR_TO_ANGSTROM);
        geometry.push(p.z / BOHR_TO_ANGSTROM);
        match atom.isotope {
            Some(iso) => {
                mass_numbers.push(iso as i64);
                any_isotope = true;
            }
            None => mass_numbers.push(-1),
        }
    }

    let connectivity = if mol.bond_count() > 0 {
        Some(
            mol.bonds()
                .map(|(_, b)| {
                    (
                        b.atom1.0 as usize,
                        b.atom2.0 as usize,
                        chematic_bond_order_to_qc(b.order),
                    )
                })
                .collect(),
        )
    } else {
        None
    };

    Ok(QcMolecule {
        schema_name: "qcschema_molecule".to_string(),
        schema_version: 1,
        symbols,
        geometry,
        molecular_charge,
        molecular_multiplicity,
        fix_com: false,
        fix_orientation: false,
        masses: None,
        real: None,
        atomic_numbers: None,
        mass_numbers: if any_isotope {
            Some(mass_numbers)
        } else {
            None
        },
        atom_labels: None,
        name: None,
        comment: None,
        connectivity,
        fragments: None,
        fragment_charges: None,
        fragment_multiplicities: None,
        fix_symmetry: None,
        provenance: None,
        id: None,
        extras: JsonObject::new(),
        unknown_fields: JsonObject::new(),
    })
}
