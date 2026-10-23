//! Additive CML Cartesian reader. Legacy graph/2D readers keep their behavior.
use super::*;

/// A source-ordered graph and optional complete Cartesian Angstrom coordinates.
pub type CmlWith3d = (Molecule, Option<Vec<[f64; 3]>>);

fn malformed(detail: impl Into<String>) -> CmlError {
    CmlError::MalformedXml(detail.into())
}

fn cartesian(attrs: &HashMap<String, String>) -> Result<Option<[f64; 3]>, CmlError> {
    if ["xyz3", "xFract", "yFract", "zFract", "xyzFract"]
        .iter()
        .any(|key| attrs.contains_key(*key))
    {
        return Err(malformed(
            "packed/fractional CML coordinates are outside the Cartesian triplet profile",
        ));
    }
    let fields = ["x3", "y3", "z3"];
    let count = fields
        .iter()
        .filter(|key| attrs.contains_key(**key))
        .count();
    if count == 0 {
        return Ok(None);
    }
    if count != 3 {
        return Err(malformed(
            "CML 3D coordinates require x3, y3 and z3 on every atom",
        ));
    }
    let mut point = [0.0; 3];
    for (slot, key) in point.iter_mut().zip(fields) {
        let value = &attrs[key];
        *slot = value
            .trim()
            .parse::<f64>()
            .map_err(|_| malformed(format!("invalid CML {key}: {value}")))?;
        if !slot.is_finite() {
            return Err(CmlError::NonFiniteCoordinate(value.clone()));
        }
    }
    Ok(Some(point))
}

/// Read one nonempty CML molecule with optional Cartesian coordinates in Angstrom.
///
/// Returns source atom order and `None` when all atoms lack x3/y3/z3. Otherwise
/// every atom must have a complete finite triplet; no zero filling or geometry
/// generation occurs. If both 2D and 3D are supplied, only 3D is returned.
/// Requires unique nonempty atom IDs and exactly one molecule. CML namespace
/// prefixes and legacy unnamespaced input are supported; foreign subtrees and
/// attributes never become graph data. XML is validated without DTD processing.
/// XML schema, opaque metadata and coordinate stereo are not interpreted.
/// Existing graph/2D readers are unchanged. Default limits are CmlParseLimits.
pub fn parse_cml_with_3d(input: &str) -> Result<CmlWith3d, CmlError> {
    parse_cml_with_3d_with_limits(input, &CmlParseLimits::default())
}

/// The byte/physical-line bounds are checked before XML normalization; normalized
/// lines/elements and atom/bond caps also apply before graph construction.
pub fn parse_cml_with_3d_with_limits(
    input: &str,
    limits: &CmlParseLimits,
) -> Result<CmlWith3d, CmlError> {
    cml_lines(input, limits)?;
    // Retain the existing normalized line budget independently of XML layout.
    let normalized = split_cml_tags(input);
    cml_lines(
        &normalized,
        &CmlParseLimits {
            max_input_bytes: normalized.len(),
            ..*limits
        },
    )?;
    let (source_atoms, source_bonds) = super::geometry_xml::records(input, limits)?;
    let mut ids = HashMap::new();
    let mut atoms = Vec::with_capacity(source_atoms.len());
    let mut points = Vec::with_capacity(source_atoms.len());
    let mut present = 0;
    for attrs in &source_atoms {
        let id = attrs
            .get("id")
            .filter(|id| !id.trim().is_empty())
            .ok_or_else(|| malformed("CML 3D reader requires nonempty atom IDs"))?;
        if ids.insert(id.clone(), atoms.len()).is_some() {
            return Err(malformed(format!("duplicate CML atom ID: {id}")));
        }
        if !attrs.contains_key("elementType") {
            return Err(malformed("CML atom requires elementType"));
        }
        let point = cartesian(attrs)?;
        present += usize::from(point.is_some());
        points.push(point);
        atoms.push(parse_cml_atom(attrs)?);
    }
    if present != 0 && present != points.len() {
        return Err(malformed(
            "CML 3D coordinates must be present on every atom or absent on all",
        ));
    }
    let mut bonds = Vec::with_capacity(source_bonds.len());
    for attrs in &source_bonds {
        let refs = attrs.get("atomRefs2").map(String::as_str).unwrap_or("");
        let endpoints = refs.split_whitespace().collect::<Vec<_>>();
        if endpoints.len() != 2 || endpoints[0] == endpoints[1] {
            return Err(CmlError::InvalidAtomRefs2(refs.into()));
        }
        if !attrs.contains_key("order") {
            return Err(malformed("CML bond requires order"));
        }
        bonds.push(parse_cml_bond(attrs, &ids)?);
    }
    let (mol, _) = build_cml(atoms, bonds)?;
    if mol.atom_count() == 0 {
        return Err(CmlError::EmptyMolecule);
    }
    let coords = (present != 0).then(|| points.into_iter().map(Option::unwrap).collect());
    Ok((mol, coords))
}
