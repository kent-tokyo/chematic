//! Additive CML Cartesian reader. Legacy graph/2D readers keep their behavior.
use super::*;
use std::collections::HashSet;

/// A source-ordered graph and optional complete Cartesian Angstrom coordinates.
pub type CmlWith3d = (Molecule, Option<Vec<[f64; 3]>>);

fn malformed(detail: impl Into<String>) -> CmlError {
    CmlError::MalformedXml(detail.into())
}

// Validate attribute boundaries and uniqueness before the legacy value decoder
// can overwrite duplicates. This is a local XML subset, not schema validation.
fn atom_attrs(line: &str) -> Result<HashMap<String, String>, CmlError> {
    let mut rest = line
        .find(|c: char| c.is_ascii_whitespace() || c == '/' || c == '>')
        .map(|end| &line[end..])
        .ok_or_else(|| malformed("invalid CML atom tag"))?;
    let mut names = HashSet::new();
    loop {
        rest = rest.trim_start();
        if rest == ">" || rest == "/>" {
            break;
        }
        let end = rest
            .find(|c: char| c.is_ascii_whitespace() || c == '=')
            .ok_or_else(|| malformed("invalid CML atom attribute"))?;
        let name = &rest[..end];
        if name.is_empty() || !names.insert(name) {
            return Err(malformed(format!(
                "duplicate or empty CML atom attribute: {name}"
            )));
        }
        rest = rest[end..].trim_start();
        rest = rest
            .strip_prefix('=')
            .ok_or_else(|| malformed("missing attribute ="))?
            .trim_start();
        let quote = rest
            .chars()
            .next()
            .filter(|c| matches!(c, '\'' | '"'))
            .ok_or_else(|| malformed("unquoted CML atom attribute"))?;
        rest = &rest[1..];
        let end = rest
            .find(quote)
            .ok_or_else(|| malformed("unclosed CML atom attribute"))?;
        if rest[..end].contains('<') {
            return Err(malformed("raw < in CML atom attribute"));
        }
        rest = &rest[end + 1..];
        if !rest.is_empty()
            && !rest.starts_with(|c: char| c.is_ascii_whitespace() || c == '/' || c == '>')
        {
            return Err(malformed("missing space between CML atom attributes"));
        }
    }
    Ok(parse_xml_attrs(line))
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
/// Requires unique nonempty atom IDs and exactly one molecule. Compact, LF and
/// CRLF forms are supported by the existing lightweight structural parser.
/// Namespace/schema, opaque metadata and coordinate stereo are not interpreted.
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
    validate_strict_cml_structure(input)?;
    let normalized = split_cml_tags(input);
    // The input-byte cap describes source bytes. Inserted tag separators must
    // not make an input at that cap fail; source size already bounds expansion.
    let normalized_limits = CmlParseLimits {
        max_input_bytes: normalized.len(),
        ..*limits
    };
    let lines = cml_lines(&normalized, &normalized_limits)?;
    let mut ids = HashSet::new();
    let mut points = Vec::new();
    let mut molecules = 0;
    let mut present = 0;
    for line in lines.iter().map(|line| line.trim()) {
        if is_element_tag(line, "molecule") {
            molecules += 1;
        }
        if !is_element_tag(line, "atom") {
            continue;
        }
        if points.len() >= limits.max_atoms {
            return Err(CmlError::ResourceLimit {
                resource: "atom elements",
                actual: points.len() + 1,
                limit: limits.max_atoms,
            });
        }
        let attrs = atom_attrs(line)?;
        let id = attrs
            .get("id")
            .filter(|id| !id.trim().is_empty())
            .ok_or_else(|| malformed("CML 3D reader requires nonempty atom IDs"))?;
        if !ids.insert(id.clone()) {
            return Err(malformed(format!("duplicate CML atom ID: {id}")));
        }
        let point = cartesian(&attrs)?;
        present += usize::from(point.is_some());
        points.push(point);
    }
    if molecules != 1 {
        return Err(malformed("CML 3D reader requires exactly one molecule"));
    }
    if present != 0 && present != points.len() {
        return Err(malformed(
            "CML 3D coordinates must be present on every atom or absent on all",
        ));
    }
    let (atoms, bonds) = collect_cml_data(lines, limits)?;
    let (mol, _) = build_cml(atoms, bonds)?;
    if mol.atom_count() == 0 {
        return Err(CmlError::EmptyMolecule);
    }
    let coords = (present != 0).then(|| points.into_iter().map(Option::unwrap).collect());
    Ok((mol, coords))
}
