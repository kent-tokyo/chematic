use super::tokens::{Token, tokenize};
use super::{CcdComponent, CcdError, CcdParseLimits, Row, fields, invalid};
use std::collections::{BTreeMap, BTreeSet, HashSet};

struct Tables {
    row_bytes: usize,
    component_id: Option<String>,
    atoms: Vec<Row>,
    bonds: Vec<Row>,
    unhandled: BTreeSet<String>,
}
impl Tables {
    fn add_rows(
        &mut self,
        tags: &[String],
        values: &[Token<'_>],
        limits: &CcdParseLimits,
    ) -> Result<(), CcdError> {
        if tags.is_empty() || !values.len().is_multiple_of(tags.len()) {
            return Err(invalid("incomplete loop row or empty loop header"));
        }
        let category = tags[0].split('.').next().unwrap();
        if tags.iter().any(|s| s.split('.').next() != Some(category)) {
            return Err(invalid("mixed-category loop"));
        }
        if category != "_chem_comp_atom" && category != "_chem_comp_bond" {
            self.unhandled.extend(tags.iter().cloned());
            return Ok(());
        }
        let (rows, resource, limit) = if category == "_chem_comp_atom" {
            (&mut self.atoms, "atoms", limits.max_atoms)
        } else {
            (&mut self.bonds, "bonds", limits.max_bonds)
        };
        if values.len() / tags.len() > limit.saturating_sub(rows.len()) {
            return Err(CcdError::ResourceLimit { resource, limit });
        }
        Self::check_storage(&mut self.row_bytes, tags, values, limits)?;
        for chunk in values.chunks_exact(tags.len()) {
            rows.push(
                tags.iter()
                    .cloned()
                    .zip(chunk.iter().map(|t| t.value.into()))
                    .collect(),
            );
        }
        Ok(())
    }
    fn check_storage(
        row_bytes: &mut usize,
        tags: &[String],
        values: &[Token<'_>],
        limits: &CcdParseLimits,
    ) -> Result<(), CcdError> {
        let limit = limits.max_input_bytes.saturating_mul(4).min(32 << 20);
        let keys = tags
            .iter()
            .try_fold(0usize, |n, tag| n.checked_add(tag.len()));
        let added = keys
            .and_then(|n| n.checked_mul(values.len() / tags.len()))
            .and_then(|n| {
                values
                    .iter()
                    .try_fold(n, |n, value| n.checked_add(value.value.len()))
            });
        *row_bytes = added
            .and_then(|n| row_bytes.checked_add(n))
            .filter(|n| *n <= limit)
            .ok_or(CcdError::ResourceLimit {
                resource: "expanded loop bytes",
                limit,
            })?;
        Ok(())
    }
    fn scalar(
        &mut self,
        tag: String,
        value: &str,
        singletons: &mut BTreeMap<String, Row>,
    ) -> Result<(), CcdError> {
        let category = tag
            .split('.')
            .next()
            .ok_or_else(|| invalid("invalid tag"))?;
        if tag == "_chem_comp.id" {
            if self.component_id.replace(value.into()).is_some() {
                return Err(invalid("duplicate component ID"));
            }
        } else if category == "_chem_comp_atom" || category == "_chem_comp_bond" {
            if singletons
                .entry(category.into())
                .or_default()
                .insert(tag, value.into())
                .is_some()
            {
                return Err(invalid("duplicate singleton tag"));
            }
        } else {
            self.unhandled.insert(tag);
        }
        Ok(())
    }
}

fn loop_table(
    tokens: &[Token<'_>],
    i: &mut usize,
    tables: &mut Tables,
    limits: &CcdParseLimits,
) -> Result<(), CcdError> {
    *i += 1;
    let mut tags = Vec::new();
    let mut seen = HashSet::new();
    while *i < tokens.len() && tokens[*i].tag() {
        let tag = tokens[*i].value.to_ascii_lowercase();
        if !seen.insert(tag.clone()) {
            return Err(invalid("duplicate loop column"));
        }
        tags.push(tag);
        *i += 1;
    }
    let start = *i;
    while *i < tokens.len() && !tokens[*i].control() {
        *i += 1;
    }
    tables.add_rows(&tags, &tokens[start..*i], limits)
}

fn read_tables(tokens: &[Token<'_>], limits: &CcdParseLimits) -> Result<Tables, CcdError> {
    let mut tables = Tables {
        row_bytes: 0,
        component_id: None,
        atoms: Vec::new(),
        bonds: Vec::new(),
        unhandled: BTreeSet::new(),
    };
    let mut singletons = BTreeMap::new();
    let mut blocks = 0;
    let mut i = 0;
    while i < tokens.len() {
        let token = tokens[i];
        if token.bare
            && token
                .value
                .get(..5)
                .is_some_and(|s| s.eq_ignore_ascii_case("data_"))
        {
            blocks += 1;
            if blocks != 1 {
                return Err(invalid(
                    "multiple data blocks; pass one component at a time",
                ));
            }
            i += 1;
        } else if blocks == 0 {
            return Err(invalid("expected data block"));
        } else if token.keyword("loop_") {
            loop_table(tokens, &mut i, &mut tables, limits)?;
        } else if token.tag() {
            let value = tokens
                .get(i + 1)
                .filter(|t| !t.control())
                .ok_or_else(|| invalid("missing singleton value"))?;
            tables.scalar(
                token.value.to_ascii_lowercase(),
                value.value,
                &mut singletons,
            )?;
            i += 2;
        } else {
            return Err(invalid(format!("unexpected token {}", token.value)));
        }
    }
    for (category, row) in singletons {
        let (rows, resource, limit) = if category == "_chem_comp_atom" {
            (&mut tables.atoms, "atoms", limits.max_atoms)
        } else {
            (&mut tables.bonds, "bonds", limits.max_bonds)
        };
        if !rows.is_empty() {
            return Err(invalid("mixed singleton/loop representation"));
        }
        if limit == 0 {
            return Err(CcdError::ResourceLimit { resource, limit });
        }
        rows.push(row);
    }
    Ok(tables)
}

fn validate_links(component: &CcdComponent) -> Result<(), CcdError> {
    let mut ids = HashSet::new();
    for atom in &component.atoms {
        if !ids.insert(atom.atom_id.as_str()) {
            return Err(invalid("duplicate atom ID"));
        }
    }
    let mut bonds = HashSet::new();
    for bond in &component.bonds {
        let (a, b) = (bond.atom_id_1.as_str(), bond.atom_id_2.as_str());
        if a == b {
            return Err(invalid("self bond"));
        }
        if !ids.contains(a) || !ids.contains(b) {
            return Err(invalid("unknown bond endpoint"));
        }
        if !bonds.insert(if a < b { (a, b) } else { (b, a) }) {
            return Err(invalid("duplicate bond"));
        }
    }
    Ok(())
}

/// Parse one CCD component; rows and explicit H atoms retain source order.
/// Coordinates are optional complete triples; stereo labels remain annotations.
pub fn parse_ccd_component(text: &str) -> Result<CcdComponent, CcdError> {
    parse_ccd_component_with_limits(text, &CcdParseLimits::default())
}
pub fn parse_ccd_component_with_limits(
    text: &str,
    limits: &CcdParseLimits,
) -> Result<CcdComponent, CcdError> {
    if text.len() > limits.max_input_bytes {
        return Err(CcdError::ResourceLimit {
            resource: "input bytes",
            limit: limits.max_input_bytes,
        });
    }
    let tokens = tokenize(text, limits)?;
    let tables = read_tables(&tokens, limits)?;
    let id = tables
        .component_id
        .filter(|s| !s.is_empty() && s != "." && s != "?")
        .ok_or_else(|| invalid("missing chem_comp.id"))?;
    if tables.atoms.is_empty() {
        return Err(invalid("missing chem_comp_atom rows"));
    }
    let component = CcdComponent {
        component_id: id.clone(),
        atoms: tables
            .atoms
            .into_iter()
            .map(|r| fields::atom(r, &id))
            .collect::<Result<_, _>>()?,
        bonds: tables
            .bonds
            .into_iter()
            .map(|r| fields::bond(r, &id))
            .collect::<Result<_, _>>()?,
        unhandled_tags: tables.unhandled.into_iter().collect(),
    };
    validate_links(&component)?;
    Ok(component)
}
