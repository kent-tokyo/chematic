//! ChemDraw abbreviation (`NodeType="Nickname"`/`"Fragment"`) expansion.
//!
//! ChemDraw stores a label such as `Boc`, `Ph` or `OTBS` as an `<n>` node
//! that holds a nested `<fragment>` with the expanded atoms. Bonds in the
//! enclosing fragment point at the outer node id; inside the nested fragment
//! one `ExternalConnectionPoint` node per attachment is bonded to the atom the
//! outer bond really reaches. With several connection points, the outer
//! node's `BondOrdering` lists the outer bond ids in connection-point order
//! (`ExternalConnectionNum`).

use std::collections::{HashMap, HashSet};

use chematic_core::AtomIdx;

use super::{CdxmlError, CdxmlParseLimits, FragAccum, ParseCtx};

const UNEXPANDABLE_ABBREVIATION_PREFIX: &str = "cannot expand CDXML abbreviation ";

/// A nested fragment read inside an abbreviation node, not yet spliced.
#[derive(Debug)]
pub(super) struct Expansion {
    node_id: String,
    label: Option<String>,
    inner: FragAccum,
}

/// One expanded ChemDraw abbreviation (for example `Boc`), so a caller can
/// show it collapsed again.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CdxmlAbbreviation {
    /// Index of the molecule (in the returned list) that holds the atoms.
    pub molecule_index: usize,
    /// The CDXML id of the abbreviation node that was replaced.
    pub node_id: String,
    /// The label text ChemDraw displays, when the node has one.
    pub label: Option<String>,
    /// The expanded atoms, in the returned molecule's numbering.
    pub atoms: Vec<AtomIdx>,
}

/// A non-fatal CDXML structure diagnostic: something the parse kept going
/// past instead of rejecting the document.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CdxmlParseDiagnostic {
    /// Index of the molecule (in the returned list) the object belonged to.
    pub molecule_index: usize,
    /// Stable machine-readable code, for example
    /// `generic_query_atom_as_wildcard`.
    pub code: &'static str,
    /// The CDXML id of the affected object, when it has one.
    pub object_id: Option<String>,
    /// Human-readable explanation.
    pub message: String,
}

/// Result of [`super::parse_cdxml_all_with_report`].
#[derive(Clone)]
pub struct CdxmlParseReport {
    /// The molecules and their 2D coordinates, as [`super::parse_cdxml_all`]
    /// returns them.
    pub molecules: Vec<(chematic_core::Molecule, Vec<(f64, f64)>)>,
    /// Expanded ChemDraw abbreviations.
    pub abbreviations: Vec<CdxmlAbbreviation>,
    /// Features kept in a reduced form or dropped instead of failing.
    pub diagnostics: Vec<CdxmlParseDiagnostic>,
}

impl CdxmlError {
    /// Returns the abbreviation label (or `node <id>` when it has none) when
    /// this error is the refusal to expand a ChemDraw abbreviation, for
    /// example because a connection point is missing.
    ///
    /// Like [`CdxmlError::unsupported_bond_order`], this keeps the public
    /// error enum source-compatible with earlier 1.x releases.
    pub fn unexpandable_abbreviation(&self) -> Option<&str> {
        match self {
            Self::InvalidDocument(message) => message
                .strip_prefix(UNEXPANDABLE_ABBREVIATION_PREFIX)
                .and_then(|rest| rest.split_once(": "))
                .map(|(label, _)| label.trim_matches('"')),
            _ => None,
        }
    }

    pub(super) fn is_unexpandable_abbreviation(message: &str) -> bool {
        message.starts_with(UNEXPANDABLE_ABBREVIATION_PREFIX)
    }
}

fn unexpandable(name: &str, reason: &str) -> CdxmlError {
    CdxmlError::InvalidDocument(format!(
        "{UNEXPANDABLE_ABBREVIATION_PREFIX}\"{name}\": {reason}"
    ))
}

impl FragAccum {
    /// Finish an `<n>` element: an abbreviation's label arrives in its `<t>`
    /// after the nested fragment.
    pub(super) fn close_node(&mut self) {
        let Some(node) = self.open_node.take() else {
            return;
        };
        let label = node.label.trim();
        if label.is_empty() {
            return;
        }
        let Some(id) = node.id else {
            return;
        };
        if let Some(expansion) = self.expansions.iter_mut().find(|e| e.node_id == id) {
            expansion.label = Some(label.to_string());
        }
    }

    /// Keep a nested fragment read inside the currently open node. Only
    /// abbreviation nodes are expanded; any other nested fragment is ignored.
    pub(super) fn attach_nested(
        &mut self,
        mut inner: FragAccum,
        limits: &CdxmlParseLimits,
        molecule_index: usize,
        ctx: &mut ParseCtx,
    ) -> Result<(), CdxmlError> {
        inner.expand_abbreviations(limits)?;
        let node = self.open_node.as_ref();
        match node.and_then(|node| node.id.clone().filter(|_| node.abbreviation)) {
            Some(node_id) => self.expansions.push(Expansion {
                node_id,
                label: None,
                inner,
            }),
            None => ctx.diagnose(
                molecule_index,
                "nested_fragment_ignored",
                node.and_then(|node| node.id.as_deref()),
                "nested fragment outside a Nickname/Fragment node was ignored".into(),
            ),
        }
        Ok(())
    }

    /// Replace every abbreviation node that has a nested fragment by the
    /// fragment's atoms, rewiring the outer bonds to the attachment atoms.
    pub(super) fn expand_abbreviations(
        &mut self,
        limits: &CdxmlParseLimits,
    ) -> Result<(), CdxmlError> {
        for expansion in std::mem::take(&mut self.expansions) {
            self.splice(expansion)?;
            if self.atoms.len() > limits.max_atoms {
                return Err(CdxmlError::ResourceLimit {
                    resource: "atoms",
                    actual: self.atoms.len(),
                    limit: limits.max_atoms,
                });
            }
        }
        Ok(())
    }

    fn splice(&mut self, expansion: Expansion) -> Result<(), CdxmlError> {
        let Expansion {
            node_id,
            label,
            mut inner,
        } = expansion;
        let Some(position) = self.atoms.iter().position(|atom| atom.id == node_id) else {
            return Ok(());
        };
        let name = label.clone().unwrap_or_else(|| format!("node {node_id}"));
        let attachments = inner.take_attachments(&name)?;
        let outer = self
            .bonds
            .iter()
            .enumerate()
            .filter(|(_, bond)| bond.begin == node_id || bond.end == node_id)
            .map(|(index, _)| index)
            .collect::<Vec<_>>();
        let ordering = self.bond_orderings.get(&node_id).cloned();
        for index in outer {
            let bond_id = self.bonds[index].id.clone();
            let target = attachment_for(
                &attachments,
                ordering.as_deref(),
                bond_id.as_deref(),
                &inner,
            )
            .ok_or_else(|| {
                unexpandable(
                    &name,
                    &format!(
                        "no external connection point for bond {}",
                        bond_id.as_deref().unwrap_or("without id")
                    ),
                )
            })?;
            let bond = &mut self.bonds[index];
            for endpoint in [&mut bond.begin, &mut bond.end] {
                if *endpoint == node_id {
                    *endpoint = target.clone();
                }
            }
        }
        self.atoms.remove(position);
        let atom_ids = inner.atoms.iter().map(|atom| atom.id.clone()).collect();
        self.abbreviations.push((node_id, label, atom_ids));
        self.abbreviations.append(&mut inner.abbreviations);
        self.atoms.append(&mut inner.atoms);
        self.bonds.append(&mut inner.bonds);
        Ok(())
    }

    /// Remove the connection points and their bonds, returning
    /// `(connection number, attached atom id)` per point.
    fn take_attachments(&mut self, name: &str) -> Result<Vec<(u32, String)>, CdxmlError> {
        let points = std::mem::take(&mut self.connection_points);
        let mut attachments = Vec::with_capacity(points.len());
        for (slot, (point, number)) in points.iter().enumerate() {
            let mut touching = self
                .bonds
                .iter()
                .filter(|bond| bond.begin == *point || bond.end == *point);
            let (Some(bond), None) = (touching.next(), touching.next()) else {
                return Err(unexpandable(
                    name,
                    &format!("external connection point {point} must have exactly one bond"),
                ));
            };
            let other = if bond.begin == *point {
                &bond.end
            } else {
                &bond.begin
            };
            if !self.atoms.iter().any(|atom| atom.id == *other) {
                return Err(unexpandable(
                    name,
                    &format!("external connection point {point} is bonded to unknown atom {other}"),
                ));
            }
            attachments.push((number.unwrap_or(slot as u32 + 1), other.clone()));
        }
        let ids = points
            .iter()
            .map(|(id, _)| id.as_str())
            .collect::<HashSet<_>>();
        self.bonds
            .retain(|bond| !ids.contains(bond.begin.as_str()) && !ids.contains(bond.end.as_str()));
        Ok(attachments)
    }

    /// Connection points outside an abbreviation mark nothing chemical:
    /// drop them and the bonds to them.
    pub(super) fn drop_top_level_connection_points(
        &mut self,
        molecule_index: usize,
        ctx: &mut ParseCtx,
    ) {
        for (id, _) in std::mem::take(&mut self.connection_points) {
            self.bonds.retain(|bond| bond.begin != id && bond.end != id);
            ctx.diagnose(
                molecule_index,
                "connection_point_dropped",
                Some(&id),
                "external connection point outside an abbreviation was dropped".into(),
            );
        }
    }

    pub(super) fn record_abbreviations(&mut self, molecule_index: usize, ctx: &mut ParseCtx) {
        let positions = self
            .atoms
            .iter()
            .enumerate()
            .map(|(index, atom)| (atom.id.as_str(), AtomIdx(index as u32)))
            .collect::<HashMap<_, _>>();
        for (node_id, label, ids) in std::mem::take(&mut self.abbreviations) {
            let atoms = ids
                .iter()
                .filter_map(|id| positions.get(id.as_str()).copied())
                .collect();
            ctx.abbreviations.push(super::CdxmlAbbreviation {
                molecule_index,
                node_id,
                label,
                atoms,
            });
        }
    }
}

/// The inner atom an outer bond attaches to.
fn attachment_for(
    attachments: &[(u32, String)],
    ordering: Option<&[String]>,
    bond_id: Option<&str>,
    inner: &FragAccum,
) -> Option<String> {
    let slot = ordering
        .zip(bond_id)
        .and_then(|(ordering, id)| ordering.iter().position(|entry| entry == id));
    if let Some(slot) = slot {
        let number = slot as u32 + 1;
        return attachments
            .iter()
            .find(|(n, _)| *n == number)
            .map(|(_, atom)| atom.clone());
    }
    match (attachments, inner.atoms.as_slice()) {
        ([(_, atom)], _) => Some(atom.clone()),
        ([], [atom]) => Some(atom.id.clone()),
        _ => None,
    }
}
