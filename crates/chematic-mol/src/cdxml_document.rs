//! Loss-preserving CDXML document envelope.
//!
//! The molecule-only parser intentionally ignores presentation objects. This
//! API keeps the original XML as the source of truth while exposing page and
//! object boundaries for editors that need to inspect or edit document state.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::cdxml::{CdxmlError, CdxmlParseLimits};
use crate::cml::parse_xml_attrs;

pub type CdxmlValue = Value;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CdxmlObject {
    pub tag: String,
    pub attributes: BTreeMap<String, CdxmlValue>,
    pub raw_xml: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CdxmlPage {
    pub id: Option<String>,
    pub attributes: BTreeMap<String, CdxmlValue>,
    pub children: Vec<CdxmlObject>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CdxmlDocument {
    pub document_attributes: BTreeMap<String, CdxmlValue>,
    pub pages: Vec<CdxmlPage>,
    raw_xml: String,
}

/// Safe, source-oriented edits for document-level CDXML consumers.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", tag = "kind")]
pub enum CdxmlEdit {
    SetPageAttribute {
        page_id: String,
        key: String,
        value: String,
    },
    ReplaceObject {
        page_id: String,
        object_index: usize,
        raw_xml: String,
    },
    SetObjectAttribute {
        page_id: String,
        object_index: usize,
        key: String,
        value: String,
    },
    InsertObject {
        page_id: String,
        object_index: usize,
        raw_xml: String,
    },
    RemoveObject {
        page_id: String,
        object_index: usize,
    },
    /// Replace an object at a parent-to-child sibling path. `[0, 1]` means
    /// the second child of the first page-level object.
    ReplaceObjectPath {
        page_id: String,
        path: Vec<usize>,
        raw_xml: String,
    },
}

impl CdxmlDocument {
    /// Parse a CDXML document without discarding unknown presentation data.
    pub fn parse(input: &str) -> Result<Self, CdxmlError> {
        Self::parse_with_limits(input, &CdxmlParseLimits::default())
    }

    pub fn parse_with_limits(input: &str, limits: &CdxmlParseLimits) -> Result<Self, CdxmlError> {
        if input.len() > limits.max_input_bytes {
            return Err(CdxmlError::ResourceLimit {
                resource: "input bytes",
                actual: input.len(),
                limit: limits.max_input_bytes,
            });
        }
        let mut document_attributes = BTreeMap::new();
        let mut pages = Vec::new();
        let mut current: Option<CdxmlPage> = None;
        for (line_no, raw) in input.lines().enumerate() {
            if line_no >= limits.max_lines {
                return Err(CdxmlError::ResourceLimit {
                    resource: "lines",
                    actual: line_no + 1,
                    limit: limits.max_lines,
                });
            }
            if raw.len() > limits.max_line_bytes {
                return Err(CdxmlError::ResourceLimit {
                    resource: "line bytes",
                    actual: raw.len(),
                    limit: limits.max_line_bytes,
                });
            }
            let line = raw.trim();
            if line.starts_with("<page") && !line.starts_with("</page") {
                if pages.len() >= limits.max_fragments {
                    return Err(CdxmlError::ResourceLimit {
                        resource: "pages",
                        actual: pages.len() + 1,
                        limit: limits.max_fragments,
                    });
                }
                let attrs = parse_xml_attrs(line);
                let id = attrs.get("id").cloned();
                current = Some(CdxmlPage {
                    id,
                    attributes: attrs
                        .into_iter()
                        .map(|(k, v)| (k, Value::String(v)))
                        .collect(),
                    children: Vec::new(),
                });
            } else if line.starts_with("</page") {
                if let Some(page) = current.take() {
                    pages.push(page);
                }
            } else if let Some(page) = current.as_mut() {
                if line.starts_with('<')
                    && !line.starts_with("</")
                    && !line.starts_with("<?")
                    && !line.starts_with("<!")
                {
                    let tag = line
                        .trim_start_matches('<')
                        .split(|c: char| c.is_whitespace() || c == '>' || c == '/')
                        .next()
                        .unwrap_or_default()
                        .to_string();
                    if page.children.len() >= limits.max_bonds.saturating_add(limits.max_atoms) {
                        return Err(CdxmlError::ResourceLimit {
                            resource: "objects",
                            actual: page.children.len() + 1,
                            limit: limits.max_bonds.saturating_add(limits.max_atoms),
                        });
                    }
                    let attrs = parse_xml_attrs(line)
                        .into_iter()
                        .map(|(k, v)| (k, Value::String(v)))
                        .collect();
                    page.children.push(CdxmlObject {
                        tag,
                        attributes: attrs,
                        raw_xml: raw.to_string(),
                    });
                }
            } else if line.starts_with("<CDXML") {
                document_attributes = parse_xml_attrs(line)
                    .into_iter()
                    .map(|(k, v)| (k, Value::String(v)))
                    .collect();
            }
        }
        if current.is_some() {
            return Err(CdxmlError::InvalidCoords("unterminated page".into()));
        }
        Ok(Self {
            document_attributes,
            pages,
            raw_xml: input.to_string(),
        })
    }

    /// Return the exact source representation, preserving unknown tags and attributes.
    pub fn write(&self) -> String {
        self.raw_xml.clone()
    }

    /// Number of document pages, including pages containing only presentation
    /// objects or otherwise unknown CDXML content.
    pub fn page_count(&self) -> usize {
        self.pages.len()
    }

    /// Return page IDs in source order. Missing IDs remain `None` so callers
    /// cannot accidentally address a page using an invented identifier.
    pub fn page_ids(&self) -> Vec<Option<&str>> {
        self.pages.iter().map(|page| page.id.as_deref()).collect()
    }

    /// Apply a JSON-encoded document edit and return the reparsed document.
    ///
    /// This is the binding-neutral command boundary. The edit is applied to
    /// the original XML, then parsed again under the normal resource limits.
    pub fn apply_json_edit(&self, edit_json: &str) -> Result<Self, CdxmlError> {
        let edit: CdxmlEdit = serde_json::from_str(edit_json).map_err(|error| {
            CdxmlError::InvalidCoords(format!("invalid CDXML edit JSON: {error}"))
        })?;
        self.apply(&edit)
    }

    /// Apply a bounded edit and reparse, so indexes/attributes stay consistent.
    pub fn apply(&self, edit: &CdxmlEdit) -> Result<Self, CdxmlError> {
        let mut lines: Vec<String> = self.raw_xml.lines().map(str::to_owned).collect();
        let page_id = page_id_for(edit);
        let page = self
            .pages
            .iter()
            .position(|p| p.id.as_deref() == Some(page_id))
            .ok_or_else(|| CdxmlError::UnknownAtomRef("unknown page id".into()))?;
        if let CdxmlEdit::ReplaceObjectPath { path, raw_xml, .. } = edit {
            if path.is_empty() {
                return Err(CdxmlError::UnknownAtomRef("object path is empty".into()));
            }
            let mut current_page = 0usize;
            let mut in_target = false;
            let mut open_paths: Vec<Vec<usize>> = Vec::new();
            let mut sibling_counts = vec![0usize];
            for i in 0..lines.len() {
                let trimmed = lines[i].trim().to_owned();
                if trimmed.starts_with("<page") && !trimmed.starts_with("</page") {
                    in_target = current_page == page;
                    open_paths.clear();
                    sibling_counts.clear();
                    sibling_counts.push(0);
                    continue;
                }
                if !in_target {
                    if trimmed.starts_with("</page") {
                        current_page += 1;
                    }
                    continue;
                }
                if trimmed.starts_with("</page") {
                    in_target = false;
                    current_page += 1;
                    continue;
                }
                if trimmed.starts_with("</") {
                    if !open_paths.is_empty() {
                        open_paths.pop();
                        sibling_counts.pop();
                    }
                    continue;
                }
                if !is_editable_object_line(&trimmed) {
                    continue;
                }
                let mut object_path = open_paths.last().cloned().unwrap_or_default();
                object_path.push(*sibling_counts.last().unwrap_or(&0));
                *sibling_counts.last_mut().unwrap_or(&mut 0) += 1;
                if object_path == *path {
                    lines[i] = raw_xml.clone();
                    return Self::parse(&format!("{}\n", lines.join("\n")));
                }
                if !trimmed.ends_with("/>") {
                    open_paths.push(object_path);
                    sibling_counts.push(0);
                }
            }
            return Err(CdxmlError::UnknownAtomRef("object path not found".into()));
        }
        let mut page_index = 0usize;
        let mut in_page = false;
        let mut object_index = 0usize;
        if matches!(
            edit,
            CdxmlEdit::InsertObject { .. } | CdxmlEdit::RemoveObject { .. }
        ) {
            let mut current_page = 0usize;
            let mut in_target = false;
            let mut seen_objects = 0usize;
            for i in 0..lines.len() {
                let trimmed = lines[i].trim().to_owned();
                if trimmed.starts_with("<page") && !trimmed.starts_with("</page") {
                    in_target = current_page == page;
                    seen_objects = 0;
                } else if in_target
                    && trimmed.starts_with("<")
                    && !trimmed.starts_with("</")
                    && !trimmed.starts_with("<?")
                    && !trimmed.starts_with("<!")
                {
                    if let CdxmlEdit::RemoveObject {
                        object_index: target,
                        ..
                    } = edit
                        && seen_objects == *target
                    {
                        lines.remove(i);
                        return Self::parse(&format!("{}\n", lines.join("\n")));
                    }
                    seen_objects += 1;
                } else if in_target && trimmed.starts_with("</page") {
                    if let CdxmlEdit::InsertObject {
                        object_index: target,
                        raw_xml,
                        ..
                    } = edit
                        && seen_objects == *target
                    {
                        lines.insert(i, raw_xml.clone());
                        return Self::parse(&format!("{}\n", lines.join("\n")));
                    }
                    current_page += 1;
                    in_target = false;
                }
            }
            return Err(CdxmlError::UnknownAtomRef(
                "object index out of range".into(),
            ));
        }
        for line in &mut lines {
            let trimmed = line.trim().to_owned();
            if trimmed.starts_with("<page") && !trimmed.starts_with("</page") {
                in_page = page_index == page;
                object_index = 0;
                if in_page && let CdxmlEdit::SetPageAttribute { key, value, .. } = edit {
                    let mut attrs = parse_xml_attrs(&trimmed);
                    attrs.insert(key.clone(), value.clone());
                    let mut rebuilt = String::from("<page");
                    for (k, v) in attrs {
                        rebuilt.push_str(&format!(" {k}=\"{}\"", xml_escape(&v)));
                    }
                    rebuilt.push('>');
                    *line = rebuilt;
                }
            }
            if in_page
                && trimmed.starts_with('<')
                && !trimmed.starts_with("</")
                && !trimmed.starts_with("<page")
                && !trimmed.starts_with("<?")
                && !trimmed.starts_with("<!")
            {
                match edit {
                    CdxmlEdit::ReplaceObject {
                        object_index: target,
                        raw_xml,
                        ..
                    } if object_index == *target => *line = raw_xml.clone(),
                    CdxmlEdit::SetObjectAttribute {
                        object_index: target,
                        key,
                        value,
                        ..
                    } if object_index == *target => {
                        let mut attrs = parse_xml_attrs(&trimmed);
                        attrs.insert(key.clone(), value.clone());
                        let tag = trimmed
                            .trim_start_matches('<')
                            .split(|c: char| c.is_whitespace() || c == '>' || c == '/')
                            .next()
                            .unwrap_or_default();
                        let self_closing = trimmed.ends_with("/>");
                        let mut rebuilt = format!("<{tag}");
                        for (k, v) in attrs {
                            rebuilt.push_str(&format!(" {k}=\"{}\"", xml_escape(&v)));
                        }
                        rebuilt.push_str(if self_closing { "/>" } else { ">" });
                        *line = rebuilt;
                    }
                    _ => {}
                }
                object_index += 1;
            }
            if trimmed.starts_with("</page") {
                in_page = false;
                page_index += 1;
            }
        }
        Self::parse(&format!("{}\n", lines.join("\n")))
    }

    /// A JSON-safe structural summary for editor and binding layers.
    pub fn to_json(&self) -> Value {
        let pages = self
            .pages
            .iter()
            .map(|p| {
                serde_json::json!({
                    "id": p.id,
                    "attributes": p.attributes,
                    "children": p.children.iter().map(|o| serde_json::json!({
                        "tag": o.tag, "attributes": o.attributes, "raw_xml": o.raw_xml
                    })).collect::<Vec<_>>()
                })
            })
            .collect::<Vec<_>>();
        serde_json::json!({ "schema": "chematic.cdxml-document.v1", "document_attributes": self.document_attributes, "pages": pages })
    }
}

fn xml_escape(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('"', "&quot;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

fn page_id_for(edit: &CdxmlEdit) -> &str {
    match edit {
        CdxmlEdit::SetPageAttribute { page_id, .. }
        | CdxmlEdit::ReplaceObject { page_id, .. }
        | CdxmlEdit::SetObjectAttribute { page_id, .. }
        | CdxmlEdit::InsertObject { page_id, .. }
        | CdxmlEdit::RemoveObject { page_id, .. }
        | CdxmlEdit::ReplaceObjectPath { page_id, .. } => page_id,
    }
}

fn is_editable_object_line(line: &str) -> bool {
    line.starts_with('<')
        && !line.starts_with("</")
        && !line.starts_with("<?")
        && !line.starts_with("<!")
        && !line.starts_with("<page")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn preserves_pages_objects_and_unknown_attributes() {
        let input = "<CDXML customRoot=\"keep\">\n<page id=\"p2\" custom=\"x\">\n<n id=\"1\" Element=\"6\" customNode=\"y\"/>\n<arrow id=\"a1\" Head3=\"yes\"/>\n</page>\n</CDXML>\n";
        let doc = CdxmlDocument::parse(input).unwrap();
        assert_eq!(doc.pages.len(), 1);
        assert_eq!(doc.pages[0].id.as_deref(), Some("p2"));
        assert_eq!(doc.pages[0].children[1].tag, "arrow");
        assert_eq!(
            doc.pages[0].children[1].attributes["Head3"],
            Value::String("yes".into())
        );
        assert_eq!(doc.write(), input);
    }

    #[test]
    fn rejects_page_budget() {
        let input = "<CDXML>\n<page id=\"p1\"></page>\n</CDXML>";
        let limits = CdxmlParseLimits {
            max_fragments: 0,
            ..Default::default()
        };
        assert!(matches!(
            CdxmlDocument::parse_with_limits(input, &limits),
            Err(CdxmlError::ResourceLimit {
                resource: "pages",
                ..
            })
        ));
    }

    #[test]
    fn applies_page_and_object_edits_without_losing_unknown_data() {
        let input = "<CDXML>\n<page id=\"p1\" keep=\"yes\">\n<arrow id=\"a1\" custom=\"z\"/>\n</page>\n</CDXML>";
        let doc = CdxmlDocument::parse(input).unwrap();
        let doc = doc
            .apply(&CdxmlEdit::SetPageAttribute {
                page_id: "p1".into(),
                key: "title".into(),
                value: "Page 1".into(),
            })
            .unwrap();
        let doc = doc
            .apply(&CdxmlEdit::ReplaceObject {
                page_id: "p1".into(),
                object_index: 0,
                raw_xml: "<text id=\"t1\" custom=\"z\"/>".into(),
            })
            .unwrap();
        let doc = doc
            .apply(&CdxmlEdit::SetObjectAttribute {
                page_id: "p1".into(),
                object_index: 0,
                key: "label".into(),
                value: "A&B".into(),
            })
            .unwrap();
        assert!(doc.write().contains("title=\"Page 1\""));
        assert!(doc.write().contains("<text"));
        assert!(doc.write().contains("custom=\"z\""));
        assert!(doc.write().contains("label=\"A&amp;B\""));

        let doc = doc
            .apply(&CdxmlEdit::InsertObject {
                page_id: "p1".into(),
                object_index: 1,
                raw_xml: "<graphic id=\"g1\"/>".into(),
            })
            .unwrap();
        assert!(doc.write().contains("<graphic id=\"g1\"/>"));
        let doc = doc
            .apply(&CdxmlEdit::RemoveObject {
                page_id: "p1".into(),
                object_index: 1,
            })
            .unwrap();
        assert!(!doc.write().contains("<graphic id=\"g1\"/>"));
    }

    #[test]
    fn replaces_nested_object_by_loss_preserving_path() {
        let input = "<CDXML>\n<page id=\"p1\">\n<group id=\"g1\" unknown=\"keep\">\n<arrow id=\"a1\" Custom=\"keep\"/>\n</group>\n</page>\n</CDXML>";
        let doc = CdxmlDocument::parse(input).unwrap();
        let doc = doc
            .apply(&CdxmlEdit::ReplaceObjectPath {
                page_id: "p1".into(),
                path: vec![0, 0],
                raw_xml: "<text id=\"t1\" Custom=\"keep\"/>".into(),
            })
            .unwrap();
        let output = doc.write();
        assert!(output.contains("<group id=\"g1\" unknown=\"keep\">"));
        assert!(output.contains("<text id=\"t1\" Custom=\"keep\"/>"));
        assert!(!output.contains("<arrow id=\"a1\""));
    }

    #[test]
    fn applies_json_command_across_multiple_pages() {
        let input = "<CDXML>\n<page id=\"p1\">\n<arrow id=\"a1\"/>\n</page>\n<page id=\"p2\">\n<text id=\"t1\"/>\n</page>\n</CDXML>";
        let doc = CdxmlDocument::parse(input).unwrap();
        assert_eq!(doc.page_count(), 2);
        assert_eq!(doc.page_ids(), vec![Some("p1"), Some("p2")]);
        let edited = doc
            .apply_json_edit(
                r#"{"kind":"set_page_attribute","page_id":"p2","key":"title","value":"Page 2"}"#,
            )
            .unwrap();
        assert!(
            edited.write().contains("title=\"Page 2\""),
            "{}",
            edited.write()
        );
        assert!(edited.write().contains("<arrow id=\"a1\"/>"));
    }
}
