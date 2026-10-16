use super::{ChemDrawError as Error, ChemDrawLimits, ChemDrawOmission, ChemDrawOptions};
use std::collections::BTreeSet;

pub(super) struct Reader<'a> {
    pub input: &'a [u8],
    pub offset: usize,
    pub options: &'a ChemDrawOptions,
    pub limits: &'a ChemDrawLimits,
    pub omissions: Vec<ChemDrawOmission>,
    items: usize,
}
impl Reader<'_> {
    pub fn invalid(&self, detail: &'static str) -> Error {
        Error::Invalid {
            offset: self.offset,
            detail,
        }
    }
    pub fn take(&mut self, n: usize) -> Result<&[u8], Error> {
        let end = self
            .offset
            .checked_add(n)
            .ok_or_else(|| self.invalid("length overflow"))?;
        let b = self
            .input
            .get(self.offset..end)
            .ok_or_else(|| self.invalid("truncated item"))?;
        self.offset = end;
        Ok(b)
    }
    pub fn u16(&mut self) -> Result<u16, Error> {
        Ok(u16::from_le_bytes(self.take(2)?.try_into().unwrap()))
    }
    fn u32(&mut self) -> Result<u32, Error> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn item(&mut self) -> Result<(), Error> {
        if self.items == self.limits.max_items {
            return Err(Error::Limit {
                resource: "items",
                limit: self.limits.max_items,
            });
        }
        self.items += 1;
        Ok(())
    }
    pub fn omit(&mut self, tag: u16, offset: usize) -> Result<(), Error> {
        if !self.options.allow_metadata_loss {
            return Err(Error::MetadataLoss { offset, tag });
        }
        self.omissions.push(ChemDrawOmission { offset, tag });
        Ok(())
    }
    fn append(&self, target: &mut String, source: &str) -> Result<(), Error> {
        if source.len() > self.limits.max_output_bytes.saturating_sub(target.len()) {
            return Err(Error::Limit {
                resource: "projected XML bytes",
                limit: self.limits.max_output_bytes,
            });
        }
        target.push_str(source);
        Ok(())
    }
    fn property(&mut self, owner: u16, tag: u16, offset: usize) -> Result<String, Error> {
        let n = self.u16()?;
        let n = if n == u16::MAX {
            self.u32()? as usize
        } else {
            usize::from(n)
        };
        if n > self.limits.max_property_bytes {
            return Err(Error::Limit {
                resource: "property bytes",
                limit: self.limits.max_property_bytes,
            });
        }
        let start = self.offset;
        self.take(n)?;
        super::properties::attribute(self, owner, tag, offset, start, n)
    }
    fn object(&mut self, tag: u16, id: u32, depth: usize, parent: u16) -> Result<String, Error> {
        if depth > self.limits.max_depth {
            return Err(Error::Limit {
                resource: "object depth",
                limit: self.limits.max_depth,
            });
        }
        let name = object_name(tag, parent, self.offset - 6)?;
        if tag == 0x8006 {
            self.omit(tag, self.offset - 6)?;
        }
        let mut attrs = format!(" id=\"{id}\"");
        let mut children = String::new();
        let mut properties = BTreeSet::new();
        let mut ids = BTreeSet::new();
        loop {
            let at = self.offset;
            let child = self.u16()?;
            if child == 0 {
                break;
            }
            self.item()?;
            if child & 0x8000 != 0 {
                let child_id = self.u32()?;
                if child_id != 0 && !ids.insert(child_id) {
                    return Err(self.invalid("duplicate child object ID"));
                }
                let xml = self.object(child, child_id, depth + 1, tag)?;
                self.append(&mut children, &xml)?;
            } else {
                if !properties.insert(child) {
                    return Err(self.invalid("duplicate property"));
                }
                let attr = self.property(tag, child, at)?;
                self.append(&mut attrs, &attr)?;
            }
        }
        required(tag, id, &properties, &attrs, self.offset)?;
        if name.is_empty() {
            return Ok(String::new());
        }
        let mut xml = String::new();
        self.append(&mut xml, &format!("<{name}{attrs}>\n"))?;
        self.append(&mut xml, &children)?;
        self.append(&mut xml, &format!("</{name}>\n"))?;
        Ok(xml)
    }
}
pub(super) fn decode(
    input: &[u8],
    options: &ChemDrawOptions,
    limits: &ChemDrawLimits,
) -> Result<(String, Vec<ChemDrawOmission>), Error> {
    // Both the 22-byte header in RDKit's ChemDraw fixture and the published
    // 28-byte header layout. Neither accepts an arbitrary reserved prefix.
    let header = if input.get(22..24) == Some(&[0, 0x80]) {
        22
    } else {
        28
    };
    if !input.starts_with(b"VjCD0100\x04\x03\x02\x01")
        || input
            .get(12..header)
            .is_none_or(|b| b.iter().any(|&v| v != 0))
    {
        return Err(Error::Invalid {
            offset: 0,
            detail: "CDX header/signature",
        });
    }
    let mut r = Reader {
        input,
        offset: header,
        options,
        limits,
        omissions: Vec::new(),
        items: 0,
    };
    if r.u16()? != 0x8000 {
        return Err(r.invalid("expected document object"));
    }
    r.item()?;
    let id = r.u32()?;
    let xml = r.object(0x8000, id, 1, 0)?;
    if r.u16()? != 0 || r.offset != input.len() {
        return Err(r.invalid("missing file terminator or trailing data"));
    }
    Ok((xml, r.omissions))
}

fn object_name(tag: u16, parent: u16, offset: usize) -> Result<&'static str, Error> {
    Ok(match (tag, parent) {
        (0x8000, 0) => "CDXML",
        (0x8001, 0x8000) => "page",
        (0x8002, 0x8001 | 0x8002) => "group",
        (0x8003, 0x8001 | 0x8002) => "fragment",
        (0x8004, 0x8003) => "n",
        (0x8005, 0x8003) => "b",
        (0x8006, 0x8001 | 0x8002 | 0x8004) => "",
        _ => return Err(Error::Unsupported { offset, tag }),
    })
}

fn required(
    tag: u16,
    id: u32,
    properties: &BTreeSet<u16>,
    attrs: &str,
    offset: usize,
) -> Result<(), Error> {
    if tag == 0x8005 {
        let attrs = crate::cml::parse_xml_attrs(attrs);
        if attrs.get("B").is_some_and(|b| Some(b) == attrs.get("E")) {
            return Err(Error::Invalid {
                offset,
                detail: "self bond",
            });
        }
    }
    let detail = if matches!(tag, 0x8004 | 0x8005) && id == 0 {
        Some("atom/bond ID must be positive")
    } else if tag == 0x8004 && !properties.contains(&0x0200) {
        Some("atom requires 2D position")
    } else if tag == 0x8005 && !(properties.contains(&0x0604) && properties.contains(&0x0605)) {
        Some("bond requires B/E endpoints")
    } else {
        None
    };
    if let Some(detail) = detail {
        Err(Error::Invalid { offset, detail })
    } else {
        Ok(())
    }
}
