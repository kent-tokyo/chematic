use super::*;
use fields::layout;
pub fn parse_gromos96(text: &str) -> Result<Gromos96Document, Gromos96Error> {
    parse_gromos96_with_limits(text, Gromos96LengthUnit::Unspecified, &Default::default())
}
pub fn parse_gromos96_with_limits(
    text: &str,
    length_unit: Gromos96LengthUnit,
    limits: &Gromos96ParseLimits,
) -> Result<Gromos96Document, Gromos96Error> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "input byte limit exceeded"));
    }
    let mut document = Gromos96Document {
        length_unit,
        blocks: Vec::new(),
        trailing: Vec::new(),
    };
    let mut active: Option<Gromos96Block> = None;
    let mut leading = Vec::new();
    for (i, line) in text.split_terminator('\n').enumerate() {
        let n = i + 1;
        let line = line.strip_suffix('\r').unwrap_or(line);
        check_line(line, n, limits)?;
        if let Some(block) = &mut active {
            if line.trim() == "END" {
                document.blocks.push(active.take().expect("active block"));
            } else {
                block.lines.push(line.into());
            }
        } else if layout(line) {
            leading.push(line.into());
        } else {
            let kind = line.trim();
            if !matches!(
                kind,
                "TITLE"
                    | "TIMESTEP"
                    | "POSITION"
                    | "POSITIONRED"
                    | "VELOCITY"
                    | "VELOCITYRED"
                    | "BOX"
            ) {
                return Err(error(n, "unsupported block or trailing data"));
            }
            if document.blocks.len() >= limits.max_frames.saturating_mul(4).saturating_add(1) {
                return Err(error(n, "block limit exceeded"));
            }
            active = Some(Gromos96Block {
                kind: kind.into(),
                leading: core::mem::take(&mut leading),
                lines: Vec::new(),
            });
        }
    }
    if active.is_some() {
        return Err(error(0, "truncated block: missing END"));
    }
    document.trailing = leading;
    document.validate_with_limits(limits)?;
    Ok(document)
}
pub fn write_gromos96(document: &Gromos96Document) -> Result<String, Gromos96Error> {
    write_gromos96_with_limits(document, &Default::default())
}
pub fn write_gromos96_with_limits(
    document: &Gromos96Document,
    limits: &Gromos96ParseLimits,
) -> Result<String, Gromos96Error> {
    document.validate_with_limits(limits)?;
    let mut output = String::new();
    let mut line = |text: &str| -> Result<(), Gromos96Error> {
        if output.len().saturating_add(text.len()).saturating_add(1) > limits.max_input_bytes {
            return Err(error(0, "text output byte limit exceeded"));
        }
        output.push_str(text);
        output.push('\n');
        Ok(())
    };
    for block in &document.blocks {
        for value in &block.leading {
            line(value)?;
        }
        line(&block.kind)?;
        for value in &block.lines {
            line(value)?;
        }
        line("END")?;
    }
    for value in &document.trailing {
        line(value)?;
    }
    Ok(output)
}
