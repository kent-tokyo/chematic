//! Logical tag lines for the lightweight graph reader, not an XML validator.

fn opaque_end(input: &str, start: usize) -> Option<usize> {
    let rest = &input[start..];
    let (prefix, closing) = if rest.starts_with("<!--") {
        (4, "-->")
    } else if rest.starts_with("<![CDATA[") {
        (9, "]]>")
    } else if rest.starts_with("<?") {
        (2, "?>")
    } else {
        return None;
    };
    // The strict structural gate rejects unterminated constructs; the lenient
    // reader may keep a graph preceding one, but must never parse its payload.
    Some(
        rest[prefix..]
            .find(closing)
            .map_or(input.len(), |end| start + prefix + end + closing.len()),
    )
}

pub(super) fn split_cml_tags(input: &str) -> String {
    let mut output = String::with_capacity(input.len() + input.len() / 8);
    let (mut in_tag, mut quote, mut previous_cr) = (false, None, false);
    let mut cursor = 0;
    while cursor < input.len() {
        if !in_tag && let Some(end) = opaque_end(input, cursor) {
            // Comments/CDATA/PIs are not graph markup, even when multiline.
            // Source-retaining document APIs operate on the original text.
            if !output.ends_with('\n') {
                output.push('\n');
            }
            cursor = end;
            continue;
        }
        let character = input[cursor..].chars().next().expect("UTF-8 boundary");
        if character == '<' && !in_tag {
            if !output.is_empty() && !output.ends_with('\n') {
                output.push('\n');
            }
            in_tag = true;
        }
        if in_tag && matches!(character, '\r' | '\n' | '\t') {
            // XML end-of-line plus CDATA attribute whitespace normalization.
            // Character references remain encoded here and are decoded later.
            if character != '\n' || !previous_cr {
                output.push(' ');
            }
            previous_cr = character == '\r';
        } else {
            output.push(character);
            previous_cr = false;
        }
        if in_tag {
            if let Some(active) = quote {
                if character == active {
                    quote = None;
                }
            } else if matches!(character, '\'' | '"') {
                quote = Some(character);
            } else if character == '>' {
                in_tag = false;
                output.push('\n');
            }
        }
        cursor += character.len_utf8();
    }
    output
}
