//! Shared lightweight XML attribute scanning, not a full XML validator.

use std::collections::HashMap;

/// Extract all `key="value"` pairs from a single XML tag line.
#[doc(hidden)]
///
/// Handles the most common cases:
/// - Double-quoted values
/// - Decimal/hexadecimal XML character references and basic XML entities: `&amp;` `&lt;` `&gt;` `&quot;` `&apos;`
/// - Ignores self-closing slash and angle brackets
pub(crate) fn parse_xml_attrs(line: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    let bytes = line.as_bytes();
    let len = bytes.len();
    let mut i = 0usize;

    while i < len {
        // Skip whitespace and punctuation between attributes.
        while i < len && (bytes[i].is_ascii_whitespace() || bytes[i] == b'/' || bytes[i] == b'>') {
            i += 1;
        }
        if i >= len {
            break;
        }

        // Read attribute key: stop at '=' or whitespace.
        let key_start = i;
        while i < len && bytes[i] != b'=' && !bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        let key = line.get(key_start..i).unwrap_or("").trim().to_string();
        if key.is_empty() {
            i += 1;
            continue;
        }

        // Expect '='
        while i < len && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= len || bytes[i] != b'=' {
            continue;
        }
        i += 1;

        // Expect opening '"' or '\''
        while i < len && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= len {
            break;
        }
        let quote = bytes[i];
        if quote != b'"' && quote != b'\'' {
            i += 1;
            continue;
        }
        i += 1;

        // Read until closing quote.
        let value_start = i;
        while i < len && bytes[i] != quote {
            i += 1;
        }
        let raw_value = line.get(value_start..i).unwrap_or("");
        let value = unescape_xml(raw_value);
        if i < len {
            i += 1;
        } // consume closing quote

        if !key.is_empty() {
            map.insert(key, value);
        }
    }
    map
}

// XML 1.0 Fifth Edition, sections 2.2 and 4.1. Unknown/invalid references
// remain literal for the legacy scanner; typed CML fields reject them.
fn reference(body: &str) -> Option<char> {
    match body {
        "amp" => return Some('&'),
        "lt" => return Some('<'),
        "gt" => return Some('>'),
        "quot" => return Some('"'),
        "apos" => return Some('\''),
        _ => {}
    }
    let (digits, radix) = if let Some(digits) = body.strip_prefix("#x") {
        (digits, 16)
    } else {
        (body.strip_prefix('#')?, 10)
    };
    if digits.is_empty()
        || !digits.bytes().all(|b| {
            if radix == 16 {
                b.is_ascii_hexdigit()
            } else {
                b.is_ascii_digit()
            }
        })
    {
        return None;
    }
    let value = u32::from_str_radix(digits, radix).ok()?;
    if !matches!(value, 0x9 | 0xa | 0xd | 0x20..=0xd7ff | 0xe000..=0xfffd | 0x10000..=0x10ffff) {
        return None;
    }
    char::from_u32(value)
}

fn unescape_xml(source: &str) -> String {
    if !source.contains('&') {
        return source.to_owned();
    }
    let mut output = String::with_capacity(source.len());
    // Split source once: decoded '&' characters are never scanned again.
    // Each segment is visited a constant number of times, even for long
    // runs of invalid ampersands followed by a single semicolon.
    for segment in source.split_inclusive(';') {
        if let Some(body) = segment.strip_suffix(';')
            && let Some(start) = body.rfind('&')
            && let Some(character) = reference(&body[start + 1..])
        {
            output.push_str(&body[..start]);
            output.push(character);
        } else {
            output.push_str(segment);
        }
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn references_are_decoded_once_with_xml_character_bounds() {
        assert_eq!(
            unescape_xml("&amp;lt; &amp;gt; &#38;#49;"),
            "&lt; &gt; &#49;"
        );
        assert_eq!(unescape_xml("&#49; &#x31; &#x3b1; &#128578;"), "1 1 α 🙂");
        assert_eq!(unescape_xml("&lt;&gt;&quot;&apos;&amp;"), "<>\"'&");
        for invalid in [
            "&#0;",
            "&#xD800;",
            "&#xFFFE;",
            "&#x110000;",
            "&#-1;",
            "&#+49;",
            "&#x+31;",
            "&#X31;",
            "&#;",
            "&#x;",
            "&unknown;",
            "&#999999999999;",
        ] {
            assert_eq!(unescape_xml(invalid), invalid);
        }
        for value in [9, 10, 13, 0x20, 0xd7ff, 0xe000, 0xfffd, 0x10000, 0x10ffff] {
            assert_eq!(
                unescape_xml(&format!("&#{value};")),
                char::from_u32(value).unwrap().to_string()
            );
        }
    }

    #[test]
    fn scanner_keeps_unicode_and_literal_reference_spelling() {
        let attrs = parse_xml_attrs("<molecule title='α &amp;lt; β &amp;gt; γ' id=\"a&#x31;\"/>");
        assert_eq!(attrs["title"], "α &lt; β &gt; γ");
        assert_eq!(attrs["id"], "a1");
        assert_eq!(unescape_xml(&"&".repeat(100_000)), "&".repeat(100_000));
        assert_eq!(
            unescape_xml(&("&".repeat(100_000) + "amp;")),
            "&".repeat(100_000)
        );
    }
}
