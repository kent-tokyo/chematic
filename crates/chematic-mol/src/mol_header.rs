//! Shared CTfile title serialization for V2000/V3000 and their SDF wrappers.

/// Append one title line, bounded to 80 UTF-8 bytes without splitting a scalar.
/// Controls and Unicode line separators occupy one ASCII space each. The first
/// scalar that would exceed the budget ends the prefix (no later suffix scan).
/// Caller-owned metadata is never changed. This deliberately conservative byte
/// budget also fits 80-character consumers; UTF-8 interoperability varies.
pub(crate) fn append_title(out: &mut String, name: &str) {
    let start = out.len();
    for ch in name.chars() {
        let ch = if ch.is_control() || matches!(ch, '\u{2028}' | '\u{2029}') {
            ' '
        } else {
            ch
        };
        if out.len() - start + ch.len_utf8() > 80 {
            break;
        }
        out.push(ch);
    }
    out.push('\n');
}
