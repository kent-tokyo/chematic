from pathlib import Path
s=Path('/private/tmp/chematic-openbabel-lines-source-20261015');r=Path('/private/tmp/chematic-openbabel-lines-20261015');p=s/'crates/chematic-mol/src/mol3000.rs';t=p.read_text();a=t.index('    let logical_lines = collect_v30_lines(&lines);',t.index('pub fn read_mol_v3000_with_diagnostics'));b=t.index('    finish_v3000_report(parser, metadata, line2)',a)
t=t[:a]+'''    let mut parser = V3000Parser::new();
    consume_v3000_lines(&mut parser, &lines)?;
'''+t[b:]
t+='\ninclude!("mol3000/line_stream.rs");\n\n#[cfg(test)]\n#[path = "mol3000/line_stream_tests.rs"]\nmod line_stream_tests;\n';p.write_text(t)
(s/'crates/chematic-mol/src/mol3000/line_stream.rs').write_text('''/// Consume ordinary V30 payloads directly from validated physical lines.
///
/// Prefix lines before the first continuation need no owned text or logical
/// line vector. The first continuation and its entire suffix retain the
/// established folding path, keeping its exact consumption and error rules.
fn consume_v3000_lines(
    parser: &mut V3000Parser,
    lines: &[(usize, &str)],
) -> Result<(), MolParseError> {
    let mut tokens = Vec::new();
    for (index, &(line_num, raw)) in lines.iter().enumerate() {
        let Some(payload) = raw.strip_prefix(V30_PREFIX) else {
            continue;
        };
        if payload.ends_with('-') {
            // Owned folded payloads must outlive their borrowed tokens. Keep
            // the suffix together and reuse one token vector within it.
            drop(tokens);
            let logical_lines = collect_v30_lines(&lines[index..]);
            let mut suffix_tokens = Vec::new();
            for line in &logical_lines {
                suffix_tokens.clear();
                suffix_tokens.extend(line.payload.split_whitespace());
                parser.consume(line, &suffix_tokens)?;
            }
            return Ok(());
        }
        tokens.clear();
        tokens.extend(payload.split_whitespace());
        parser.consume(
            &LogicalLine {
                line_num,
                payload: std::borrow::Cow::Borrowed(payload),
            },
            &tokens,
        )?;
    }
    Ok(())
}
''')
old=(s/'crates/chematic-mol/src/mol3000/token_buffer_tests.rs').read_text();start=old.index('fn compare(');end=old.index('\nfn chain(',start)
reference='''//! Direct-line consumption versus the frozen logical-vector reader.
use super::*;

fn frozen_logical_vector(input: &str) -> Result<MolReadReport, MolParseError> {
    let lines = validate_v3000_physical_lines(input)?;
    let (metadata, line2) = v3000_metadata(&lines)?;
    let logical_lines = collect_v30_lines(&lines);
    let mut parser = V3000Parser::new();
    let mut tokens = Vec::new();
    for line in &logical_lines {
        tokens.clear();
        tokens.extend(line.payload.split_whitespace());
        parser.consume(line, &tokens)?;
    }
    finish_v3000_report(parser, metadata, line2)
}

'''
comparison=old[start:end].replace('frozen_per_line_tokens','frozen_logical_vector').replace('token storage changed','direct line consumption changed')
fixtures=old[old.index('\nfn chain('):]
(s/'crates/chematic-mol/src/mol3000/line_stream_tests.rs').write_text(reference+comparison+fixtures+'''
#[test]
fn continuation_at_every_logical_position_keeps_full_report() {
    let input = chain(8, " ", false);
    let physical = input.lines().collect::<Vec<_>>();
    for (index, &line) in physical.iter().enumerate() {
        let Some(payload) = line.strip_prefix(V30_PREFIX) else { continue; };
        for (split, _) in payload.char_indices().filter(|&(_, c)| c == ' ') {
            let mut lines = physical.iter().map(|s| s.to_string()).collect::<Vec<_>>();
            lines[index] = format!("M  V30 {}-\\nM  V30 {}", &payload[..split], &payload[split + 1..]);
            compare(&(lines.join("\\n") + "\\n"));
        }
    }
}

#[test]
fn malformed_continuations_and_physical_errors_keep_precedence() {
    let input = chain(3, " ", false);
    for suffix in ["-", "--", "---", " -", " -\\nM  V30 -", "-\\n\\n", "-\\nM  END"] {
        for needle in ["BEGIN CTAB", "COUNTS 3 2 0 0 0", "1 C 0 -0.0 0 0", "END ATOM", "END CTAB"] {
            compare(&input.replacen(needle, &format!("{needle}{suffix}"), 1));
        }
    }
    let invalid_atom = input.replace("1 C 0 -0.0 0 0", "1 C NaN -0.0 0 0");
    compare(&invalid_atom.replace("M  V30 END CTAB", "invalid physical line"));
    compare(&input.replace("M  END", "M  END\\nM  END"));
    for spacing in ["\\u{000b}", "\\u{000c}", "\\u{0085}", "\\u{00a0}", "\\u{2003}", "\\u{2028}", "\\u{2029}", "\\u{3000}"] {
        compare(&chain(3, spacing, false));
        compare(&chain(3, spacing, true));
    }
    compare(&input.replace('\\n', "\\r\\n"));
    compare(&input.trim_end_matches('\\n'));
}
''')
p=s/'CHANGELOG.md';t=p.read_text().replace('### Performance\n','### Performance\n\n- Consume ordinary V3000 payloads directly from validated physical lines,\n  avoiding the logical-line vector until a continuation requires the existing\n  folding path. Preserve diagnostics, Unicode tokens and error precedence.\n',1);p.write_text(t)
print('Prepared direct-line candidate and frozen-reader tests')

p=s/'crates/chematic-mol/src/mol3000/line_stream_tests.rs'
p.write_text(p.read_text().replace("compare(&input.trim_end_matches('\\n'));", "compare(input.trim_end_matches('\\n'));"))
