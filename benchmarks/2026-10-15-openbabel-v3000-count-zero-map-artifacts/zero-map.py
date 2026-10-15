from pathlib import Path
s=Path('/private/tmp/chematic-openbabel-count-zero-map-source-20261015')
p=s/'crates/chematic-mol/src/mol3000/writer.rs';t=p.read_text();old='        write!(out, " 0.0000 {}", atom.atom_map.unwrap_or(0))\n            .expect("writing to String cannot fail");';assert old in t;t=t.replace(old,'        if atom.atom_map.unwrap_or(0) == 0 {\n            out.push_str(" 0.0000 0");\n        } else {\n            write!(out, " 0.0000 {}", atom.atom_map.unwrap_or(0))\n                .expect("writing to String cannot fail");\n        }',1);p.write_text(t)
p=s/'CHANGELOG.md';p.write_text(p.read_text().replace('  state transitions, errors and diagnostics.','  state transitions, errors and diagnostics. Write absent or zero V3000 2D\n  atom maps directly, retaining the original nonzero-map formatting.'))
print('Zero atom-map direct suffix written')
