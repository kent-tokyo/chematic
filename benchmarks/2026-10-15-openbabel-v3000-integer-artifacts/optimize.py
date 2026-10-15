from pathlib import Path
s=Path('/private/tmp/chematic-openbabel-integer-source-20261015')
p=s/'crates/chematic-mol/src/mol2000/writer.rs';t=p.read_text();t=t.replace('pub(super) fn push_right_aligned_u32','pub(crate) fn push_right_aligned_u32',1).replace('if width == 3 && value < 1000 {','if matches!(width, 0 | 3) && value < 1000 {',1);needle='        out.push_str(TEXT[value as usize]);\n        return;';replacement='''        let field = TEXT[value as usize];
        if width == 0 {
            let start = if value < 10 { 2 } else if value < 100 { 1 } else { 0 };
            out.push_str(&field[start..]);
        } else {
            out.push_str(field);
        }
        return;''';assert needle in t;t=t.replace(needle,replacement,1).replace('// V2000 counts and bond fields are almost always three columns wide.','// V2000 fields use three columns; V3000 ids use the same digits unpadded.');p.write_text(t)
p=s/'crates/chematic-mol/src/mol2000.rs';t=p.read_text();needle='#[cfg(test)]\nuse writer::{push_right_aligned_i16, push_right_aligned_u32};';assert needle in t;t=t.replace(needle,'pub(crate) use writer::push_right_aligned_u32;\n#[cfg(test)]\nuse writer::push_right_aligned_i16;',1);t+='\n#[cfg(test)]\n#[path = "mol2000/v3000_integer_tests.rs"]\nmod v3000_integer_tests;\n';p.write_text(t)
p=s/'crates/chematic-mol/src/mol3000/writer.rs';t=p.read_text();needle='''        write!(
            out,
            "M  V30 {} {} ",
            idx.0 + 1,
            v3000_atom_symbol(mol, idx, atom)
        )
        .expect("writing to String cannot fail");''';replacement='''        out.push_str("M  V30 ");
        crate::mol2000::push_right_aligned_u32(out, idx.0 + 1, 0);
        out.push(' ');
        out.push_str(v3000_atom_symbol(mol, idx, atom));
        out.push(' ');''';assert needle in t;t=t.replace(needle,replacement,1)
needle='''        write!(
            out,
            "M  V30 {bond_id} {} {atom1} {atom2}{stereo}",
            v3000_bond_order(written)
        )
        .expect("writing to String cannot fail");''';replacement='''        out.push_str("M  V30 ");
        crate::mol2000::push_right_aligned_u32(out, bond_id, 0);
        out.push(' ');
        out.push(char::from(b'0' + v3000_bond_order(written)));
        out.push(' ');
        crate::mol2000::push_right_aligned_u32(out, atom1, 0);
        out.push(' ');
        crate::mol2000::push_right_aligned_u32(out, atom2, 0);
        out.push_str(stereo);''';assert needle in t;t=t.replace(needle,replacement,1);p.write_text(t)
p=s/'crates/chematic-mol/src/mol2000/v3000_integer_tests.rs';p.write_text('''//! Decimal formatter oracle for padded V2000 and unpadded V3000 fields.
use super::*;

#[test]
fn decimal_table_and_fallback_match_integer_formatter() {
    let mut values: Vec<u32> = (0..=1024).collect();
    for power in [10u32,100,1000,10000,100000,1000000,10000000,100000000,1000000000] {
        values.extend([power-1,power,power+1]);
    }
    values.extend([u32::MAX-1,u32::MAX]);
    let mut seed = 0x72646974u32;
    for _ in 0..10000 {
        seed = seed.wrapping_mul(1664525).wrapping_add(1013904223);
        values.push(seed);
    }
    for width in 0..=12 {
        for &value in &values {
            let mut actual = String::from("prefix:");
            push_right_aligned_u32(&mut actual,value,width);
            assert_eq!(actual,format!("prefix:{value:>width$}"));
        }
    }
}
''')
p=s/'CHANGELOG.md';p.write_text(p.read_text().replace('### Performance\n','### Performance\n\n- Reuse the decimal field table for unpadded V3000 2D atom and bond ids,\n  appending bond-type digits directly while preserving exact output bytes.\n',1));print('Direct integer writer candidate written; original Unicode reader retained')
