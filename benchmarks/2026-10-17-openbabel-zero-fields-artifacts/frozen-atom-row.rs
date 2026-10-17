fn write_atom_row(
    out: &mut String,
    mol: &Molecule,
    idx: AtomIdx,
    atom: &Atom,
    coords: Option<(f64, f64)>,
) {
    let symbol = mdl_atom_symbol(mol, idx, atom);
    let charge = encode_charge(atom.charge);
    let mass = encode_mass_difference(atom.element, atom.isotope).unwrap_or(0);
    let atom_map = atom.atom_map.unwrap_or(0);
    let valence = valence_field(mol, idx, atom);
    if let Some((x, y)) = coords {
        push_fixed4(out, x, 10);
        push_fixed4(out, y, 10);
        out.push_str("    0.0000 ");
    } else {
        out.push_str("    0.0000    0.0000    0.0000 ");
    }
    out.push_str(symbol);
    for _ in symbol.len()..3 {
        out.push(' ');
    }
    if mass == 0 && charge == 0 {
        out.push_str(" 0  0");
    } else {
        push_right_aligned_i16(out, mass, 2);
        push_right_aligned_u32(out, charge as u32, 3);
    }
    out.push_str("  0  0  0");
    push_right_aligned_u32(out, u32::from(valence), 3);
    out.push_str("  0  0  0");
    push_right_aligned_u32(out, atom_map as u32, 3);
    out.push_str("  0\n");
}
