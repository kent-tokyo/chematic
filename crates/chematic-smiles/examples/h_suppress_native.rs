fn main() {
    for s in std::env::args().skip(1) {
        let m = chematic_smiles::parse(&s).unwrap();
        let h = chematic_smiles::rdkit_hydrogen_suppressed(&m).unwrap();
        println!(
            "{s}\tnative in: {}\tnative out: {}\trdkit out: {:?}",
            chematic_smiles::write(&m),
            chematic_smiles::write(&h),
            chematic_smiles::rdkit_canonical_smiles(&h)
        );
        for (i, a) in h.atoms() {
            if a.chirality != chematic_core::Chirality::None {
                println!(
                    "  atom {} {:?} order {:?}",
                    i.0,
                    a.chirality,
                    h.stereo_neighbor_order(i)
                );
            }
        }
    }
}
