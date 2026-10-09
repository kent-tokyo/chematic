fn main() {
    for s in std::env::args().skip(1) {
        let m = chematic_smiles::parse(&s).unwrap();
        println!(
            "{s} {:?}",
            chematic_chem::descriptors::crippen_types_rdkit_model(&m)
        );
    }
}
