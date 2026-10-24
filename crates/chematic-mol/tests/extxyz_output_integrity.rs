use chematic_core::Element;
use chematic_mol::{
    ExtxyzWriter, XyzAtom, XyzFrame, XyzProperty, XyzPropertyKind as K, XyzValue as V,
    parse_extxyz, write_extxyz,
};
fn frame() -> XyzFrame {
    XyzFrame {
        atoms: vec![XyzAtom {
            element: Element::C,
            x: -0.,
            y: 1.25,
            z: 2.5,
        }],
        comment: String::new(),
        lattice: None,
        properties: vec![],
        info: vec![],
    }
}
fn rejects_without_writing(f: &XyzFrame) {
    let mut sink = b"prefix".to_vec();
    let error = ExtxyzWriter::new(&mut sink).write_frame(f).unwrap_err();
    assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
    assert_eq!(sink, b"prefix");
    assert!(write_extxyz(f).is_err());
}
fn property(kind: K, count: usize, values: Vec<Vec<V>>) -> XyzProperty {
    XyzProperty {
        name: "data".into(),
        kind,
        count,
        values,
    }
}
#[test]
fn invalid_shapes_types_tokens_and_numbers_fail_before_the_sink() {
    let bad = [
        property(K::Real, 1, vec![]),
        property(K::Real, 1, vec![vec![V::Real(1.)], vec![V::Real(2.)]]),
        property(K::Real, 2, vec![vec![V::Real(1.)]]),
        property(K::Real, 1, vec![vec![V::Integer(1)]]),
        property(K::Real, 1, vec![vec![V::Real(f64::NAN)]]),
        property(K::Real, usize::MAX, vec![vec![]]),
        property(K::String, 1, vec![vec![V::Str("two words".into())]]),
        property(K::String, 1, vec![vec![V::Str("".into())]]),
        property(K::Logical, 1, vec![vec![V::Str("T".into())]]),
    ];
    for p in bad {
        let mut f = frame();
        f.properties.push(p);
        rejects_without_writing(&f);
    }
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut f = frame();
        f.atoms[0].x = value;
        rejects_without_writing(&f);
        let mut f = frame();
        f.lattice = Some([value; 9]);
        rejects_without_writing(&f);
    }
}
#[test]
fn duplicates_reserved_names_and_stale_comments_are_refused() {
    let p = property(K::Real, 1, vec![vec![V::Real(1.)]]);
    let mut f = frame();
    f.properties = vec![p.clone(), p.clone()];
    rejects_without_writing(&f);
    for name in ["species", "pos", "positions", "symbols", "numbers", "Z"] {
        let mut f = frame();
        let mut p = p.clone();
        p.name = name.into();
        f.properties.push(p);
        rejects_without_writing(&f);
    }
    let mut f = frame();
    let mut p1 = p.clone();
    p1.name = "charge".into();
    let mut p2 = p.clone();
    p2.name = "charges".into();
    f.properties = vec![p1, p2];
    rejects_without_writing(&f);
    let mut f = frame();
    f.info = vec![("note".into(), "a".into()), ("note".into(), "b".into())];
    rejects_without_writing(&f);
    let mut f = frame();
    f.info = vec![("empty".into(), "".into())];
    rejects_without_writing(&f);
    for comment in [
        "line\nbreak",
        "x=y",
        "Properties=species:S:1:pos:R:3:data:R:1",
    ] {
        let mut f = frame();
        f.comment = comment.into();
        rejects_without_writing(&f);
    }
    let mut f = frame();
    f.comment = "Properties=species:S:1:pos:R:3".into();
    assert!(write_extxyz(&f).is_ok());
}
#[test]
fn typed_columns_and_binary64_values_survive_roundtrip() {
    let mut f = frame();
    f.atoms[0].y = f64::from_bits(1);
    f.atoms[0].z = f64::MAX;
    f.lattice = Some([-0., 1., 2., 3., 4., 5., 6., 7., 8.]);
    f.properties = vec![property(
        K::Real,
        2,
        vec![vec![V::Real(-0.), V::Real(f64::from_bits(1))]],
    )];
    for (name, kind, value) in [
        ("integer", K::Integer, V::Integer(i64::MIN)),
        ("string", K::String, V::Str("C.ar".into())),
        ("logical", K::Logical, V::Logical(true)),
    ] {
        f.properties.push(XyzProperty {
            name: name.into(),
            kind,
            count: 1,
            values: vec![vec![value]],
        });
    }
    f.info = vec![
        ("path".into(), "alpha\\beta".into()),
        ("note".into(), "O'Brien".into()),
    ];
    let parsed = parse_extxyz(&write_extxyz(&f).unwrap()).unwrap();
    assert_eq!(parsed.atoms[0].x.to_bits(), f.atoms[0].x.to_bits());
    assert_eq!(parsed.atoms[0].y.to_bits(), 1);
    assert_eq!(parsed.atoms[0].z, f64::MAX);
    assert_eq!(
        parsed.lattice.unwrap().map(f64::to_bits),
        f.lattice.unwrap().map(f64::to_bits)
    );
    assert_eq!(parsed.properties, f.properties);
    assert_eq!(parsed.info, f.info);
}
#[test]
fn valid_sink_failures_still_propagate_and_zero_atom_columns_are_valid() {
    struct Failing;
    impl std::io::Write for Failing {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("sink failed"))
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    assert_eq!(
        ExtxyzWriter::new(Failing)
            .write_frame(&frame())
            .unwrap_err()
            .kind(),
        std::io::ErrorKind::Other
    );
    let mut f = frame();
    f.atoms.clear();
    f.properties.push(property(K::Real, 2, vec![]));
    assert_eq!(
        parse_extxyz(&write_extxyz(&f).unwrap()).unwrap().properties[0].count,
        2
    );
}
