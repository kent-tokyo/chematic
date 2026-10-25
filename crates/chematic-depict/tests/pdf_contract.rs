#![cfg(feature = "pdf")]

use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};
use chematic_depict::{
    RenderOptions, depict_pdf, depict_pdf_opts, depict_reaction_svg, depict_svg, depict_svg_grid,
    svg_to_pdf,
};
use chematic_rxn::parse_reaction;
use chematic_smiles::parse;
use lopdf::{Document, Object, content::Content};

fn document(bytes: &[u8]) -> Document {
    assert!(bytes.starts_with(b"%PDF-"));
    let pdf = Document::load_mem(bytes).expect("PDF should be readable");
    assert_eq!(pdf.get_pages().len(), 1, "depiction should have one page");
    pdf
}

fn drawing_operators(pdf: &Document) -> Vec<String> {
    let mut operators: Vec<String> = pdf
        .get_pages()
        .values()
        .flat_map(|page_id| {
            pdf.get_and_decode_page_content(*page_id)
                .expect("page content should decode")
                .operations
                .into_iter()
                .map(|operation| operation.operator)
        })
        .collect();
    operators.extend(
        pdf.objects
            .values()
            .filter_map(|object| match object {
                Object::Stream(stream) => {
                    let is_form = stream
                        .dict
                        .get(b"Subtype")
                        .ok()
                        .and_then(|value| value.as_name().ok())
                        == Some(b"Form".as_slice());
                    if is_form { Some(stream) } else { None }
                }
                _ => None,
            })
            .filter_map(|stream| stream.decompressed_content().ok())
            .filter_map(|bytes| Content::decode(&bytes).ok())
            .flat_map(|content| {
                content
                    .operations
                    .into_iter()
                    .map(|operation| operation.operator)
            })
            .collect::<Vec<_>>(),
    );
    operators
}

#[test]
fn molecule_pdf_keeps_vector_geometry_and_atom_labels() {
    for smiles in ["c1ccccc1", "CCO", "[NH4+]"] {
        let molecule = parse(smiles).unwrap();
        let pdf = document(&depict_pdf(&molecule));
        let operators = drawing_operators(&pdf);
        assert!(
            operators
                .iter()
                .any(|operator| matches!(operator.as_str(), "S" | "f" | "f*")),
            "{smiles}: expected vector drawing in PDF content, got {operators:?}"
        );
        if smiles == "c1ccccc1" {
            assert!(operators.iter().any(|operator| operator == "S"));
        } else {
            assert!(
                operators.iter().any(|operator| operator == "Tj"),
                "{smiles}: atom labels should remain PDF text"
            );
        }
    }
}

#[test]
fn query_bond_keeps_a_dash_pattern() {
    let mut builder = MoleculeBuilder::new();
    let first = builder.add_atom(Atom::new(Element::C));
    let second = builder.add_atom(Atom::new(Element::C));
    builder
        .add_bond(first, second, BondOrder::QueryAny)
        .unwrap();
    let pdf = document(&depict_pdf(&builder.build()));
    let page_id = *pdf.get_pages().values().next().unwrap();
    let content = pdf.get_and_decode_page_content(page_id).unwrap();
    assert!(content.operations.iter().any(|operation| {
        operation.operator == "d"
            && operation
                .operands
                .first()
                .and_then(|value| value.as_array().ok())
                .is_some_and(|pattern| !pattern.is_empty())
    }));
}

#[test]
fn svg_conversion_matches_molecule_entry_point() {
    let molecule = parse("CCO").unwrap();
    let from_molecule = document(&depict_pdf(&molecule));
    let from_svg = document(&svg_to_pdf(&depict_svg(&molecule)));
    assert_eq!(from_molecule.get_pages().len(), from_svg.get_pages().len());
    assert!(!drawing_operators(&from_svg).is_empty());

    let benzene = parse("c1ccccc1").unwrap();
    let grid = depict_svg_grid(&[&molecule, &benzene], 2);
    let reaction = depict_reaction_svg(&parse_reaction("CC>>CCO").unwrap());
    for svg in [grid, reaction] {
        let pdf = document(&svg_to_pdf(&svg));
        assert!(
            drawing_operators(&pdf)
                .iter()
                .any(|operator| operator == "S")
        );
    }
}

#[test]
fn explicit_dimensions_survive_pdf_conversion() {
    let molecule = parse("CCO").unwrap();
    let opts = RenderOptions {
        width: Some(320),
        height: Some(240),
        ..Default::default()
    };
    let pdf = document(&depict_pdf_opts(&molecule, &opts));
    let page_id = *pdf.get_pages().values().next().unwrap();
    let page = pdf.get_object(page_id).unwrap().as_dict().unwrap();
    let media_box = page.get(b"MediaBox").unwrap().as_array().unwrap();
    let size: Vec<f32> = media_box
        .iter()
        .map(|value| value.as_float().unwrap())
        .collect();
    assert_eq!(size, vec![0.0, 0.0, 320.0, 240.0]);
}
