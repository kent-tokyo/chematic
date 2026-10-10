use crate::{MancudeContext, budget::CipBudget, digraph::CipDigraph};
use chematic_core::AtomIdx;
#[test]
fn digraph_dumps_preserve_node_identity_parentage_and_fractional_labels() {
    for source in ["CC(=O)N", "C1CCCCC1", "N#CC", "c1ccncc1", "c1cc[nH]c1"] {
        let mut mol = chematic_smiles::parse(source).unwrap();
        chematic_perception::kekulize_inplace(&mut mol).unwrap();
        let context = MancudeContext::compute(&mol);
        for root in 0..mol.atom_count() {
            let mut graph = CipDigraph::new_with_mancude(
                &mol,
                AtomIdx(root as u32),
                CipBudget::default_budget(),
                &context,
            )
            .unwrap();
            graph.expand_all(graph.root()).unwrap();
            let json = crate::debug::to_json(&graph);
            let value: serde_json::Value = serde_json::from_str(&json).unwrap();
            assert_eq!(value["root_atom"], root);
            let nodes = value["nodes"].as_array().unwrap();
            assert_eq!(nodes.len(), graph.nodes().len());
            for (actual, node) in nodes.iter().zip(graph.nodes()) {
                assert_eq!(actual["id"], node.id.0);
                assert_eq!(actual["depth"], node.depth);
                assert_eq!(actual["atomic_number"], node.atomic_number.to_string());
                assert_eq!(
                    actual["parent"],
                    node.parent
                        .map(|p| serde_json::json!(p.0))
                        .unwrap_or(serde_json::Value::Null)
                );
            }
            let tree = crate::debug::to_tree_string(&graph);
            assert_eq!(tree.lines().count(), nodes.len());
            assert!(tree.starts_with(&format!("atom={root}\n")));
            assert_eq!(tree, crate::debug::to_tree_string(&graph));
            assert_eq!(json, crate::debug::to_json(&graph));
        }
    }
}

#[test]
fn comparison_traces_explain_rankings_without_changing_the_decision() {
    use crate::{CompareContext, DigraphExpander, compare_ligands, trace::ComparisonTrace};
    for source in [
        "CC(F)Cl",
        "CC(=O)N",
        "C1CCCCC1",
        "c1ccncc1",
        "[13CH3]C(C)(F)Cl",
    ] {
        let mut mol = chematic_smiles::parse(source).unwrap();
        chematic_perception::kekulize_inplace(&mut mol).unwrap();
        for root in 0..mol.atom_count() {
            let mut graph =
                CipDigraph::new(&mol, AtomIdx(root as u32), CipBudget::default_budget()).unwrap();
            let children = graph.children(graph.root()).unwrap();
            for (i, &left) in children.iter().enumerate() {
                for &right in &children[i + 1..] {
                    let expected =
                        compare_ligands(&mut graph, left, right, &mut CompareContext::default())
                            .unwrap();
                    let mut trace = ComparisonTrace::new(left, right);
                    let actual = compare_ligands(
                        &mut graph,
                        left,
                        right,
                        &mut CompareContext::with_trace(&mut trace),
                    )
                    .unwrap();
                    assert_eq!(actual, expected);
                    assert!(!trace.decisions.is_empty());
                    let output = trace.to_string();
                    assert!(
                        output.starts_with(&format!(
                            "comparing root {} vs root {}:",
                            left.0, right.0
                        ))
                    );
                    for step in &trace.decisions {
                        assert!(step.left_node.0 < graph.nodes().len() as u32);
                        assert!(step.right_node.0 < graph.nodes().len() as u32);
                        assert!(output.contains(&step.to_string()));
                        assert!(step.to_string().contains(&step.rule));
                    }
                }
            }
            assert_eq!(graph.edges().len() + 1, graph.nodes().len());
            for (index, edge) in graph.edges().iter().enumerate() {
                assert_eq!(edge.id.0 as usize, index);
                assert_eq!(graph.node(edge.child).parent, Some(edge.parent));
                assert_eq!(
                    graph.node(edge.child).depth,
                    graph.node(edge.parent).depth + 1
                );
            }
        }
    }
    let mol = chematic_smiles::parse("C(F)Cl").unwrap();
    let mut graph = CipDigraph::new(&mol, AtomIdx(0), CipBudget::default_budget()).unwrap();
    let children = graph.children(graph.root()).unwrap();
    let mut ctx = CompareContext::default();
    ctx.max_recursive_calls = 0;
    let error = compare_ligands(&mut graph, children[0], children[1], &mut ctx).unwrap_err();
    assert!(error.to_string().contains("budget"));
}
