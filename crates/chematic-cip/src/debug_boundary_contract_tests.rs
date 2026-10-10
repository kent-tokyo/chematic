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
