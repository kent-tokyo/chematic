use super::*;
use chematic_core::{Atom, BondOrder, validate_valence};

#[test]
fn valence_warnings_do_not_replace_structure_change_status_or_mutate_input() {
    // Deliberately invalid valence on an otherwise valid public graph:
    // a validation warning must remain visible through the pipeline.
    let mut builder = MoleculeBuilder::new();
    let carbon = builder.add_atom(Atom::new(Element::C));
    for _ in 0..5 {
        let fluorine = builder.add_atom(Atom::new(Element::F));
        builder
            .add_bond(carbon, fluorine, BondOrder::Single)
            .unwrap();
    }
    let mut sodium = Atom::new(Element::NA);
    sodium.charge = 1;
    builder.add_atom(sodium);
    let input = builder.build();
    let before = MoleculeSnapshot::from_mol(&input);
    assert_eq!(validate_valence(&input).len(), 1);
    for largest_fragment_only in [false, true] {
        let pipeline = StandardizationPipeline::new(StandardizeOptions {
            canonical_tautomer: false,
            neutralize_charges: false,
            remove_explicit_h: false,
            largest_fragment_only,
            zwitterion_handling: ZwitterionHandling::Keep,
        });
        assert_eq!(
            pipeline.options().largest_fragment_only,
            largest_fragment_only
        );
        let (output, report) = pipeline.run(&input);
        assert_eq!(report.input, before);
        assert_eq!(report.output, MoleculeSnapshot::from_mol(&output));
        assert_eq!(report.changed(), largest_fragment_only);
        assert_eq!(
            report.status,
            if largest_fragment_only {
                PipelineStatus::CompletedWithWarnings
            } else {
                PipelineStatus::Unchanged
            }
        );
        assert_eq!(
            output.atom_count(),
            if largest_fragment_only { 6 } else { 7 }
        );
        assert_eq!(output.bond_count(), 5);
        assert_eq!(validate_valence(&output).len(), 1);
        assert_eq!(report.warnings[0].code, "input_valence_validation_failed");
        assert!(
            report.warnings[0]
                .message
                .contains("1 valence validation issue")
        );
        let staged: Vec<_> = report
            .warnings
            .iter()
            .filter(|w| w.code == "valence_validation_failed")
            .collect();
        assert_eq!(staged.len(), if largest_fragment_only { 2 } else { 1 });
        for warning in staged {
            assert!(
                warning
                    .message
                    .contains("produced 1 valence validation issue")
            );
        }
        assert_eq!(MoleculeSnapshot::from_mol(&input), before);
    }
}
