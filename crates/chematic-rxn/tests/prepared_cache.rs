//! Load behaviour of the process-wide prepared-template cache
//! (`PreparedReaction::shared`). One test per binary run, so the global miss
//! counter is not shared with other tests.

use chematic_rxn::PreparedReaction;
use chematic_rxn::transform::prepared_cache_misses;

fn template(i: usize) -> String {
    // Distinct, valid SMIRKS: an isotope label on the product atom.
    format!("[C:1]>>[{}C:1]", i + 1)
}

#[test]
fn hot_templates_survive_cold_traffic_and_threads() {
    let mol = chematic_smiles::parse("CC").unwrap();

    // A hot set of 300 templates used every round, against 5,000 one-off
    // templates (more than twice the 2,048 capacity). A cache cleared whole
    // when full would prepare the hot set again after every 2,048 new
    // templates; the generational shards prepare each template once.
    let before = prepared_cache_misses();
    let hot: Vec<String> = (0..300).map(template).collect();
    let mut cold = (1000..6000).map(template);
    for _round in 0..100 {
        for t in &hot {
            let p = PreparedReaction::shared(t).unwrap();
            assert!(!p.run_reactants(&[&mol]).unwrap().is_empty());
        }
        for t in cold.by_ref().take(50) {
            PreparedReaction::shared(&t).unwrap();
        }
    }
    assert_eq!(prepared_cache_misses() - before, 300 + 5000);

    // Eight threads over one shared set of 1,000 templates (within capacity):
    // every result is right, and each template is prepared about once (two
    // threads can both prepare a template neither found).
    let before = prepared_cache_misses();
    let shared: Vec<String> = (10_000..11_000).map(template).collect();
    std::thread::scope(|scope| {
        for offset in 0..8 {
            let shared = &shared;
            scope.spawn(move || {
                let mol = chematic_smiles::parse("CC").unwrap();
                for round in 0..5 {
                    for k in 0..shared.len() {
                        let t = &shared[(k + offset * 125 + round) % shared.len()];
                        let products = PreparedReaction::shared(t)
                            .unwrap()
                            .run_reactants(&[&mol])
                            .unwrap();
                        assert!(!products.is_empty());
                    }
                }
            });
        }
    });
    let misses = prepared_cache_misses() - before;
    assert!(
        (1000..=1100).contains(&misses),
        "{misses} misses for 1,000 templates"
    );
}
