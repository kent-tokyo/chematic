from pathlib import Path
import hashlib,json
s=Path('/private/tmp/chematic-openbabel-count-capacity-source-20261015');r=Path('/private/tmp/chematic-openbabel-count-capacity-20261015');p=s/'crates/chematic-mol/src/mol3000.rs';old=(r/'baseline-mol3000.rs').read_text()
a=old.index('    fn in_ctab(');b=old.index('    fn in_atom_block(',a);reference=old[a:b].replace('fn in_ctab(&mut self,','fn reference_in_ctab(parser: &mut V3000Parser,').replace('self.', 'parser.')
t=old.replace('    expected_atoms: usize,','    allocation_limit: usize,\n    expected_atoms: usize,',1)
t=t.replace('            expected_atoms: 0,','            allocation_limit: 0,\n            expected_atoms: 0,',1)
needle='        } else if is_marker(tokens, "BEGIN", "ATOM") {'
t=t.replace(needle,needle+'''
            // Allocate once when entering ATOM, after the final COUNTS row.
            // A declared count is only an allocation hint: bound it by the
            // actual logical row count before trusting potentially huge input.
            let atoms = self.expected_atoms.min(self.allocation_limit);
            let bonds = self.expected_bonds.min(self.allocation_limit);
            self.builder = MoleculeBuilder::with_capacity(atoms, bonds);
            self.atom_idx_map.entries.reserve_exact(atoms);
            self.bond_idx_map.entries.reserve_exact(bonds);
            self.coords.reserve_exact(atoms);
            self.raw_z.reserve_exact(atoms);''',1)
t=t.replace('let mut parser = V3000Parser::new();','let mut parser = V3000Parser::new();\n    parser.allocation_limit = logical_lines.len();',1)
t+='\n#[cfg(test)]\n#[path = "mol3000/count_capacity_tests.rs"]\nmod count_capacity_tests;\n';p.write_text(t)
recipe={'id':'v3000-count-capacity-v1','count_lexemes':['0','1','2','7','9','11','64','1024','4294967295','18446744073709551615','18446744073709551616','-1','+1','','bad'],'allocation_limits':[0,1,5,64],'source_lines':[1,42],'prime':['COUNTS','7','9','11'],'cartesian':'All atom/bond/SGROUP count-lexeme triples for each allocation limit and source line, after the prime COUNTS row; compare original mutation and error precedence','markers':[[],['COUNTS'],['COUNTS','1','2'],['BEGIN','ATOM'],['END','CTAB'],['OTHER']]}
relative='validation/benchmark_corpora/v3000-count-capacity-v1.recipe.json';payload=(json.dumps(recipe,indent=2,sort_keys=True)+'\n').encode();(s/relative).write_bytes(payload)
manifest={'id':'v3000-count-capacity-v1','path':relative,'provenance':{'source':'Original synthetic V3000 count-token and parser-state recipe; no new molecular or third-party data.','retrieved_or_generated':'2026-10-15','transformation':recipe['cartesian']+'; exercise empty, incomplete and block-marker rows separately.'},'license':'MIT OR Apache-2.0','redistribution':{'status':'allowed','rationale':'Original project-authored syntax and bounded-allocation cases.'},'minimization':'15 integer lexemes cover zero, ordinary counts, platform limits, overflow, signs and malformed data. Four row budgets check zero and bounded allocation. A repeated COUNTS prime checks partial mutations and retained error precedence.','oracle':{'name':'Frozen original V3000Parser::in_ctab','version':'chematic 1.1.4 plus nine accepted source candidates','settings':'Exact typed error Debug, state discriminant and all expected count fields; original code is frozen in the test module.','mol3000_sha256':hashlib.sha256(old.encode()).hexdigest()},'expected':{'status':'All count/state/error outcomes equal the original method; reserved vector capacities never exceed the logical-row budget; repeated COUNTS retains empty storage.','absolute_tolerance':0,'units':'record counts and vector capacity','cartesian_cases':15**3*4*2},'sha256':hashlib.sha256(payload).hexdigest()}
(s/'validation/benchmark_corpora/v3000-count-capacity-v1.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
test='''//! Frozen COUNTS transition oracle; see the count-capacity fixture manifest.
use super::*;

'''+reference+'''
fn recipe() -> serde_json::Value {
    serde_json::from_str(include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/v3000-count-capacity-v1.recipe.json"))).unwrap()
}

fn compare_counts(tokens: &[&str], limit: usize, line: usize) {
    let mut expected = V3000Parser::new();
    let mut actual = V3000Parser::new();
    expected.allocation_limit = limit;
    actual.allocation_limit = limit;
    let prime = ["COUNTS", "7", "9", "11"];
    reference_in_ctab(&mut expected, &prime, 1).unwrap();
    actual.in_ctab(&prime, 1).unwrap();
    let left = reference_in_ctab(&mut expected, tokens, line);
    let right = actual.in_ctab(tokens, line);
    assert_eq!(format!(left_debug = "{left:?}"), format!("{right:?}"));
    assert_eq!(actual.expected_atoms, expected.expected_atoms);
    assert_eq!(actual.expected_bonds, expected.expected_bonds);
    assert_eq!(actual.expected_sgroups, expected.expected_sgroups);
    assert_eq!(std::mem::discriminant(&actual.state), std::mem::discriminant(&expected.state));
    assert!(actual.coords.capacity() <= limit);
    assert!(actual.raw_z.capacity() <= limit);
    assert!(actual.atom_idx_map.entries.capacity() <= limit);
    assert!(actual.bond_idx_map.entries.capacity() <= limit);
    assert_eq!(actual.builder.atom_count(), 0);
}

#[test]
fn count_transitions_errors_and_allocation_bounds_match_original() {
    let recipe = recipe();
    let values: Vec<_> = recipe["count_lexemes"].as_array().unwrap().iter().map(|v| v.as_str().unwrap()).collect();
    let mut cases = 0;
    for limit in recipe["allocation_limits"].as_array().unwrap() {
        for line in recipe["source_lines"].as_array().unwrap() {
            for &a in &values {
                for &b in &values {
                    for &sg in &values {
                        compare_counts(&["COUNTS",a,b,sg], limit.as_u64().unwrap() as usize, line.as_u64().unwrap() as usize);
                        cases += 1;
                    }
                }
            }
        }
    }
    assert_eq!(cases,27000);
}

#[test]
fn incomplete_counts_and_block_markers_match_original() {
    let recipe = recipe();
    for row in recipe["markers"].as_array().unwrap() {
        let tokens: Vec<_> = row.as_array().unwrap().iter().map(|v| v.as_str().unwrap()).collect();
        // consume() skips empty rows before dispatching to in_ctab().
        if !tokens.is_empty() { compare_counts(&tokens, 5, 42); }
    }
}
'''
test=test.replace('format!(left_debug = "{left:?}")','format!("{left:?}")')
test+='''
#[test]
fn allocation_waits_for_atom_block_and_is_bounded_by_real_rows() {
    let recipe = recipe();
    for limit in recipe["allocation_limits"].as_array().unwrap() {
        let limit = limit.as_u64().unwrap() as usize;
        let mut parser = V3000Parser::new();
        parser.allocation_limit = limit;
        let large = usize::MAX.to_string();
        for _ in 0..16 {
            parser.in_ctab(&["COUNTS", &large, &large, "0"], 1).unwrap();
            assert_eq!(parser.coords.capacity(), 0);
            assert_eq!(parser.raw_z.capacity(), 0);
            assert_eq!(parser.atom_idx_map.entries.capacity(), 0);
            assert_eq!(parser.bond_idx_map.entries.capacity(), 0);
        }
        parser.in_ctab(&["BEGIN", "ATOM"], 2).unwrap();
        assert_eq!(parser.coords.capacity(), limit);
        assert_eq!(parser.raw_z.capacity(), limit);
        assert_eq!(parser.atom_idx_map.entries.capacity(), limit);
        assert_eq!(parser.bond_idx_map.entries.capacity(), limit);
    }
}
'''
(s/'crates/chematic-mol/src/mol3000/count_capacity_tests.rs').write_text(test)
p=s/'CHANGELOG.md';p.write_text((r/'baseline-CHANGELOG.md').read_text().replace('### Performance\n','### Performance\n\n- Reserve V3000 molecule, coordinate and source-ID storage from declared\n  counts, bounded by the actual logical-row count, while preserving parser\n  state transitions, errors and diagnostics.\n',1))
print('Bounded COUNTS-capacity candidate written')
