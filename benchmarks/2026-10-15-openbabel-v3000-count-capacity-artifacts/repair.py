from pathlib import Path
r=Path('/private/tmp/chematic-openbabel-count-capacity-20261015');s=Path('/private/tmp/chematic-openbabel-count-capacity-source-20261015');(r/'pre-repair-mol3000.rs').write_bytes((s/'crates/chematic-mol/src/mol3000.rs').read_bytes());(r/'native-library-tests.log').rename(r/'native-library-tests-constructor-error.log');(r/'validation-pipeline.log').rename(r/'validation-constructor-error.log')
p=r/'optimize.py';t=p.read_text().replace("old=p.read_text();assert old==(r/'baseline-mol3000.rs').read_text()","old=(r/'baseline-mol3000.rs').read_text()")
t=t.replace("t=t.replace('impl V3000Parser {\\n    fn new() -> Self {','impl V3000Parser {\\n    fn new(allocation_limit: usize) -> Self {',1)\n",'')
t=t.replace("'            allocation_limit,\\n            expected_atoms: 0,'","'            allocation_limit: 0,\\n            expected_atoms: 0,'")
t=t.replace("needle='            self.expected_sgroups = Some(parse_v3000_count(tokens[3], line, \"SGROUP\")?);'","needle='        } else if is_marker(tokens, \"BEGIN\", \"ATOM\") {'")
t=t.replace('COUNTS precedes every atom/bond row, so these stores are empty.','Allocate once when entering ATOM, after the final COUNTS row.')
t=t.replace("t=t.replace('let mut parser = V3000Parser::new();','let mut parser = V3000Parser::new(logical_lines.len());',1)","t=t.replace('let mut parser = V3000Parser::new();','let mut parser = V3000Parser::new();\\n    parser.allocation_limit = logical_lines.len();',1)")
t=t.replace('V3000Parser::new(limit);','V3000Parser::new();')
t=t.replace('    let prime = ["COUNTS", "7", "9", "11"];','    expected.allocation_limit = limit;\n    actual.allocation_limit = limit;\n    let prime = ["COUNTS", "7", "9", "11"];')
t=t.replace("p=s/'CHANGELOG.md';p.write_text(p.read_text().replace", "p=s/'CHANGELOG.md';p.write_text((r/'baseline-CHANGELOG.md').read_text().replace")
# Add a direct capacity assertion that exercises the transition where reservation occurs.
needle="(s/'crates/chematic-mol/src/mol3000/count_capacity_tests.rs').write_text(test)"
extra="""test+='''
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
"""
t=t.replace(needle,extra+needle);p.write_text(t)
# Recreate only the intended private candidate files from the frozen baseline.
import subprocess
subprocess.run(['python3',str(p)],check=True)
v=(r/'validate.py').read_text();prefix=v[:v.index('while not')];commands=v[v.index('commands='):];commands=commands.replace("('optimize.log',['python3',str(r/'optimize.py')]),",'')
(r/'resume-validation.py').write_text(prefix+commands)
print('Constructor compatibility restored; allocation deferred to a single ATOM transition')
