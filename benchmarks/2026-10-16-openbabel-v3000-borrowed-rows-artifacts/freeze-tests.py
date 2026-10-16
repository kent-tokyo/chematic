from pathlib import Path
import hashlib, json, re
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
OUT=Path(__file__).resolve().parent
original=(OUT/'baseline-source/crates/chematic-mol/src/mol3000.rs').read_text()
name='read_mol_v3000_with_diagnostics'
start=original.index('pub fn '+name+'(');at=original.index('{',start);depth=1;end=at+1
while depth:depth+=(original[end]=='{')-(original[end]=='}');end+=1
frozen=original[start:end].replace('pub fn '+name,'fn frozen_collected_rows',1)
previous=(ROOT/'crates/chematic-mol/src/mol3000/token_buffer_tests.rs').read_text()
begin=previous.index('fn compare(')
tests=previous[begin:].replace('frozen_per_line_tokens(input)','frozen_collected_rows(input)').replace('assert_eq!(old.to_string(), new.to_string());','assert_eq!(old, new);')
tests=tests.replace('fn compare(input: &str) -> bool {','fn compare(input: &str) -> bool {')
tests=tests.replace('long_short_empty_and_continued_rows_keep_full_reports','plain_and_continued_rows_keep_full_reports').replace('token storage changed parse success/error boundary','borrowed rows changed parse success/error boundary')
extra=r'''

#[test]
fn all_pinned_v3000_rows_preserve_full_reports() {
    use flate2::read::GzDecoder;
    use std::io::Read;
    let corpora: [(&[u8], usize); 3] = [
        (include_bytes!("../../../../validation/benchmark_corpora/openbabel-v3000-speed-slice0-1000-v1.molblocks.gz"),1000),
        (include_bytes!("../../../../validation/benchmark_corpora/openbabel-v3000-speed-slice5000-1000-v1.molblocks.gz"),1000),
        (include_bytes!("../../../../validation/benchmark_corpora/openbabel-v3000-speed-chains64-v1.molblocks.gz"),64),
    ];
    let mut rows=0;let mut plain=0;
    for (corpus, expected) in corpora {
        let mut input=String::new();GzDecoder::new(corpus).read_to_string(&mut input).unwrap();
        let mut count=0;
        for block in input.split("$$$$\n").filter(|block| !block.trim().is_empty()) {
            let lines=validate_v3000_physical_lines(block).unwrap();
            if let Some(n)=plain_v3000_line_count(&lines) {
                assert_eq!(n,collect_v30_lines(&lines).len());plain+=1;
            }
            for _ in 0..2 { assert!(compare(block)); }
            count+=1;rows+=1;
        }
        assert_eq!(count,expected);
    }
    assert_eq!(rows,2064);assert!(plain>2000);
}

#[test]
fn allocation_hints_and_continuation_selection_match_original_rows() {
    for continuation in [false,true] {
        let input=chain(3," ",continuation);
        let lines=validate_v3000_physical_lines(&input).unwrap();
        assert_eq!(plain_v3000_line_count(&lines).is_some(),!continuation);
        let mut old=V3000Parser::new();let logical=collect_v30_lines(&lines);
        old.allocation_limit=logical.len();
        for line in &logical {let tokens=line.payload.split_whitespace().collect::<Vec<_>>();old.consume(line,&tokens).unwrap();}
        let mut new=V3000Parser::new();consume_v3000_rows(&mut new,&lines).unwrap();
        assert_eq!(new.allocation_limit,old.allocation_limit);
        assert_eq!(new.coords.capacity(),old.coords.capacity());
        assert_eq!(new.raw_z.capacity(),old.raw_z.capacity());
        assert_eq!(new.atom_idx_map.entries.capacity(),old.atom_idx_map.entries.capacity());
        assert_eq!(new.bond_idx_map.entries.capacity(),old.bond_idx_map.entries.capacity());
    }
    let input=chain(3," ",false);
    for header in ["M  V30 title", "M  V30 title-", "M  V30 ", "plain title-"] {
        let (_,tail)=input.split_once('\n').unwrap();compare(&format!("{header}\n{tail}"));
    }
    for counts in ["0 0 0 0 0","999999999 999999999 0 0 0","3 2 0 0 0"] {
        compare(&input.replace("COUNTS 3 2 0 0 0",&format!("COUNTS {counts}")));
    }
}
'''
dest=ROOT/'crates/chematic-mol/src/mol3000/record_stream_tests.rs'
dest.write_text('//! Frozen main reader: full reports and exact errors must be preserved.\nuse super::*;\n\n'+frozen+'\n\n'+tests+extra)
(OUT/'frozen-reader-proof.json').write_text(json.dumps({'base_revision':json.loads((OUT/'measurement-plan.json').read_text())['base_revision'],'original_function_sha256':hashlib.sha256(original[start:end].encode()).hexdigest(),'only_name_and_visibility_changed':True},indent=2)+'\n')
print('Frozen full original reader and complete report/error comparison')
