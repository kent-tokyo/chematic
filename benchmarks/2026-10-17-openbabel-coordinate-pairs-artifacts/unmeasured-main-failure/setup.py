from pathlib import Path
import shutil,subprocess,json,hashlib
R=Path('/private/tmp/chematic-openbabel-writer-20261017');O=Path(__file__).resolve().parent
base=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip()
B=O/'baseline-source';B.mkdir()
for d in ['crates','tools']:shutil.copytree(R/d,B/d)
# When replaying, start from the named base before adding the candidate files.
for f in ['Cargo.toml','Cargo.lock','README.md','LICENSE','LICENSE-MIT','LICENSE-APACHE']:
 if (R/f).exists():shutil.copy2(R/f,B/f)
p=R/'crates/chematic-mol/src/mol2000/writer.rs';s=p.read_text();start=s.index('#[inline]\n/// Append `x`');end=s.index('\npub(crate) fn push_right_aligned_u32',start)
frozen=s[start:end].replace('pub(crate) fn push_fixed4','fn frozen_push_fixed4',1)
p.write_text(s[:start]+'#[path = "fixed4.rs"]\nmod fixed4;\npub(crate) use fixed4::push_fixed4;\n'+s[end:])
tests='''use super::push_fixed4;
FROZEN
fn check(x: f64) {
    for width in [0, 1, 6, 10, 24, 30] {
        let mut actual = String::from("prefix:");
        let mut frozen = actual.clone();
        push_fixed4(&mut actual, x, width);
        frozen_push_fixed4(&mut frozen, x, width);
        assert_eq!(actual, frozen, "bits={:016x}, width={width}", x.to_bits());
        assert_eq!(actual, format!("prefix:{x:>width$.4}"), "bits={:016x}", x.to_bits());
    }
}
#[test]
fn decimal_pairs_match_frozen_and_standard_formatter() {
    for n in 0..10_000 {
        for x in [n as f64 / 10_000.0, -(n as f64 / 10_000.0), n as f64 + 0.0001] { check(x); }
    }
    let mut bits=0x5678_1234_abcd_9876u64;
    for _ in 0..100_000 {
        bits ^= bits << 13; bits ^= bits >> 7; bits ^= bits << 17;
        check(f64::from_bits(bits));
    }
}
#[test]
fn rounding_edges_signs_and_padding_match() {
    for x in [0.0, -0.0, f64::INFINITY, f64::NEG_INFINITY, f64::NAN, -f64::NAN,
              f64::MIN_POSITIVE, f64::MAX, -f64::MAX, f64::from_bits(1)] {check(x);}
    for n in [0u64,1,9,10,99,100,999,1000,9999,10000,99999,1_000_000,99_999_999_999] {
        for fraction in [0.0,0.00005,0.00995,0.09995,0.99995] {
            let x=n as f64 + fraction;
            for value in [x, x.next_up(), x.next_down(), -x, (-x).next_up(), (-x).next_down()] {check(value);}
        }
    }
}
'''.replace('FROZEN',frozen)
(R/'crates/chematic-mol/src/mol2000/fixed4_tests.rs').write_text(tests)
(O/'source-proof.json').write_text(json.dumps({'base_commit':base,'frozen_helper_sha256':hashlib.sha256(frozen.replace('fn frozen_push_fixed4','pub(crate) fn push_fixed4',1).encode()).hexdigest(),'transformation':'Convert the fixed4 integer and fractional digits in decimal pairs; include common field padding in the same stack buffer; unchanged rounding and fallback'},indent=2)+'\n')
for name in ['build.py','native.py','prepare.py','formal.py']:
 text=(Path('/private/tmp/chematic-stereo-shape-evidence-20261017')/name).read_text().replace('/private/tmp/chematic-stereo-shape-filter-20261017',str(R))
 if name=='native.py':text=text.replace("packages=['chematic-perception','chematic-mol','chematic-chem','chematic-fp']","packages=['chematic-mol']")
 if name=='prepare.py':
  text=text.replace("'Move unchanged local output-shape rejection before expensive substituent signatures in stereo_centers'","'Decimal pairs and stack-buffer padding in shared fixed4 coordinate formatter'")
  text=text.replace("['v3000_3d-branched-spiral','v3000_3d-branched-offset']","['sdf-0','sdf-5000']")
  text=text.replace('Exact source snapshots and rewrites, frozen classifier matches on 14072 graphs and all 6144 five-atom/H-override graphs; both primary parse paired CI lower>1;', 'Exact source snapshots and rewrites, frozen helper and standard formatter match decimal pairs, random f64 bits and directed boundaries; both primary write paired CI lower>1;')
 if name=='formal.py':text=text.replace("mode=='parse': primary[key]=ratio","mode=='write': primary[key]=ratio")
 (O/name).write_text(text)
print(base)
