from pathlib import Path
import json,hashlib
s=Path('/private/tmp/chematic-openbabel-count-zero-map-source-20261015');r=Path('/private/tmp/chematic-openbabel-count-zero-map-20261015');p=s/'crates/chematic-mol/src/mol3000.rs';t=p.read_text();t=t.replace('self.expected_atoms.min(self.allocation_limit)','v3000_storage_capacity(self.expected_atoms, self.allocation_limit)').replace('self.expected_bonds.min(self.allocation_limit)','v3000_storage_capacity(self.expected_bonds, self.allocation_limit)')
pos=t.index('fn parse_v3000_count(')
t=t[:pos]+'''// Preserve normal Vec growth capacities where the real input budget allows.
// Storage shape after parsing can affect subsequent prepared operations.
fn v3000_storage_capacity(count: usize, limit: usize) -> usize {
    let bounded = count.min(limit);
    if bounded == 0 {
        0
    } else {
        bounded
            .checked_next_power_of_two()
            .unwrap_or(limit)
            .max(4)
            .min(limit)
    }
}

'''+t[pos:];p.write_text(t)
p=s/'crates/chematic-mol/src/mol3000/count_capacity_tests.rs';p.write_text(p.read_text()+'''
#[test]
fn rounded_capacity_matches_original_growth_with_sufficient_row_budget() {
    for count in 0..=1024 {
        let mut original = Vec::<(f64, f64)>::new();
        for _ in 0..count { original.push((0.0,0.0)); }
        assert_eq!(v3000_storage_capacity(count, count*2+7), original.capacity());
    }
    assert_eq!(v3000_storage_capacity(usize::MAX, 0), 0);
    assert_eq!(v3000_storage_capacity(usize::MAX, 5), 5);
    assert_eq!(v3000_storage_capacity(usize::MAX, usize::MAX), usize::MAX);
}
''')
p=s/'validation/benchmark_corpora/v3000-count-capacity-v2.recipe.json';recipe=json.loads(p.read_text());recipe['capacity_oracle']={'counts_inclusive':[0,1024],'row_budget':'2 * count + 7','reference':'Push count zero XY tuples into an initially empty Vec and compare its capacity; synthetic storage-only values'};p.write_text(json.dumps(recipe,indent=2,sort_keys=True)+'\n')
p=s/'validation/benchmark_corpora/v3000-count-capacity-v2.manifest.json';m=json.loads(p.read_text());m['expected']['capacity_growth_cases']=1025;m['sha256']=hashlib.sha256((s/m['path']).read_bytes()).hexdigest();p.write_text(json.dumps(m,indent=2)+'\n')
p=s/'CHANGELOG.md';p.write_text(p.read_text().replace('counts, bounded by the actual logical-row count, while preserving parser','counts with geometric-growth capacities bounded by the actual logical-row\n  count, while preserving parser'))
print('Round reservations toward existing growth capacities; chemistry unchanged')
