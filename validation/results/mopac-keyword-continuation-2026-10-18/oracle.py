"""Original synthetic keyword fixtures, normalized Open Babel geometry controls."""
from pathlib import Path
import json
import subprocess

root = Path(__file__).resolve().parent
cases = [json.loads(line.split('MOPAC_ORACLE ', 1)[1])
         for line in (root / 'after.log').read_text().splitlines()
         if line.startswith('MOPAC_ORACLE ')]
assert len(cases) == 36
summary = dict(cases=36, normalized_geometry_reads=0, coordinate_values_checked=0,
               direct_ordinary_reads=0, direct_continuation_reads=0,
               direct_continuation_failures=0, direct_continuation_failures_with_success_exit=0,
               mismatches=0)
observations = []

def obabel(text):
    return subprocess.run(['/opt/homebrew/bin/obabel', '-imopcrt', '-oxyz'], input=text,
                          capture_output=True, text=True, timeout=15)

def verify_xyz(text, positions):
    lines = text.splitlines()
    assert int(lines[0]) == 2
    assert [line.split()[0] for line in lines[2:]] == ['C', 'O']
    observed = [[float(x) for x in line.split()[1:]] for line in lines[2:]]
    assert len(observed) == 2
    for actual, expected in zip(observed, positions):
        assert all(abs(a-b) <= 1e-12 for a,b in zip(actual,expected)), (actual,expected)

for case in cases:
    count = case['keyword_lines']
    assert case['input'].splitlines()[:count] == case['keywords'].split('\n')
    assert case['output'].splitlines()[:count] == case['keywords'].split('\n')
    for field in ('input', 'output'):
        lines = case[field].splitlines()
        assert lines[count:count+2] == ['original title', 'retained comment']
        # These authored fixtures have no quoted keyword values. This is an
        # explicit control projection, not a parser for general MOPAC keywords.
        flat = ' '.join(token for line in lines[:count] for token in line.split() if token != '+')
        normalized = flat + '\n' + '\n'.join(lines[count:]) + '\n'
        converted = obabel(normalized)
        assert converted.returncode == 0 and converted.stdout, converted.stderr
        verify_xyz(converted.stdout, case['positions'])
        summary['normalized_geometry_reads'] += 1
        summary['coordinate_values_checked'] += 6
        direct = obabel(case[field])
        if count == 1:
            assert direct.returncode == 0 and direct.stdout
            verify_xyz(direct.stdout, case['positions'])
            summary['direct_ordinary_reads'] += 1
        elif not direct.stdout:
            summary['direct_continuation_failures'] += 1
            summary['direct_continuation_failures_with_success_exit'] += int(direct.returncode == 0)
            observations.append(dict(keyword_lines=count, plus_position=case['plus_position'], field=field,
                                     returncode=direct.returncode, stdout=direct.stdout, stderr=direct.stderr))
        else:
            verify_xyz(direct.stdout, case['positions'])
            summary['direct_continuation_reads'] += 1
summary['openbabel_version'] = subprocess.run(['/opt/homebrew/bin/obabel','-V'], capture_output=True, text=True, check=True).stdout.strip()
(root/'oracle-summary.json').write_text(json.dumps(summary, indent=2)+'\n')
(root/'direct-openbabel-diagnostics.json').write_text(json.dumps(observations, indent=2)+'\n')
print(json.dumps(summary, indent=2))
