"""Independent observer for the bounded Cartesian ORCA writer profile."""
from __future__ import annotations

import math
import re
import struct


def require(condition, detail):
    if not condition:
        raise ValueError(detail)


def observe(text):
    """Observe emitted settings/rows, without using either chemistry engine."""
    lines = text.splitlines()
    result = dict(comments=[], keywords=[], blocks=[])
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        i += 1
        if not line:
            continue
        if line.startswith('#'):
            result['comments'].append(line)
        elif line.startswith('!'):
            result['keywords'].extend(line[1:].split('#')[0].split())
        elif line.startswith('%'):
            first = line[1:].split(maxsplit=1)
            name, rest = first[0], first[1] if len(first) == 2 else ''
            if rest == 'end':
                result['blocks'].append(dict(name=name, raw='', has_end=True))
            elif rest:
                result['blocks'].append(dict(name=name, raw=rest, has_end=False))
            else:
                window = []
                while i < len(lines) and not lines[i].lstrip().startswith(('%', '!', '*')):
                    window.append(lines[i]); i += 1
                # Empty physical lines between sections are insignificant here.
                while window and not window[-1].strip():
                    window.pop()
                require(window and window[-1].strip() == 'end', 'invalid block end')
                result['blocks'].append(dict(name=name, raw='\n'.join(window[:-1]).strip(), has_end=True))
        else:
            require(line.startswith('* xyz '), 'unexpected section')
            result['coords'], i = geometry(lines, i - 1)
            require(all(not row.strip() for row in lines[i:]), 'extra data after geometry')
            break
    require('coords' in result, 'missing geometry')
    unit_words = [w.lower() for w in result['keywords'] if w.lower() in ('angs', 'bohrs')]
    require(len(set(unit_words)) <= 1, 'unit conflict')
    result['units'] = 'bohr' if 'bohrs' in unit_words else 'angstrom'
    return result


def geometry(lines, start):
    header = lines[start].split()
    require(len(header) == 4 and header[:2] == ['*', 'xyz'], 'invalid xyz header')
    require(all(re.fullmatch(r'[+-]?\d+', v) for v in header[2:]), 'invalid state')
    charge, spin = map(int, header[2:])
    require(-2**31 <= charge < 2**31 and 0 < spin < 2**32, 'state range')
    atoms = []
    for i in range(start + 1, len(lines)):
        line = lines[i].strip()
        if line == '*':
            require(atoms, 'empty geometry')
            return dict(kind='xyz', charge=charge, multiplicity=spin, atoms=atoms), i + 1
        fields = line.split()
        require(len(fields) >= 4 and re.fullmatch(r'[A-Z][a-z]?', fields[0]), 'invalid atom row')
        values, frozen = [], []
        for token in fields[1:4]:
            frozen.append(token.endswith('$'))
            number = token[:-1] if frozen[-1] else token
            require(re.fullmatch(r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?', number), 'invalid number')
            value = float(number)
            require(math.isfinite(value), 'nonfinite coordinate')
            values.append(value)
        atoms.append(dict(element=fields[0], x=values[0], y=values[1], z=values[2], frozen=frozen,
                          extra=' '.join(fields[4:]) if len(fields) > 4 else None))
    raise ValueError('unterminated geometry')


def outcome(case, result, common_only=False):
    require(type(result) is dict and type(result.get('ok')) is bool, 'invalid envelope')
    if not result['ok']:
        require('output' not in result and isinstance(result.get('error'), str), 'partial refusal')
        return dict(status='refused')
    try:
        got = observe(result['output'])
    except (ValueError, KeyError, TypeError) as error:
        return dict(status='unreadable', error=str(error))
    want = case['fields']
    atoms, expected = got['coords']['atoms'], want['coords']['atoms']
    metadata = (len(atoms) == len(expected) and got['coords']['charge'] == want['coords']['charge']
                and got['coords']['multiplicity'] == want['coords']['multiplicity'])
    unit = 'bohr' if any(k.lower() == 'bohrs' for k in want['keywords']) else 'angstrom'
    metadata &= got['units'] == unit
    for actual, original in zip(atoms, expected):
        metadata &= all(actual[k] == original[k] for k in ('element', 'frozen', 'extra'))
    if not common_only:
        metadata &= all(got[k] == want[k] for k in ('comments', 'keywords', 'blocks'))
    actual_bits = [[struct.pack('>d', a[k]).hex() for k in ('x', 'y', 'z')] for a in atoms]
    numeric = actual_bits == case['coordinate_bits']
    maximum = max((abs(a[k] - b[k]) for a, b in zip(atoms, expected) for k in ('x', 'y', 'z')), default=0)
    return dict(status='fields_match' if metadata and numeric else 'changed', metadata_match=metadata,
                numeric_match=numeric, max_coordinate_error=maximum)
