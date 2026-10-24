"""Independent observer for a bounded numeric CIF/STAR development profile.

Only one data block, scalar items and flat loops; no semicolon text fields,
quoted reserved words, uncertainty suffixes, dictionary/schema validation,
fractional-coordinate reconstruction or general CIF conformance claim.
Missing numeric columns are losses, never parser defaults.
"""
from __future__ import annotations
import math
import shlex
import struct


def require(condition, message):
    if not condition:
        raise ValueError(message)


def tables(text):
    require(not any(line.startswith(';') for line in text.splitlines()), 'text field outside observer profile')
    lex = shlex.shlex(text, posix=True)
    lex.whitespace_split = True
    lex.commenters = '#'
    tokens = list(lex)
    scalars, loops, blocks, i = {}, [], 0, 0
    while i < len(tokens):
        token = tokens[i].lower()
        if token.startswith('data_'):
            blocks += 1
            require(blocks == 1 and i == 0, 'one data block required')
            i += 1
        elif token == 'loop_':
            i += 1
            headers = []
            while i < len(tokens) and tokens[i].startswith('_'):
                headers.append(tokens[i].lower())
                i += 1
            require(headers and len(set(headers)) == len(headers), 'missing or duplicate loop columns')
            values = []
            while i < len(tokens) and not tokens[i].startswith('_') and tokens[i].lower() != 'loop_':
                require(not tokens[i].lower().startswith(('data_', 'save_', 'global_', 'stop_')), 'unexpected block boundary')
                values.append(tokens[i])
                i += 1
            require(len(values) % len(headers) == 0, 'incomplete loop row')
            loops.append((headers, [dict(zip(headers, values[j:j + len(headers)], strict=True))
                                    for j in range(0, len(values), len(headers))]))
        elif token.startswith('_'):
            require(i + 1 < len(tokens) and token not in scalars, 'missing or duplicate scalar')
            scalars[token] = tokens[i + 1]
            i += 2
        else:
            raise ValueError('unexpected token')
    require(blocks == 1, 'missing data block')
    atom_loops = [rows for headers, rows in loops if any(h.startswith('_atom_site.') for h in headers)]
    require(len(atom_loops) == 1 and atom_loops[0], 'one nonempty atom loop required')
    return scalars, atom_loops[0]


def bits(values):
    result = []
    for value in values:
        require(type(value) in (str, int, float), 'not numeric')
        number = float(value)
        require(math.isfinite(number), 'not finite')
        result.append(struct.pack('>d', number).hex())
    return result


def observe(expected, output):
    fields = dict(atom_identity=False, cartesian=False, occupancy=False, b_iso=False,
                  cell=None if expected['cell'] is None else False)
    result = dict(outcome='no_output', fields=fields, preserved=False)
    if not isinstance(output, str) or not output.strip():
        return result
    result['outcome'] = 'invalid_output'
    try:
        scalars, atoms = tables(output)
    except ValueError as error:
        result['observer_error'] = str(error)
        return result
    result['outcome'] = 'cif_output'
    wanted = expected['atoms']
    fields['atom_identity'] = len(atoms) == len(wanted) and all(
        a.get('_atom_site.id') == str(b['serial']) and
        a.get('_atom_site.type_symbol', '').upper() == b['element'].upper()
        for a, b in zip(atoms, wanted))
    channels = dict(cartesian=(('cartn_x', 'cartn_y', 'cartn_z'), ('x', 'y', 'z')),
                    occupancy=(('occupancy',), ('occupancy',)),
                    b_iso=(('b_iso_or_equiv',), ('b_iso',)))
    missing, actual_bits = [], {}
    for channel, (tags, keys) in channels.items():
        tags = ['_atom_site.' + tag for tag in tags]
        absent = [tag for tag in tags if any(tag not in atom or atom[tag] in ('.', '?') for atom in atoms)]
        missing.extend(absent)
        try:
            actual = bits(a[tag] for a in atoms for tag in tags)
            actual_bits[channel] = actual
            fields[channel] = fields['atom_identity'] and not absent and actual == bits(a[key] for a in wanted for key in keys)
        except (ValueError, KeyError, OverflowError):
            actual_bits[channel] = None
    if expected['cell'] is not None:
        keys = ('a', 'b', 'c', 'alpha', 'beta', 'gamma')
        tags = ['_cell.' + ('length_' + key if key in 'abc' else 'angle_' + key) for key in keys]
        missing.extend(tag for tag in tags if tag not in scalars or scalars[tag] in ('.', '?'))
        try:
            actual_bits['cell'] = bits(scalars[tag] for tag in tags)
            fields['cell'] = actual_bits['cell'] == bits(expected['cell'][key] for key in keys)
        except (ValueError, KeyError, OverflowError):
            actual_bits['cell'] = None
    result.update(actual_bits=actual_bits, missing_fields=sorted(set(missing)),
                  preserved=all(value for value in fields.values() if value is not None))
    return result
