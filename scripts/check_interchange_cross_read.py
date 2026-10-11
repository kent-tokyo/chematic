"""Read the archived common-lane CJSON candidates with the opposite engine.

Run as: python -m scripts.check_interchange_cross_read --paired-report ...
--output NEW.json. Also measures Open Babel reading its own CJSON as a control.
Use --include-chematic-mol2 to additionally test CheMatic reading both engines'
CJSON candidates into MOL2; the original nine-row mode stays available.
Every read attempt, including failure, has a terminal row and archived pipes.
"""
import argparse
from collections import Counter
from pathlib import Path

from scripts.check_interchange_preservation import (
    bounded_run, evaluate, identity, publish_report, read_mol2, require, strict_json,
)


def verified_bytes(item):
    data = Path(item['path']).read_bytes()
    actual = identity(item['path'], data)
    require(actual['sha256'] == item['sha256'] and actual['size_bytes'] == item['size_bytes'],
            'archived artifact differs from selected report')
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--paired-report', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--include-chematic-mol2', action='store_true')
    args = parser.parse_args()
    try:
        require(not args.output.exists(), 'output already exists')
        pair_bytes = args.paired_report.read_bytes()
        pair = strict_json(pair_bytes)
        require(pair['schema'] == 'chematic.paired-interchange-development.v1', 'unsupported paired report')
        protocol = strict_json(verified_bytes(pair['protocol']))
        engines = pair['engines']
        for engine in engines.values():
            verified_bytes(engine)
        for runtime in engines['openbabel']['runtime_artifacts']:
            verified_bytes(runtime)
        cli, obabel = engines['chematic']['path'], engines['openbabel']['path']
        artifacts = args.output.with_suffix(args.output.suffix + '.artifacts')
        artifacts.mkdir()
        rows = []
        for case in protocol['cases']:
            if case['profile'] != 'common':
                continue
            source_rows = {r['engine']: r for r in pair['rows']
                           if r['case'] == case['id'] and r['route'] == 'mol2-cjson'}
            native = strict_json(verified_bytes(source_rows['chematic']['stdout']))['output'].encode('utf8')
            foreign = verified_bytes(source_rows['openbabel']['stdout'])
            directions = [
                ('chematic_to_openbabel', native, [obabel, '-icjson', '-omol2']),
                ('openbabel_self_control', foreign, [obabel, '-icjson', '-omol2']),
                ('openbabel_to_chematic_source', foreign, [cli, 'convert-with-report', '--input-format', 'cjson', '--output-format', 'cjson', '--policy', 'strict']),
            ]
            if args.include_chematic_mol2:
                reverse = [cli, 'convert-with-report', '--input-format', 'cjson', '--output-format', 'mol2', '--policy', 'allow_loss']
                directions.extend([
                    ('chematic_self_control_mol2', native, reverse),
                    ('openbabel_to_chematic_mol2', foreign, reverse),
                ])
            for direction, data, command in directions:
                run = bounded_run(command, data, timeout=protocol['limits']['timeout_seconds'],
                                  max_output=protocol['limits']['max_output_bytes'])
                row = {'case': case['id'], 'direction': direction, 'command': command,
                       'resource_status': run['resource_status'], 'exit_code': run['exit_code']}
                for stream, content in [('input.cjson', data), ('stdout', run['stdout']), ('stderr', run['stderr'])]:
                    path = artifacts / f"{case['id']}-{direction}.{stream}"
                    path.write_bytes(content)
                    row[stream] = identity(path, content)
                try:
                    require(run['resource_status'] == 'completed' and run['exit_code'] == 0, 'read process failed or exceeded budget')
                    if direction == 'openbabel_to_chematic_source':
                        report = strict_json(run['stdout'])
                        equal = report['status'] == 'completed' and report['output'].encode('utf8') == data
                        row['outcome'] = 'exact_source_retained' if equal else 'source_retention_failed'
                        row['opaque_semantics_interpreted'] = False
                    else:
                        candidate = run['stdout'].decode('utf8')
                        if direction in ('chematic_self_control_mol2', 'openbabel_to_chematic_mol2'):
                            report = strict_json(run['stdout'])
                            require(report['status'] in ('completed', 'completed_with_loss') and isinstance(report['output'], str), 'missing MOL2 candidate')
                            row['conversion_status'] = report['status']
                            row['candidate_graph_verified'] = report['candidate_graph_verified']
                            candidate = report['output']
                        check = evaluate(case['expected'], read_mol2(candidate),
                                         coordinate_tolerance=protocol['coordinate_abs_tolerance']['mol2'],
                                         charge_tolerance=protocol['partial_charge_abs_tolerance']['mol2'],
                                         mapping_budget=protocol['limits']['mapping_search_nodes'])
                        row.update(check)
                        statuses = check['fields'].values()
                        row['outcome'] = 'output_with_loss' if 'lost' in statuses else 'output_unverified' if 'unverified' in statuses else 'preserved_required_fields'
                except (ValueError, KeyError, TypeError, IndexError) as error:
                    row['outcome'] = run['resource_status'] if run['resource_status'] != 'completed' else 'parse_write_failure'
                    row['reason'] = str(error)
                rows.append(row)
        for engine in engines.values():
            verified_bytes(engine)
        for runtime in engines['openbabel']['runtime_artifacts']:
            verified_bytes(runtime)
        expected_count = sum(c['profile'] == 'common' for c in protocol['cases']) * (5 if args.include_chematic_mol2 else 3)
        require(len(rows) == expected_count, 'incomplete cross-read accounting')
        result = {'schema': 'chematic.interchange-cross-read-development.v1',
                  'scope': 'Exposed common-lane cross reads plus Open Babel self-read control; no holdout or superiority claim.',
                  'paired_report': identity(args.paired_report, pair_bytes), 'runner': identity(__file__),
                  'engines': engines, 'protocol': pair['protocol'],
                  'expected_terminal_count': expected_count, 'terminal_count': len(rows),
                  'include_chematic_mol2': args.include_chematic_mol2,
                  'outcomes': dict(Counter(r['outcome'] for r in rows)), 'rows': rows,
                  'limitations': ['Open Babel uses default CJSON-read/MOL2-write options; names, typing and charge generation can differ from its CJSON writer.', 'Exact source retention does not interpret opaque semantics. Optional CheMatic MOL2 projection checks the bounded graph and typed name/partial-charge fields; generated atom types/names/IDs and unsupported metadata remain losses.', 'The same declared per-field tolerances and tiny independent reader are used; this is development evidence only.']}
        publish_report(args.output, result)
        print(result['outcomes'])
    except (ValueError, KeyError, TypeError, OSError) as error:
        parser.exit(2, f'error: {error}\n')


if __name__ == '__main__':
    main()
