#!/usr/bin/env python3
"""Paired, bounded MOL2 preservation measurements with explicit fixture oracles.

POSIX CLI development lane only. No engine is the oracle, no timing/ranking or
independent holdout claim. Output losses and safe refusals remain separate.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
import os
import platform
from pathlib import Path
import re
import selectors
import signal
import struct
import subprocess
import tempfile
import sys
import time


class ContractError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise ContractError(message)


def strict_json(data):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, f"duplicate JSON key: {key}")
            result[key] = value
        return result

    def invalid(value):
        raise ContractError(f"nonfinite JSON literal: {value}")

    return json.loads(data, object_pairs_hook=unique, parse_constant=invalid)


def identity(path, data=None):
    path = Path(path)
    if data is None:
        data = path.read_bytes()
    return {"path": str(path.resolve()), "sha256": hashlib.sha256(data).hexdigest(), "size_bytes": len(data)}


def bounded_run(command, data=b"", *, timeout=10, max_output=1048576):
    """Bound both pipes together while they are being read, including stdin time.

    Kill the process group on timeout/overflow so inherited child pipes cannot
    keep the run alive. communicate() would accumulate unlimited output first.
    """
    require(os.name == "posix", "bounded CLI runner currently requires POSIX")
    require(timeout > 0 and max_output > 0, "invalid process budget")
    deadline = time.monotonic() + timeout
    proc = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, start_new_session=True)
    output = {"stdout": bytearray(), "stderr": bytearray()}
    status = "completed"
    finished = False
    position = 0
    selector = selectors.DefaultSelector()
    pipes = [proc.stdin, proc.stdout, proc.stderr]
    try:
        for pipe in pipes:
            os.set_blocking(pipe.fileno(), False)
        for name, pipe in [("stdout", proc.stdout), ("stderr", proc.stderr)]:
            selector.register(pipe, selectors.EVENT_READ, name)
        if data:
            selector.register(proc.stdin, selectors.EVENT_WRITE, "stdin")
        else:
            proc.stdin.close()
        while selector.get_map():
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                status = "timeout"
                break
            for key, _ in selector.select(remaining):
                if key.data == "stdin":
                    try:
                        position += os.write(key.fd, data[position:position + 4096])
                    except BlockingIOError:
                        continue
                    except BrokenPipeError:
                        position = len(data)
                    if position == len(data):
                        selector.unregister(key.fileobj)
                        key.fileobj.close()
                else:
                    try:
                        chunk = os.read(key.fd, 65536)
                    except BlockingIOError:
                        continue
                    if not chunk:
                        selector.unregister(key.fileobj)
                        key.fileobj.close()
                        continue
                    available = max_output - sum(map(len, output.values()))
                    output[key.data].extend(chunk[:available])
                    if len(chunk) > available:
                        status = "output_limit"
                        break
            if status != "completed":
                break
        if status == "completed":
            try:
                proc.wait(timeout=max(0, deadline - time.monotonic()))
                finished = True
            except subprocess.TimeoutExpired:
                status = "timeout"
    finally:
        if not finished:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        proc.wait()
        selector.close()
        for pipe in pipes:
            pipe.close()
    return {"resource_status": status, "exit_code": proc.returncode,
            **{name: bytes(value) for name, value in output.items()}}


def finite(value):
    require(type(value) in (int, float) and math.isfinite(value), "expected finite number")
    return float(value)


def integer(value):
    require(type(value) is int, "expected integer")
    return value


ELEMENTS = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "P": 15,
            "S": 16, "Cl": 17, "Br": 35, "I": 53}
ORDERS = {"1": "single", "2": "double", "3": "triple", "am": "single", "ar": "aromatic"}


def read_mol2(text, max_atoms=64):
    sections = {}
    active = None
    for line in text.splitlines():
        if line.startswith("@<TRIPOS>"):
            active = line[len("@<TRIPOS>"):].strip()
            require(active not in sections, "repeated MOL2 section or record")
            sections[active] = []
        elif active:
            sections[active].append(line)
    header = sections.get("MOLECULE", [])
    require(len(header) >= 4, "missing MOL2 header")
    counts = header[1].split()
    require(len(counts) >= 2, "missing MOL2 counts")
    n, m = map(int, counts[:2])
    require(0 <= n <= max_atoms and 0 <= m <= max_atoms * (max_atoms - 1) // 2, "MOL2 counts outside budget")
    atoms, ids = [], {}
    for line in sections.get("ATOM", []):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        require(len(fields) >= 6, "short MOL2 atom")
        source_id = int(fields[0])
        require(source_id > 0 and source_id not in ids, "invalid/duplicate MOL2 atom ID")
        symbol = fields[5].split(".")[0]
        require(symbol in ELEMENTS, "atom type outside fixture reader profile")
        partial = finite(float(fields[8])) if len(fields) > 8 else None
        ids[source_id] = len(atoms)
        atoms.append({"source_id": source_id, "element": ELEMENTS[symbol],
                      "formal_charge": 1 if fields[5] == "N.4" else 0,
                      "partial_charge": partial if header[3].strip() != "NO_CHARGES" else None,
                      "coordinates": [finite(float(v)) for v in fields[2:5]],
                      "name": fields[1], "atom_type": fields[5]})
    bonds, edges, bond_ids = [], set(), set()
    for line in sections.get("BOND", []):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        require(len(fields) >= 4, "short MOL2 bond")
        bond_id, a, b = map(int, fields[:3])
        require(bond_id > 0 and bond_id not in bond_ids, "invalid/duplicate MOL2 bond ID")
        require(a in ids and b in ids and a != b, "invalid MOL2 bond reference")
        edge = tuple(sorted((ids[a], ids[b])))
        require(edge not in edges and fields[3] in ORDERS, "duplicate/unsupported MOL2 bond")
        bond_ids.add(bond_id)
        edges.add(edge)
        bonds.append({"atoms": list(edge), "order": ORDERS[fields[3]]})
    require((len(atoms), len(bonds)) == (n, m), "MOL2 declared/actual counts differ")
    return {"atoms": atoms, "bonds": bonds, "name": header[0],
            "charge_method": header[3].strip(), "text": text}


def read_cjson(text, max_atoms=64):
    value = strict_json(text)
    require(isinstance(value, dict), "CJSON root must be object")
    numbers = value["atoms"]["elements"]["number"]
    require(isinstance(numbers, list) and len(numbers) <= max_atoms, "CJSON atom budget")
    n = len(numbers)
    atoms_value = value["atoms"]
    charges = atoms_value.get("formalCharges", [0] * n)
    require(isinstance(charges, list) and len(charges) == n, "CJSON charge length")
    coords = atoms_value.get("coords", {}).get("3d")
    require(coords is None or isinstance(coords, list) and len(coords) == 3 * n, "CJSON coordinate length")
    partials = value.get("partialCharges", {})
    require(isinstance(partials, dict) and len(partials) <= 1, "CJSON partial-charge methods outside profile")
    method, partial = next(iter(partials.items()), (None, None))
    require(partial is None or isinstance(partial, list) and len(partial) == n, "CJSON partial-charge length")
    atoms = [{"element": integer(z), "formal_charge": integer(charges[i]),
              "partial_charge": finite(partial[i]) if partial is not None else None,
              "coordinates": [finite(v) for v in coords[3*i:3*i+3]] if coords is not None else None,
              "source_id": None, "name": None, "atom_type": None}
             for i, z in enumerate(numbers)]
    b = value.get("bonds", {})
    indices = b.get("connections", {}).get("index", [])
    require(isinstance(indices, list) and len(indices) % 2 == 0, "CJSON index pairs")
    orders = b.get("order", [1] * (len(indices) // 2))
    require(isinstance(orders, list) and len(orders) == len(indices) // 2, "CJSON bond order length")
    known = {1: "single", 2: "double", 3: "triple", 1.5: "aromatic", 4: "aromatic"}
    bonds, edges = [], set()
    for i in range(len(orders)):
        a, z = map(integer, indices[2*i:2*i+2])
        edge = tuple(sorted((a, z)))
        require(0 <= a < n and 0 <= z < n and a != z and edge not in edges, "invalid CJSON edge")
        require(type(orders[i]) in (int, float) and orders[i] in known, "unsupported CJSON order")
        edges.add(edge)
        bonds.append({"atoms": list(edge), "order": known[orders[i]]})
    return {"atoms": atoms, "bonds": bonds, "name": value.get("properties", {}).get("name", value.get("name")),
            "charge_method": method, "text": text}


def correspondence(expected, candidate, budget=10000):
    """Unique element + unlabelled-topology bijection; never fit charges/coords/orders.

    Symmetry and search exhaustion are explicit unknowns, not arbitrary matches.
    """
    atoms, written = expected["atoms"], candidate["atoms"]
    if len(atoms) != len(written) or len(expected["bonds"]) != len(candidate["bonds"]):
        return None, "count_mismatch"
    n = len(atoms)
    def adjacency(record):
        result = [set() for _ in range(n)]
        for bond in record["bonds"]:
            a, b = bond["atoms"]
            result[a].add(b)
            result[b].add(a)
        return result
    left, right = adjacency(expected), adjacency(candidate)
    options = [[j for j in range(n) if a["element"] == written[j]["element"] and len(left[i]) == len(right[j])]
               for i, a in enumerate(atoms)]
    order = sorted(range(n), key=lambda i: len(options[i]))
    found, mapping, used = [], {}, set()
    nodes, exhausted = 0, False
    def visit(depth):
        nonlocal nodes, exhausted
        nodes += 1
        if nodes > budget:
            exhausted = True
            return
        if depth == n:
            found.append([mapping[i] for i in range(n)])
            return
        a = order[depth]
        for b in options[a]:
            if b in used or any((i in left[a]) != (j in right[b]) for i, j in mapping.items()):
                continue
            mapping[a] = b
            used.add(b)
            visit(depth + 1)
            used.remove(b)
            del mapping[a]
            if exhausted or len(found) >= 2:
                break
    visit(0)
    if exhausted:
        return None, "search_limit"
    if len(found) == 1:
        return found[0], "unique"
    return None, "ambiguous" if found else "not_isomorphic"


def evaluate(expected, candidate, *, coordinate_tolerance, charge_tolerance, mapping_budget=10000):
    mapping, reason = correspondence(expected, candidate, mapping_budget)
    fields = {}
    def check(field, condition):
        fields[field] = "preserved" if condition else "lost"
    check("atom_count", len(expected["atoms"]) == len(candidate["atoms"]))
    check("bond_count", len(expected["bonds"]) == len(candidate["bonds"]))
    check("elements", Counter(a["element"] for a in expected["atoms"]) == Counter(a["element"] for a in candidate["atoms"]))
    fields["connectivity"] = "preserved" if reason in ("unique", "ambiguous") else "unverified" if reason == "search_limit" else "lost"
    for field, key in [("formal_charges", "formal_charge"), ("coordinates", "coordinates"),
                       ("partial_charges", "partial_charge"), ("atom_names", "name"),
                       ("atom_types", "atom_type"), ("source_atom_ids", "source_id")]:
        if not expected["atoms"]:
            fields[field] = "not_applicable"
        elif not candidate["atoms"] or any(a[key] is None for a in candidate["atoms"]):
            fields[field] = "lost"
        elif mapping is None:
            fields[field] = "unverified"
        else:
            pairs = [(a[key], candidate["atoms"][mapping[i]][key]) for i, a in enumerate(expected["atoms"])]
            if key == "coordinates":
                ok = all(abs(x-y) <= coordinate_tolerance for a, b in pairs for x, y in zip(a, b))
            elif key == "partial_charge":
                ok = all(abs(a-b) <= charge_tolerance for a, b in pairs)
            else:
                ok = all(a == b for a, b in pairs)
            check(field, ok)
    if mapping is None:
        fields["bond_orders"] = "unverified"
    else:
        orders = {tuple(sorted(b["atoms"])): b["order"] for b in candidate["bonds"]}
        check("bond_orders", all(orders[tuple(sorted(mapping[i] for i in b["atoms"]))] == b["order"] for b in expected["bonds"]))
    check("molecule_name", expected["name"] == candidate["name"])
    check("charge_method", expected["charge_method"] == candidate["charge_method"])
    markers = expected.get("opaque_markers", [])
    fields["opaque_marker_retention"] = "not_applicable" if not markers else "preserved" if all(m in candidate["text"] for m in markers) else "lost"
    exact = None
    if mapping is not None and fields["coordinates"] == "preserved":
        exact = all(struct.pack('>d', x) == struct.pack('>d', y)
                    for i, a in enumerate(expected["atoms"])
                    for x, y in zip(a["coordinates"], candidate["atoms"][mapping[i]]["coordinates"]))
    return {"fields": fields, "correspondence": reason, "source_to_candidate": mapping,
            "coordinate_bits_equal_diagnostic": exact}


def load_protocol(path, data=None):
    if data is None:
        data = path.read_bytes()
    require(len(data) <= 1048576, "protocol exceeds byte budget")
    protocol = strict_json(data)
    require(protocol["schema"] == "chematic.interchange-development-protocol.v1", "unsupported protocol")
    limits = protocol["limits"]
    for key in ["max_input_bytes", "max_output_bytes", "max_atoms", "mapping_search_nodes"]:
        require(type(limits[key]) is int and limits[key] > 0, f"invalid limit: {key}")
    require(limits["max_input_bytes"] <= 1048576 and limits["max_output_bytes"] <= 16777216 and limits["max_atoms"] <= 64 and limits["mapping_search_nodes"] <= 100000, "protocol exceeds hard budget")
    require(0 < finite(limits["timeout_seconds"]) <= 60, "timeout outside budget")
    for kind in ["coordinate_abs_tolerance", "partial_charge_abs_tolerance"]:
        require(all(0 <= finite(protocol[kind][target]) <= 0.001 for target in ["mol2", "cjson"]), "invalid numeric tolerance")
    require(0 < len(protocol["cases"]) <= 128, "corpus size outside budget")
    cases, seen = [], set()
    for case in protocol["cases"]:
        require(re.fullmatch(r'[a-z0-9-]+', case["id"]) and case["id"] not in seen, "invalid/duplicate case ID")
        require(case["profile"] in ("common", "stress", "negative") and type(case["valid_in_declared_profile"]) is bool, "invalid fixture profile")
        require((case["profile"] == "negative") == (not case["valid_in_declared_profile"]), "fixture validity/profile conflict")
        seen.add(case["id"])
        fixture = (path.parent / case["file"]).resolve()
        require(fixture.is_relative_to(path.parent.resolve()), "fixture escapes protocol directory")
        source = fixture.read_bytes()
        require(len(source) <= limits["max_input_bytes"], "fixture exceeds input budget")
        try:
            parsed = read_mol2(source.decode('utf8'), limits["max_atoms"])
        except (ValueError, KeyError, IndexError) as error:
            require(not case["valid_in_declared_profile"], f"valid fixture cannot be read: {case['id']}: {error}")
        else:
            require(case["valid_in_declared_profile"], f"negative fixture accepted by independent reader: {case['id']}")
            check = evaluate(case["expected"], parsed, coordinate_tolerance=0, charge_tolerance=0,
                             mapping_budget=limits["mapping_search_nodes"])
            require(all(v in ("preserved", "not_applicable") for v in check["fields"].values()), f"fixture disagrees with explicit oracle: {case['id']}: {check}")
        cases.append((case, source, identity(fixture, source)))
    require(cases, "empty corpus")
    require(protocol["routes"] == [{"id": "mol2-self", "target": "mol2", "chematic_policy": "strict"}, {"id": "mol2-cjson", "target": "cjson", "chematic_policy": "allow_loss"}], "unsupported protocol routes")
    return protocol, cases, identity(path, data)


def measure(engine, command, source, case, route, protocol, artifacts):
    run = bounded_run(command, source, timeout=protocol["limits"]["timeout_seconds"],
                      max_output=protocol["limits"]["max_output_bytes"])
    row = {"case": case["id"], "profile": case["profile"], "route": route["id"], "engine": engine,
           "command": command, "resource_status": run["resource_status"], "exit_code": run["exit_code"]}
    for stream in ["stdout", "stderr"]:
        path = artifacts / f"{case['id']}-{route['id']}-{engine}.{stream}"
        path.write_bytes(run[stream])
        row[stream] = identity(path, run[stream])
    if run["resource_status"] != "completed":
        row["outcome"] = run["resource_status"]
        return row
    report = None
    text = run["stdout"].decode('utf8', errors='replace')
    try:
        if engine == "chematic":
            if not text.strip():
                row["outcome"] = "safe_rejection" if not case["valid_in_declared_profile"] and run["exit_code"] != 0 else "parse_write_failure"
                return row
            report = strict_json(text)
            require(report["schema"] == "chematic.conversion-report.v1", "unexpected report schema")
            if report["status"] == "refused":
                require(report["output"] is None and run["exit_code"] == 2, "invalid refusal envelope/exit")
                row["outcome"] = "safe_refusal"
                return row
            require(report["status"] in ("completed", "completed_with_loss") and isinstance(report["output"], str), "invalid output envelope")
            text = report["output"]
            row["reported_observations"] = report["observations"]
        require(run["exit_code"] == 0, "process failed")
        candidate = (read_mol2 if route["target"] == "mol2" else read_cjson)(text, protocol["limits"]["max_atoms"])
        if not case["valid_in_declared_profile"]:
            row["outcome"] = "accepted_outside_profile"
            return row
        checks = evaluate(case["expected"], candidate,
                          coordinate_tolerance=protocol["coordinate_abs_tolerance"][route["target"]],
                          charge_tolerance=protocol["partial_charge_abs_tolerance"][route["target"]],
                          mapping_budget=protocol["limits"]["mapping_search_nodes"])
        row.update(checks)
        fields = checks["fields"].values()
        row["outcome"] = "output_with_loss" if "lost" in fields else "output_unverified" if "unverified" in fields else "preserved_required_fields"
        row["machine_loss_report_available"] = report is not None
        if report is not None and report["candidate_graph_verified"]:
            scoped = ["atom_count", "bond_count", "elements", "connectivity", "formal_charges", "bond_orders"]
            row["independent_graph_check_agrees"] = all(checks["fields"][f] == "preserved" for f in scoped)
    except (ValueError, KeyError, IndexError, TypeError) as error:
        row["outcome"] = "safe_rejection" if not case["valid_in_declared_profile"] and not run["stdout"].strip() else "parse_write_failure"
        row["reason"] = str(error)
    return row


FIELD_NAMES = ("atom_count", "bond_count", "elements", "connectivity", "formal_charges",
               "coordinates", "partial_charges", "atom_names", "atom_types", "source_atom_ids",
               "bond_orders", "molecule_name", "charge_method", "opaque_marker_retention")


def summarize(rows):
    cells = defaultdict(list)
    for row in rows:
        cells[(row["profile"], row["route"], row["engine"])].append(row)
    return [{"profile": profile, "route": route, "engine": engine, "input_count": len(items),
             "terminal_count": len(items), "outcomes": dict(Counter(r["outcome"] for r in items)),
             "fields": {field: dict(Counter(r.get("fields", {}).get(field, "no_candidate") for r in items))
                        for field in FIELD_NAMES}}
            for (profile, route, engine), items in sorted(cells.items())]


def publish_report(path, value):
    """Atomically publish without overwriting an existing verification record."""
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, prefix='.' + path.name, delete=False) as stream:
        temp = Path(stream.name)
        try:
            json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
            os.link(temp, path)
        finally:
            temp.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol', type=Path, default=Path(__file__).resolve().parents[2] / 'validation/interchange/protocol.json')
    parser.add_argument('--chematic', type=Path, required=True)
    parser.add_argument('--chematic-evidence', type=Path, required=True)
    parser.add_argument('--obabel', type=Path, required=True)
    parser.add_argument('--expected-obabel-version', required=True)
    parser.add_argument('--expected-obabel-sha256', required=True)
    parser.add_argument('--obabel-runtime', type=Path, action='append', default=[], help='Explicit dylib/plugin artifacts; repeat as needed.')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        require(not args.output.exists(), 'output already exists; choose a new verification record')
        protocol_bytes = args.protocol.read_bytes()
        runner_bytes = Path(__file__).read_bytes()
        protocol, cases, protocol_id = load_protocol(args.protocol, protocol_bytes)
        cli, obabel = args.chematic.resolve(), args.obabel.resolve()
        require(cli.is_file() and obabel.is_file(), 'missing engine executable')
        cli_id, obabel_id = identity(cli), identity(obabel)
        require(obabel_id["sha256"] == args.expected_obabel_sha256, "Open Babel executable differs from pinned SHA-256")
        evidence_bytes = args.chematic_evidence.read_bytes()
        evidence = strict_json(evidence_bytes)
        require(evidence["schema"] == "chematic.development-verification.v1", "unsupported build evidence schema")
        require(any(a['sha256'] == cli_id['sha256'] and a['size_bytes'] == cli_id['size_bytes'] for a in evidence['artifacts']), 'chematic executable does not match selected build evidence')
        runtimes = [identity(path) for path in args.obabel_runtime]
        versions = {}
        metadata = {"protocol": protocol_bytes, "runner": runner_bytes}
        for name, command in [('chematic', [str(cli), '--version']), ('openbabel', [str(obabel), '-V']), ('openbabel_formats', [str(obabel), '-L', 'formats'])]:
            result = bounded_run(command)
            require(result['resource_status'] == 'completed' and result['exit_code'] == 0, f'{name} metadata command failed')
            text = result['stdout'].decode('utf8')
            metadata[name] = result['stdout']
            versions[name] = text.strip()
        require(re.search(r'Open Babel ' + re.escape(args.expected_obabel_version) + r'(?=\s|$)', versions['openbabel']), 'Open Babel version differs from pinned target')
        require(versions['chematic'] == 'chematic ' + evidence['workspace_version'], 'chematic version differs from selected evidence')
        artifacts = args.output.with_suffix(args.output.suffix + '.artifacts')
        artifacts.mkdir(parents=False, exist_ok=False)
        metadata_ids = []
        for name, data in metadata.items():
            path = artifacts / (name + '.txt')
            path.write_bytes(data)
            metadata_ids.append(identity(path, data))
        rows = []
        for case, source, source_id in cases:
            path = artifacts / (case['id'] + '.input.mol2')
            path.write_bytes(source)
            for route in protocol['routes']:
                for engine, command in [
                    ('chematic', [str(cli), 'convert-with-report', '--input-format', 'mol2', '--output-format', route['target'], '--policy', route['chematic_policy']]),
                    ('openbabel', [str(obabel), '-imol2', '-o' + route['target']]),
                ]:
                    row = measure(engine, command, source, case, route, protocol, artifacts)
                    row['input'] = source_id
                    row['archived_input'] = identity(path, source)
                    rows.append(row)
        for artifact in [cli_id, obabel_id, *runtimes]:
            require(identity(artifact['path']) == artifact, 'engine/runtime artifact changed during measurement')
        expected_count = len(cases) * len(protocol['routes']) * 2
        require(len(rows) == expected_count, 'incomplete terminal row accounting')
        result = {
            'schema': 'chematic.paired-interchange-development.v1',
            'created_at_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
            'scope': protocol['scope'], 'protocol': protocol_id,
            'runner': identity(__file__, runner_bytes), 'oracle': protocol['oracle'],
            'platform': platform.platform(), 'python_runtime': sys.version,
            'openbabel_environment': {key: os.environ.get(key) for key in ['BABEL_LIBDIR', 'BABEL_DATADIR']},
            'policy_note': 'CheMatic checked source-retention/allow_loss policies versus Open Babel default CLI; no equivalent Open Babel policy is asserted. No charge generation or geometry generation flags.',
            'engines': {'chematic': {**cli_id, 'version': versions['chematic'], 'lane': 'archived source-built CLI', 'build_evidence': identity(args.chematic_evidence, evidence_bytes)},
                        'openbabel': {**obabel_id, 'version': versions['openbabel'], 'lane': 'installed CLI', 'runtime_artifacts': runtimes}},
            'metadata': metadata_ids, 'limits': protocol['limits'],
            'row_accounting': {'input_count': len(cases), 'route_count': len(protocol['routes']), 'engine_count': 2, 'expected_terminal_count': expected_count, 'terminal_count': len(rows)},
            'cells': summarize(rows), 'rows': rows,
            'limitations': ['Exposed development fixtures, not an independent holdout or superiority gate.', 'Local MOL2 N.4 formal-charge projection only; generic MOL2 charge interpretation is unverified.', 'Tiny independent readers cover only this declared fixture profile.', 'Unlabelled topology plus elements determine correspondence; symmetry/search limits remain unverified.', 'Opaque marker retention checks bytes only, not editable or interpreted semantics.', 'OS/system libraries, environment and full runtime dependency closure are not archived; listed explicit artifacts only.', 'Failures, refusals and negative cases remain in denominators. No comparative speed or statistical superiority claim.'],
        }
        publish_report(args.output, result)
        print(json.dumps({'report': str(args.output), 'terminal_rows': len(rows), 'outcomes': dict(Counter(r['outcome'] for r in rows))}))
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(2, f'error: {error}\n')


if __name__ == '__main__':
    main()
