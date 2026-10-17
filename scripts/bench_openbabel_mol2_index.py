"""Fixed source-only MOL2 follow-up. Set CHEMATIC_MOL2_EVIDENCE to its explicit journal directory."""

from pathlib import Path
import argparse, datetime, gzip, hashlib, itertools, json, os, random, shutil, statistics, subprocess

E = Path(os.environ.get("CHEMATIC_MOL2_EVIDENCE", ".")).resolve()
R = Path(__file__).resolve().parents[1]
MODES = ("parse", "write", "roundtrip")
os.environ.update(
    CHEMATIC_BENCH_WARMUP_NS="250000000",
    CHEMATIC_BENCH_TIMED_NS="100000000",
    PYTHONDONTWRITEBYTECODE="1",
)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def atomic(p, v):
    part = p.with_suffix(p.suffix + ".partial")
    part.write_text(json.dumps(v, indent=2, allow_nan=False) + "\n")
    part.replace(p)


def run(arm, name, mode, count):
    return subprocess.run(
        [str(E / (arm + "-mol2")), str(E / (name + ".mol2s")), mode, "3", str(count)],
        capture_output=True,
        timeout=240,
    )


def timing(arm, name, mode, count):
    p = run(arm, name, mode, count)
    assert p.returncode == 0, (arm, name, mode, p.stderr.decode()[:300])
    row = json.loads(p.stdout)
    reps = row["repeats"]
    assert type(reps) == int and reps >= 3 and row["records"] == count * reps
    assert (
        row["minimum_warmup_ns"] == 250000000
        and row["warmup_elapsed_ns"] >= 250000000
        and row["warmup_repeats"] > 0
    )
    assert row["minimum_timed_ns"] == 100000000 and row["elapsed_ns"] >= 100000000
    return row


def ci(values, seed=0):
    rng = random.Random(seed)
    n = len(values)
    boot = sorted(statistics.median(rng.choices(values, k=n)) for _ in range(10000))
    return {
        "median": statistics.median(values),
        "ci95": [boot[249], boot[9749]],
        "winning_blocks": sum(v > 1 for v in values),
        "blocks": n,
    }


def ratio(blocks, a, b, mode):
    return ci(
        [
            (x[a][mode]["elapsed_ns"] / x[a][mode]["repeats"])
            / (x[b][mode]["elapsed_ns"] / x[b][mode]["repeats"])
            for x in blocks
        ],
        0 if a == "baseline" else 1,
    )


def fields(text):
    sections = {}
    current = None
    for line in text.splitlines():
        t = line.strip()
        if t.startswith("@<TRIPOS>"):
            current = t[9:]
            sections.setdefault(current, [])
            continue
        if current and t and not t.startswith("#"):
            sections[current].append(t)
    atoms = [line.split() for line in sections.get("ATOM", [])]
    ids = {int(a[0]): i for i, a in enumerate(atoms)}
    bonds = [line.split() for line in sections.get("BOND", [])]
    return {
        "name": sections.get("MOLECULE", [""])[0],
        "molecule_type": sections.get("MOLECULE", ["", "", ""])[2:3],
        "charge_type": sections.get("MOLECULE", ["", "", "", ""])[3:4],
        "atom_ids": [int(a[0]) for a in atoms],
        "bond_ids": [int(b[0]) for b in bonds],
        "atom_names": [a[1] for a in atoms],
        "atom_types": [a[5] for a in atoms],
        "coords": [[float(x) for x in a[2:5]] for a in atoms],
        "substructures": [a[6:8] for a in atoms],
        "partial_charges": [float(a[8]) if len(a) > 8 else None for a in atoms],
        "status_bits": [a[9:] for a in atoms],
        "bonds": [[ids[int(b[1])], ids[int(b[2])], b[3], b[4:]] for b in bonds],
        "extra_sections": {
            k: v for k, v in sections.items() if k not in ("MOLECULE", "ATOM", "BOND")
        },
    }


def graph(text):
    from rdkit import Chem

    m = Chem.MolFromMol2Block(text, sanitize=False, removeHs=False)
    if m is None:
        return None
    return {
        "atoms": [
            [a.GetAtomicNum(), a.GetFormalCharge(), a.GetIsotope(), a.GetIsAromatic()]
            for a in m.GetAtoms()
        ],
        "bonds": [
            [
                b.GetBeginAtomIdx(),
                b.GetEndAtomIdx(),
                str(b.GetBondType()),
                b.GetIsAromatic(),
            ]
            for b in m.GetBonds()
        ],
    }


def preflight():
    assert not (E / "preflight.json").exists()
    inputs = json.loads((E / "inputs.json").read_text())
    results = {}
    for name, info in inputs.items():
        raw = {}
        row = {"records": info["records"], "external_errors": []}
        for arm in ("baseline", "candidate"):
            for mode in ("snapshot", "dump"):
                p = run(arm, name, mode, info["records"])
                assert p.returncode == 0, (name, arm, mode, p.stderr.decode())
                raw[arm, mode] = p.stdout
                (E / f"preflight-{name}-{arm}-{mode}.gz").write_bytes(
                    gzip.compress(p.stdout, mtime=0)
                )
        assert raw["baseline", "snapshot"] == raw["candidate", "snapshot"], name
        assert raw["baseline", "dump"] == raw["candidate", "dump"], name
        row["source_snapshots_equal"] = True
        row["source_rewrites_equal"] = True
        row["snapshot_sha256"] = hashlib.sha256(
            raw["candidate", "snapshot"]
        ).hexdigest()
        row["rewrite_sha256"] = hashlib.sha256(raw["candidate", "dump"]).hexdigest()
        ob = run("openbabel", name, "dump", info["records"])
        row["openbabel_read_success"] = ob.returncode == 0
        (E / f"preflight-{name}-openbabel-stderr.txt").write_bytes(ob.stderr)
        if ob.returncode == 0:
            (E / f"preflight-{name}-openbabel-dump.gz").write_bytes(
                gzip.compress(ob.stdout, mtime=0)
            )
            source = (E / (name + ".mol2s")).read_text().split("\n$$$$\n")
            outputs = ob.stdout.decode().split("\n$$$$\n")
            row["external_errors"] = []
            assert len(source) == len(outputs) == info["records"]
            for i, (a, b) in enumerate(zip(source, outputs, strict=True)):
                try:
                    fa, fb = fields(a), fields(b)
                    for k in fa:
                        if k == "coords":
                            if len(fa[k]) != len(fb[k]) or any(
                                abs(x - y) > 1e-4
                                for ca, cb in zip(fa[k], fb[k])
                                for x, y in zip(ca, cb)
                            ):
                                row["external_errors"].append([i, k])
                        elif fa[k] != fb[k]:
                            row["external_errors"].append([i, k])
                    ga, gb = graph(a), graph(b)
                    if ga is None or gb is None or ga != gb:
                        row["external_errors"].append([i, "indexed_graph"])
                except Exception as err:
                    row["external_errors"].append([i, type(err).__name__])
        else:
            row["external_errors"].append(["profile", "OpenBabel_read_failure"])
        # Parse-only external ranking still requires all retained fields to match.
        row["external_eligible"] = ob.returncode == 0 and not row["external_errors"]
        arms = ["baseline", "candidate"]
        row["external_timing"] = False
        row["instrumentation"] = {
            arm: {mode: timing(arm, name, mode, info["records"]) for mode in MODES}
            for arm in arms
        }
        perms = list(itertools.permutations(arms))
        orders = perms * (24 // len(perms))
        random.Random(20261018).shuffle(orders)
        row["orders"] = [list(p) for p in orders]
        results[name] = row
        print(
            "preflight",
            name,
            "eligible",
            row["external_eligible"],
            "errors",
            row["external_errors"][:5],
            flush=True,
        )
    atomic(
        E / "preflight.json",
        {
            "passed": True,
            "records": sum(x["records"] for x in inputs.values()),
            "formal_samples": 0,
            "profiles": results,
            "rdkit_version": __import__("rdkit").__version__,
        },
    )


def integrity():
    seal = json.loads((E / "seal.json").read_text())
    for p, h in seal["files_sha256"].items():
        assert sha(p) == h, p
    for p, h in seal["source_files_sha256"].items():
        assert sha(R / p) == h, p
    return sha(E / "seal.json")


def formal():
    assert (
        not list(E.glob("paired-*.json")) and not (E / "decision.json").exists()
    ), "no timing reruns"
    pre = json.loads((E / "preflight.json").read_text())
    assert pre["passed"]
    assert json.loads((E / "validation.json").read_text())["passed"]
    seal = integrity()
    primary = json.loads((E / "measurement-plan.json").read_text())["source_gate"][
        "primary"
    ]
    summaries = {}
    floor = []
    negative = []
    for name, info in pre["profiles"].items():
        assert integrity() == seal
        blocks = []
        out = E / ("paired-" + name + ".json")
        for i, order in enumerate(info["orders"]):
            assert (
                shutil.disk_usage(E).free > 500 * 1024**2
            ), "capacity guard before any next invocation"
            block = {
                "index": i,
                "order": order,
                "load_average_start": list(os.getloadavg()),
            }
            for arm in order:
                block[arm] = {}
                for mode in MODES:
                    row = timing(arm, name, mode, info["records"])
                    atomic(E / f"journal-{name}-{i:02}-{arm}-{mode}.json", row)
                    block[arm][mode] = row
            blocks.append(block)
            atomic(
                out,
                {
                    "complete": False,
                    "profile": name,
                    "blocks": blocks,
                    "seal_sha256": seal,
                },
            )
            print(name, i + 1, 24, flush=True)
        summary = {}
        for mode in MODES:
            lane = name + "/" + mode
            source = ratio(blocks, "baseline", "candidate", mode)
            if source["median"] < 0.99:
                floor.append(lane)
            if source["ci95"][1] < 1:
                negative.append(lane)
            external = (
                ratio(blocks, "openbabel", "candidate", mode)
                if "openbabel" in blocks[0]
                else None
            )
            summary[mode] = {
                "source": source,
                "openbabel": external,
                "external_eligible": info["external_eligible"],
                "strict_external_win": bool(
                    info["external_eligible"]
                    and external
                    and external["ci95"][0] > 1
                    and external["winning_blocks"] == 24
                ),
            }
            summaries[lane] = summary[mode]
        assert integrity() == seal
        atomic(
            out,
            {
                "complete": True,
                "profile": name,
                "blocks": blocks,
                "summary": summary,
                "seal_sha256": seal,
            },
        )
    passed = (
        all(
            summaries[l]["source"]["ci95"][0] > 1
            and summaries[l]["source"]["winning_blocks"] == 24
            for l in primary
        )
        and not floor
    )
    atomic(
        E / "decision.json",
        {
            "complete": True,
            "source_gate_passed": passed,
            "source_floor_failures": floor,
            "source_negative_intervals": negative,
            "primary": primary,
            "lanes": summaries,
            "universal_superiority": False,
        },
    )
    print("complete source gate", passed, flush=True)


def prepare_inputs(packet_path):
    packet = json.loads(gzip.decompress(packet_path.read_bytes()))
    assert packet["label"] == "index-and-fixed4"
    E.mkdir(parents=True, exist_ok=True)
    for name, info in packet["inputs"].items():
        if name == "ethanol":
            texts = [(R / "benchmarks/fixtures/ethanol.mol2").read_text()]
        elif name == "mol2-stereo-fixtures":
            texts = [(R / path).read_text() for path in info["fixture_paths"]]
        else:
            n = info["atoms"]
            kind = (
                "sparse"
                if name.startswith("sparse")
                else "reversed" if name.startswith("reversed") else "dense"
            )
            ids = [i + 1 for i in range(n)]
            if kind == "sparse":
                ids = [4294967295 - 13 * i for i in range(n)]
            if kind == "reversed":
                ids = list(reversed(ids))
            atoms = "".join(
                f"{id} C{i+1} 0 0 0 C.3 1 LIG1 0\n" for i, id in enumerate(ids)
            )
            bonds = "".join(f"{i+1} {ids[i]} {ids[i+1]} 1\n" for i in range(n - 1))
            texts = [
                f"@<TRIPOS>MOLECULE\nchain-{n}-{kind}\n{n} {n-1} 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n{atoms}@<TRIPOS>BOND\n{bonds}"
            ]
        path = E / (name + ".mol2s")
        path.write_text("\n$$$$\n".join(texts))
        assert sha(path) == info["sha256"], name
    atomic(E / "inputs.json", packet["inputs"])
    atomic(
        E / "measurement-plan.json",
        {
            "source_gate": packet["source_gate"],
            "external_timing": False,
            "fresh_reproduction": True,
            "packet_sha256": sha(packet_path),
        },
    )
    print(
        "prepared all nine exact input hashes; new timings must use a new empty journal directory"
    )


def seal_reproduction(validation_log):
    import importlib.metadata, re

    assert not list(E.glob("paired-*.json")) and not (E / "seal.json").exists()
    test_text = validation_log.read_text()
    rows = re.findall(
        r"test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored", test_text
    )
    assert (
        rows
        and sum(int(x[0]) for x in rows) >= 1448
        and all(x[1] == "0" and x[2] == "0" for x in rows)
    )
    assert "all doctests ran in" in test_text
    atomic(
        E / "validation.json",
        {
            "passed": True,
            "fresh_test_log_sha256": sha(validation_log),
            "test_suites": len(rows),
        },
    )
    paths = subprocess.check_output(["git", "ls-files"], cwd=R, text=True).splitlines()
    source = {
        path: sha(R / path)
        for path in paths
        if (
            path.startswith(("crates/", "tools/"))
            or path in ("Cargo.toml", "Cargo.lock")
        )
        and (R / path).is_file()
    }
    files = {str(p): sha(p) for p in E.iterdir() if p.is_file()}
    files[str(Path(__file__).resolve())] = sha(Path(__file__))
    prefix = Path(
        subprocess.check_output(
            ["pkg-config", "--variable=prefix", "openbabel-3"], text=True
        ).strip()
    )
    for p in prefix.rglob("*"):
        if p.is_file():
            files[str(p)] = sha(p)
    for name in ("rdkit", "numpy"):
        dist = importlib.metadata.distribution(name)
        for relative in dist.files or []:
            p = Path(dist.locate_file(relative)).resolve()
            if (
                p.is_file()
                and p.suffix not in (".pyc", ".pyo")
                and "__pycache__" not in p.parts
            ):
                files[str(p)] = sha(p)
    atomic(
        E / "seal.json",
        {
            "source_files_sha256": source,
            "files_sha256": files,
            "scope": "reproduction source/binaries/inputs and Open Babel, RDKit, NumPy package contents; system/transitive libraries excluded",
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--prepare-inputs", type=Path)
    parser.add_argument("--seal-after-tests", type=Path)
    args = parser.parse_args()
    if "CHEMATIC_MOL2_EVIDENCE" not in os.environ:
        parser.error("set CHEMATIC_MOL2_EVIDENCE to a new explicit journal directory")
    if args.prepare_inputs:
        prepare_inputs(args.prepare_inputs)
    elif args.seal_after_tests:
        seal_reproduction(args.seal_after_tests)
    elif args.preflight:
        preflight()
    else:
        formal()
