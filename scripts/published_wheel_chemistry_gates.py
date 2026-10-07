#!/usr/bin/env python3
"""Run the #734/#754 chemistry gates on one installed published chematic wheel.

The gates are the ones the v1.0.36 Linux record ran on the PyPI CPython 3.9
wheel; this driver runs them on any platform so the macOS and Windows wheels
get the same numbers (``.github/workflows/published-wheel-chemistry-gates.yml``).

* 83 reaction fixtures, checked provenance (graph, atom origin, template map)
* SMARTS opt-in profile on the pinned 310,000-cell archived oracle
* accurate CIP, hybridization, MMFF94 atom types and MOL-writer stereo on the
  exposed 10k and ChEMBL 5k corpora
* optionally the BioTransformer rule corpus (974 rules x 400 reactants), when
  ``--rules-dir`` holds the three pinned rule tables (they are not vendored:
  part of them is CC BY-NC-SA)

chematic runs under ``--chematic-python`` (the interpreter the wheel is
installed in) and RDKit 2026.03.6 under ``--rdkit-python``; they may be the
same interpreter. Every count is compared with ``--expected`` (the published
Linux numbers) and the driver exits 1 on any difference, so a platform that
types, labels or matches differently fails instead of being averaged in.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
CORPORA = {
    "exposed_10k": ROOT / "validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi",
    "chembl_5k": ROOT / "scripts/chembl_accuracy_corpus_4999.smi",
}
REACTANTS_SHA256 = "ac68df7bc3133fa20d39bdf3384a38df0edd09e28032fdc2015322bf9d41119f"
RULE_FILES = {
    "database_metabolicReactions.json": "9bbfa5b8c0fd3d06539381714ca95fc5ff73fc3da821ce96d3f7bfb58915c04d",
    "database_ENVMICRO_metabolicReactions.json": "f586471d0f3674cda4577d95af73c19f5c779b9267fd663484be34ece3521836",
    "database_standardizationReactions.json": "749c16a0a0471b18a00cebda01b460839e5d7e10e5e7a8e7d7447f20f51804a8",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def finished(output: Path | None) -> bool:
    """Whether a step's JSON output is already complete (``--resume``)."""
    if output is None or not output.is_file():
        return False
    try:
        json.loads(output.read_text(encoding="utf-8"))
    except ValueError:
        return False
    return True


def run(cmd: list[str], log: Path, output: Path | None = None, resume: bool = False) -> None:
    if resume and finished(output):
        with log.open("a", encoding="utf-8") as out:
            out.write(f"# resume: {output} is complete, skipping $ {' '.join(cmd)}\n")
        return
    t0 = time.time()
    with log.open("a", encoding="utf-8") as out:
        out.write("$ " + " ".join(cmd) + "\n")
        out.flush()
        proc = subprocess.run(cmd, cwd=ROOT, stdout=out, stderr=subprocess.STDOUT)
        out.write(f"# exit {proc.returncode} after {time.time() - t0:.0f}s\n")
    if proc.returncode != 0:
        raise SystemExit(f"step failed ({proc.returncode}): {' '.join(cmd)}; see {log}")


def module_root(python: str) -> Path:
    out = subprocess.run(
        [python, "-c", "import sys; print(sys.prefix)"], check=True, capture_output=True, text=True
    )
    return Path(out.stdout.strip())


def reactants(out: Path) -> Path:
    """Every 25th ChEMBL 5k row then every 50th exposed-10k row (the record's 400)."""
    rows: list[str] = []
    for path, every in ((CORPORA["chembl_5k"], 25), (CORPORA["exposed_10k"], 50)):
        lines = [l.split()[0] for l in path.read_text(encoding="utf-8").splitlines() if l.strip()]
        rows.extend(lines[::every])
    out.write_bytes("".join(r + "\n" for r in rows).encode("utf-8"))
    if sha256(out) != REACTANTS_SHA256:
        raise SystemExit(f"reactant set {out} does not hash to the record's {REACTANTS_SHA256}")
    return out


def merge_corpus_shards(shards: list[tuple[Path, Path]], output: Path, rows: Path) -> None:
    """Add up the counts of reactant shards of the BioTransformer corpus."""
    from collections import Counter

    counts: Counter = Counter()
    per_rule: dict = {}
    merged: dict = {}
    reactant_count = 0
    elapsed = 0.0
    with rows.open("w", encoding="utf-8") as rows_out:
        for summary_path, rows_path in shards:
            part = json.loads(summary_path.read_text(encoding="utf-8"))
            if not merged:
                merged = {k: v for k, v in part.items() if k not in {"counts", "per_rule"}}
            # rules_run and the per-rule skip counts repeat in every shard.
            for key, value in part["counts"].items():
                if ":" in key:
                    counts[key] += value
                else:
                    counts[key] = value
            for rule, stats in part["per_rule"].items():
                if not all(type(v) is int for v in stats.values()):
                    per_rule[rule] = stats
                    continue
                into = per_rule.setdefault(rule, {})
                for key, value in stats.items():
                    into[key] = into.get(key, 0) + value
            reactant_count += part["reactants"]["count"]
            elapsed += part["elapsed_seconds"]
            rows_out.write(rows_path.read_text(encoding="utf-8"))
    merged["reactants"] = {**merged["reactants"], "count": reactant_count}
    merged["elapsed_seconds"] = round(elapsed, 1)
    merged["shards"] = len(shards)
    merged["counts"] = dict(sorted(counts.items()))
    merged["per_rule"] = per_rule
    output.write_text(json.dumps(merged, indent=1) + "\n", encoding="utf-8")


def summarize(out_dir: Path, ran_corpus: bool) -> dict:
    s: dict = {}
    g = json.loads((out_dir / "reaction-checked-provenance-83.json").read_text())
    s["reaction_83"] = g["accounting"]["outcomes"]
    m = json.loads((out_dir / "smarts-optin-310k.json").read_text())
    s["smarts_310k"] = {"counts": m["counts"], "unexpected_total": m["unexpected_total"]}
    for name in CORPORA:
        c = json.loads((out_dir / f"chemistry-{name}.json").read_text())
        s[name] = {
            "cip_accurate": c["cip_accurate"],
            "hybridization": c["hybridization"],
            "mmff94": c["mmff94"],
            "mol_writer": c["mol_writer"],
        }
    if ran_corpus:
        b = json.loads((out_dir / "biotransformer-corpus.json").read_text())
        s["biotransformer_corpus"] = b["counts"]
    return s


def diff(expected, got, path=""):
    if isinstance(expected, dict) and isinstance(got, dict):
        out = []
        for k in sorted(set(expected) | set(got)):
            out += diff(expected.get(k), got.get(k), f"{path}.{k}" if path else k)
        return out
    return [] if expected == got else [f"{path}: expected {expected!r}, got {got!r}"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--wheel", type=Path, required=True, help="the published wheel file installed in --chematic-python")
    ap.add_argument("--chematic-python", default=sys.executable)
    ap.add_argument("--rdkit-python", default=sys.executable)
    ap.add_argument("--rules-dir", type=Path, help="directory with the three BioTransformer rule tables")
    ap.add_argument("--expected", type=Path, required=True, help="validation/published-wheel-chemistry-gates-expected-v<version>.json")
    ap.add_argument("--scope", required=True, help="label, e.g. published_pypi_v1.0.36_macos_arm64_cp312")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--resume", action="store_true",
                    help="skip steps whose output in --out-dir is already complete; rerun the same "
                         "command on hosts that end long-lived processes until every step is done")
    ap.add_argument("--corpus-shards", type=int, default=1,
                    help="run the BioTransformer corpus as this many reactant shards (each its own "
                         "process and resumable output), then add their counts up")
    args = ap.parse_args()
    if args.corpus_shards < 1:
        ap.error("--corpus-shards must be at least 1")
    # Steps run from the repository root; anchor the caller's paths first.
    args.wheel = args.wheel.resolve()
    args.expected = args.expected.resolve()
    if args.rules_dir:
        args.rules_dir = args.rules_dir.resolve()

    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    log = out / "gates.log"
    if not (args.resume and log.exists()):
        log.write_text("", encoding="utf-8")
    resume = args.resume
    cp, rp = args.chematic_python, args.rdkit_python
    version = subprocess.run(
        [cp, "-c", "import chematic; print(chematic.__version__)"], check=True, capture_output=True, text=True
    ).stdout.strip()

    run([cp, str(SCRIPTS / "reaction_python_checked_candidates.py"), "--wheel", str(args.wheel),
         "--published-from", f"https://pypi.org/project/chematic/{version}/",
         "--output", str(out / "reaction-candidates-83.json")], log,
        out / "reaction-candidates-83.json", resume)
    run([rp, str(SCRIPTS / "reaction_python_checked_provenance_gate.py"),
         "--candidates", str(out / "reaction-candidates-83.json"),
         "--output", str(out / "reaction-checked-provenance-83.json")], log,
        out / "reaction-checked-provenance-83.json", resume)
    run([cp, str(SCRIPTS / "check_python_smarts_parity_310k.py"), "--wheel", str(args.wheel),
         "--module-root", str(module_root(cp)), "--archived-oracle", "--scope", args.scope,
         "--output", str(out / "smarts-optin-310k.json")], log, out / "smarts-optin-310k.json", resume)
    for name, corpus in CORPORA.items():
        result = out / f"chemistry-{name}.json"
        if resume and finished(result):
            continue
        run([cp, str(SCRIPTS / "chematic_chemistry_dump.py"), "--corpus", str(corpus),
             "--output", str(out / f"dump-{name}.jsonl")], log)
        run([rp, str(SCRIPTS / "compare_chemistry_dump_rdkit.py"), "--dump", str(out / f"dump-{name}.jsonl"),
             "--output", str(result)], log)
        (out / f"dump-{name}.jsonl").unlink()
    ran_corpus = False
    if args.rules_dir:
        rules = []
        for fname, digest in RULE_FILES.items():
            p = args.rules_dir / fname
            if sha256(p) != digest:
                raise SystemExit(f"{p} does not hash to the pinned {digest}")
            rules.append(str(p))
        reactant_file = reactants(out / "bt_reactants_400.smi")
        n_reactants = len(reactant_file.read_text(encoding="utf-8").splitlines())
        if args.corpus_shards == 1:
            run([rp, str(SCRIPTS / "biotransformer_rule_corpus.py"), "--rules", *rules,
                 "--reactants", str(reactant_file),
                 "--chematic-python", cp,
                 "--output", str(out / "biotransformer-corpus.json"),
                 "--rows", str(out / "biotransformer-corpus-rows.jsonl")], log,
                out / "biotransformer-corpus.json", resume)
        else:
            shard_dir = out / "biotransformer-shards"
            shard_dir.mkdir(exist_ok=True)
            n = args.corpus_shards
            shards = []
            for k in range(n):
                start, end = k * n_reactants // n, (k + 1) * n_reactants // n
                summary_path = shard_dir / f"shard-{k + 1}-of-{n}.json"
                rows_path = shard_dir / f"shard-{k + 1}-of-{n}-rows.jsonl"
                run([rp, str(SCRIPTS / "biotransformer_rule_corpus.py"), "--rules", *rules,
                     "--reactants", str(reactant_file), "--reactants-slice", f"{start}:{end}",
                     "--chematic-python", cp,
                     "--output", str(summary_path), "--rows", str(rows_path)], log, summary_path, resume)
                shards.append((summary_path, rows_path))
            merge_corpus_shards(shards, out / "biotransformer-corpus.json",
                                out / "biotransformer-corpus-rows.jsonl")
        ran_corpus = True

    got = summarize(out, ran_corpus)
    expected = json.loads(args.expected.read_text())["gates"]
    if not ran_corpus:
        expected = {k: v for k, v in expected.items() if k != "biotransformer_corpus"}
    problems = diff(expected, got)
    report = {
        "schema": "published-wheel-chemistry-gates/v1",
        "scope": args.scope,
        "chematic_version": version,
        "wheel": args.wheel.name,
        "wheel_sha256": sha256(args.wheel),
        "host": {"platform": platform.platform(), "machine": platform.machine(),
                 "python": platform.python_version()},
        "expected": str(args.expected.relative_to(ROOT)) if args.expected.is_relative_to(ROOT) else str(args.expected),
        "gates": got,
        "differences": problems,
        "pass": not problems,
    }
    (out / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"scope": args.scope, "pass": not problems, "differences": problems}, indent=2))
    return 0 if not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())
