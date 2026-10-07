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

``--phase chematic`` runs only chematic's side and writes its outputs (the
reaction candidates, the SMARTS cells, the chemistry dumps and, with
``--rules-dir``, the BioTransformer responses) to ``--chematic-outputs``;
``--phase rdkit`` compares those outputs with RDKit on another host. The
workflow runs the first on the macOS and Windows runners (RDKit 2026.03.6
ships no macOS x86-64 wheel) and the second on Linux against one pinned
RDKit build. ``--phase all`` (the default) runs both here.
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
REQUESTS = ROOT / "validation/biotransformer-requests-400.tsv"
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


def rule_paths(rules_dir: Path) -> list[str]:
    rules = []
    for fname, digest in RULE_FILES.items():
        p = rules_dir / fname
        if sha256(p) != digest:
            raise SystemExit(f"{p} does not hash to the pinned {digest}")
        rules.append(str(p))
    return rules


def chematic_phase(args, cp: str, out: Path, log: Path, resume: bool) -> dict:
    """chematic's side of every gate; RDKit is not imported."""
    version = subprocess.run(
        [cp, "-c", "import chematic; print(chematic.__version__)"], check=True, capture_output=True, text=True
    ).stdout.strip()
    run([cp, str(SCRIPTS / "reaction_python_checked_candidates.py"), "--wheel", str(args.wheel),
         "--published-from", f"https://pypi.org/project/chematic/{version}/",
         "--output", str(out / "reaction-candidates-83.json")], log,
        out / "reaction-candidates-83.json", resume)
    run([cp, str(SCRIPTS / "check_python_smarts_parity_310k.py"), "--wheel", str(args.wheel),
         "--module-root", str(module_root(cp)), "--archived-oracle", "--scope", args.scope,
         "--output", str(out / "smarts-optin-310k.json")], log, out / "smarts-optin-310k.json", resume)
    for name, corpus in CORPORA.items():
        dump = out / f"dump-{name}.jsonl"
        if resume and (out / f"dump-{name}.done").is_file():
            continue
        run([cp, str(SCRIPTS / "chematic_chemistry_dump.py"), "--corpus", str(corpus),
             "--output", str(dump)], log)
        (out / f"dump-{name}.done").write_text("", encoding="utf-8")
    ran_corpus = False
    shards = 0
    if args.rules_dir:
        # The rules in --corpus-shards slices, each its own process and
        # resumable file.
        n_rules = int(subprocess.run(
            [cp, "-c", "import sys; sys.path.insert(0, sys.argv[1]); from biotransformer_rules import "
             "load_rules; from pathlib import Path; print(len(load_rules([Path(p) for p in sys.argv[2:]])[0]))",
             str(SCRIPTS), *rule_paths(args.rules_dir)],
            check=True, capture_output=True, text=True).stdout)
        shards = args.corpus_shards
        for k in range(shards):
            done = out / f"biotransformer-responses-{k + 1}-of-{shards}.done"
            if resume and done.is_file():
                continue
            start, end = k * n_rules // shards, (k + 1) * n_rules // shards
            run([cp, str(SCRIPTS / "biotransformer_chematic_responses.py"),
                 "--rules", *rule_paths(args.rules_dir), "--requests", str(REQUESTS),
                 "--rules-slice", f"{start}:{end}",
                 "--output", str(out / f"biotransformer-responses-{k + 1}-of-{shards}.jsonl.gz")], log)
            done.write_text("", encoding="utf-8")
        ran_corpus = True
    manifest = {
        "schema": "published-wheel-chematic-outputs/v1",
        "scope": args.scope,
        "chematic_version": version,
        "wheel": args.wheel.name,
        "wheel_sha256": sha256(args.wheel),
        "host": {"platform": platform.platform(), "machine": platform.machine(),
                 "python": platform.python_version()},
        "biotransformer_corpus": ran_corpus,
        "biotransformer_response_shards": shards,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def rdkit_phase(args, rp: str, src: Path, out: Path, log: Path, resume: bool) -> int:
    """RDKit's side: compare chematic's outputs in ``src`` and check the counts."""
    manifest = json.loads((src / "manifest.json").read_text(encoding="utf-8"))
    run([rp, str(SCRIPTS / "reaction_python_checked_provenance_gate.py"),
         "--candidates", str(src / "reaction-candidates-83.json"),
         "--output", str(out / "reaction-checked-provenance-83.json")], log,
        out / "reaction-checked-provenance-83.json", resume)
    if src != out:
        (out / "smarts-optin-310k.json").write_bytes((src / "smarts-optin-310k.json").read_bytes())
    for name in CORPORA:
        result = out / f"chemistry-{name}.json"
        run([rp, str(SCRIPTS / "compare_chemistry_dump_rdkit.py"), "--dump", str(src / f"dump-{name}.jsonl"),
             "--output", str(result)], log, result, resume)
    ran_corpus = manifest["biotransformer_corpus"]
    if ran_corpus:
        if not args.rules_dir:
            raise SystemExit("the chematic outputs include the BioTransformer corpus; pass --rules-dir")
        reactant_file = reactants(out / "bt_reactants_400.smi")
        n = manifest["biotransformer_response_shards"]
        run([rp, str(SCRIPTS / "biotransformer_rule_corpus.py"), "--rules", *rule_paths(args.rules_dir),
             "--reactants", str(reactant_file),
             "--replay", *(str(src / f"biotransformer-responses-{k + 1}-of-{n}.jsonl.gz") for k in range(n)),
             "--requests", str(REQUESTS),
             "--output", str(out / "biotransformer-corpus.json"),
             "--rows", str(out / "biotransformer-corpus-rows.jsonl")], log,
            out / "biotransformer-corpus.json", resume)
    got = summarize(out, ran_corpus)
    expected = json.loads(args.expected.read_text())["gates"]
    if not ran_corpus:
        expected = {k: v for k, v in expected.items() if k != "biotransformer_corpus"}
    problems = diff(expected, got)
    report = {
        "schema": "published-wheel-chemistry-gates/v1",
        "scope": args.scope,
        "chematic_version": manifest["chematic_version"],
        "wheel": manifest["wheel"],
        "wheel_sha256": manifest["wheel_sha256"],
        "host": manifest["host"],
        "rdkit_host": {"platform": platform.platform(), "machine": platform.machine(),
                       "python": platform.python_version()},
        "expected": str(args.expected.relative_to(ROOT)) if args.expected.is_relative_to(ROOT) else str(args.expected),
        "gates": got,
        "differences": problems,
        "pass": not problems,
    }
    (out / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"scope": args.scope, "pass": not problems, "differences": problems}, indent=2))
    return 0 if not problems else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--phase", choices=("all", "chematic", "rdkit"), default="all")
    ap.add_argument("--wheel", type=Path, help="the published wheel file installed in --chematic-python "
                    "(phases all and chematic)")
    ap.add_argument("--chematic-python", default=sys.executable)
    ap.add_argument("--rdkit-python", default=sys.executable)
    ap.add_argument("--rules-dir", type=Path, help="directory with the three BioTransformer rule tables")
    ap.add_argument("--expected", type=Path, help="validation/published-wheel-chemistry-gates-expected-v<version>.json "
                    "(phases all and rdkit)")
    ap.add_argument("--scope", required=True, help="label, e.g. published_pypi_v1.0.36_macos_arm64_cp312")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--chematic-outputs", type=Path,
                    help="where --phase chematic wrote its outputs (default: --out-dir)")
    ap.add_argument("--resume", action="store_true",
                    help="skip steps whose output in --out-dir is already complete; rerun the same "
                         "command on hosts that end long-lived processes until every step is done")
    ap.add_argument("--corpus-shards", type=int, default=1,
                    help="record chematic's BioTransformer responses in this many rule slices (each "
                         "its own process and resumable file)")
    args = ap.parse_args()
    if args.corpus_shards < 1:
        ap.error("--corpus-shards must be at least 1")
    if args.phase in ("all", "chematic") and not args.wheel:
        ap.error("--wheel is required for this phase")
    if args.phase in ("all", "rdkit") and not args.expected:
        ap.error("--expected is required for this phase")
    # Steps run from the repository root; anchor the caller's paths first.
    if args.wheel:
        args.wheel = args.wheel.resolve()
    if args.expected:
        args.expected = args.expected.resolve()
    if args.rules_dir:
        args.rules_dir = args.rules_dir.resolve()

    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    log = out / "gates.log"
    if not (args.resume and log.exists()):
        log.write_text("", encoding="utf-8")
    src = args.chematic_outputs.resolve() if args.chematic_outputs else out
    if args.phase in ("all", "chematic"):
        chematic_phase(args, args.chematic_python, src, log, args.resume)
    if args.phase == "chematic":
        return 0
    return rdkit_phase(args, args.rdkit_python, src, out, log, args.resume)


if __name__ == "__main__":
    raise SystemExit(main())
