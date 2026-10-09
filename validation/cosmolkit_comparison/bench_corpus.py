#!/usr/bin/env python3
"""Time the same operations in chematic, COSMolKit and RDKit (Python bindings).

Each engine runs in its own interpreter. A block runs every engine once, in
an order that rotates between blocks; per operation the median over blocks is
reported. Molecules are re-parsed (untimed) before each operation so that a
result cached on a molecule object cannot be reused. Calls return each
library's native result object (no normalisation), and SMARTS queries are
compiled once where the library offers a compiled query.

Timings depend on the host; the report records it. Speed is only meaningful
next to the correctness summary from ``compare_corpus.py``: an operation whose
outputs differ is a different amount of work.

    python bench_corpus.py --corpus CORPUS.smi --output bench.json --blocks 5
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import statistics
import subprocess
import sys
import time
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
OPS = ["parse", "canonical_smiles", "formula", "mol_wt", "tpsa", "logp", "qed",
       "morgan2_2048", "maccs", "cip", "smarts_31", "molblock", "inchi",
       "smarts_write", "murcko", "stereoisomers"]
# Operations the COSMolKit 0.5 table offers (``--ops all05``): every op both
# chematic 1.0.40 and COSMolKit 0.5 support with RDKit-equivalent results.
OPS_05 = ["parse", "canonical_smiles", "smiles_kekule", "formula", "mol_wt", "exact_mw", "tpsa",
          "logp", "qed", "hba", "hbd", "rotatable_bonds", "morgan2_2048", "morgan3_2048",
          "morgan2_chiral", "maccs", "fp_atom_pair", "fp_torsion", "fp_pattern", "fp_layered",
          "fp_rdkit", "cip", "smarts_31", "molblock", "smarts_write", "cx_smarts", "stereoisomers",
          "stereoisomer_count", "chi1v", "kappa2", "labute_asa", "slogp_vsa", "mqn", "num_rings",
          "add_hs_smiles", "distance_matrix", "tautomer_canonical", "tautomer_set",
          "rxn_amide_amine", "mmff_energy_gradient", "uff_energy_gradient"]
# Ops timed on explicit-H molecules carrying RDKit-embedded coordinates
# (prepared untimed): fn(prepared) where prepared = prep3d(mol, coords).
THREE_D = {"mmff_energy_gradient", "uff_energy_gradient"}
RXN_AMIDE = "[N;!H0;!$(NC=O);!$(N=*);!$(N#*);!$(Nc);!$([N+]):1].[C:2](=[O:3])[OH]>>[N:1][C:2]=[O:3]"


def table(engine: str):
    from corpus_engines import SMARTS_QUERIES
    if engine == "rdkit":
        from rdkit import Chem, RDLogger
        from rdkit.Chem import QED, Crippen, Descriptors, MACCSkeys, rdCIPLabeler, rdDepictor, rdFingerprintGenerator
        from rdkit.Chem import rdMolDescriptors as D
        from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
        from rdkit.Chem.Scaffolds import MurckoScaffold
        RDLogger.DisableLog("rdApp.*")
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        qs = [Chem.MolFromSmarts(q) for q in SMARTS_QUERIES]

        def molblock_2d(m):
            # 2D layout + MOL block, as chematic's rdkit_mol_block_2d and
            # COSMolKit's to_2d_sdf_string do (MolToMolBlock alone writes zeros).
            rdDepictor.Compute2DCoords(m)
            return Chem.MolToMolBlock(m)

        from rdkit.Chem import AllChem, GraphDescriptors
        from rdkit.Chem.MolStandardize import rdMolStandardize
        from rdkit.Chem.EnumerateStereoisomers import GetStereoisomerCount
        gen3 = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048)
        genc = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
        taut = rdMolStandardize.TautomerEnumerator()
        rxn = AllChem.ReactionFromSmarts(RXN_AMIDE)
        acid = Chem.MolFromSmiles("CC(=O)O")

        def run_rxn(m):
            out = []
            for ps in rxn.RunReactants((m, acid)):
                for p in ps:
                    Chem.SanitizeMol(p)
                out.append([Chem.MolToSmiles(p) for p in ps])
            return out

        def prep3d(m, xyz):
            mh = Chem.AddHs(m)
            conf = Chem.Conformer(mh.GetNumAtoms())
            for i, p in enumerate(xyz):
                conf.SetAtomPosition(i, p)
            mh.AddConformer(conf, assignId=True)
            return mh

        def mmff(mh):
            ff = AllChem.MMFFGetMoleculeForceField(mh, AllChem.MMFFGetMoleculeProperties(mh))
            return ff.CalcEnergy(), ff.CalcGrad()

        def uff(mh):
            ff = AllChem.UFFGetMoleculeForceField(mh)
            return ff.CalcEnergy(), ff.CalcGrad()

        extra = {
            "smiles_kekule": lambda m: Chem.MolToSmiles(m, kekuleSmiles=True),
            "exact_mw": Descriptors.ExactMolWt, "hba": D.CalcNumHBA, "hbd": D.CalcNumHBD,
            "rotatable_bonds": Descriptors.NumRotatableBonds, "morgan3_2048": gen3.GetFingerprint,
            "morgan2_chiral": genc.GetFingerprint,
            "fp_atom_pair": lambda m: D.GetHashedAtomPairFingerprintAsBitVect(m, nBits=2048),
            "fp_torsion": lambda m: D.GetHashedTopologicalTorsionFingerprintAsBitVect(m, nBits=2048),
            "fp_pattern": lambda m: Chem.PatternFingerprint(m, fpSize=2048),
            "fp_layered": lambda m: Chem.LayeredFingerprint(m, fpSize=2048),
            "fp_rdkit": lambda m: Chem.RDKFingerprint(m, fpSize=2048),
            "cx_smarts": Chem.MolToCXSmarts, "stereoisomer_count": GetStereoisomerCount,
            "chi1v": D.CalcChi1v, "kappa2": D.CalcKappa2, "labute_asa": D.CalcLabuteASA,
            "slogp_vsa": D.SlogP_VSA_, "mqn": D.MQNs_, "num_rings": D.CalcNumRings,
            "add_hs_smiles": lambda m: Chem.MolToSmiles(Chem.AddHs(m)),
            "distance_matrix": Chem.GetDistanceMatrix,
            "tautomer_canonical": lambda m: Chem.MolToSmiles(taut.Canonicalize(m)),
            "tautomer_set": lambda m: sorted(Chem.MolToSmiles(t) for t in taut.Enumerate(m)),
            "rxn_amide_amine": run_rxn, "mmff_energy_gradient": mmff, "uff_energy_gradient": uff,
            "__prep3d": prep3d,
        }
        return Chem.MolFromSmiles, extra | {
            "canonical_smiles": Chem.MolToSmiles, "formula": D.CalcMolFormula,
            "mol_wt": Descriptors.MolWt, "tpsa": D.CalcTPSA, "logp": Crippen.MolLogP,
            "qed": QED.qed, "morgan2_2048": gen.GetFingerprint, "maccs": MACCSkeys.GenMACCSKeys,
            "cip": rdCIPLabeler.AssignCIPLabels,
            "smarts_31": lambda m: [m.GetSubstructMatches(q, maxMatches=100000) for q in qs],
            "molblock": molblock_2d, "inchi": Chem.MolToInchi,
            "smarts_write": Chem.MolToSmarts, "murcko": MurckoScaffold.GetScaffoldForMol,
            "stereoisomers": lambda m: sorted(Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m)),
        }
    if engine == "chematic":
        import chematic as c
        acid = c.from_smiles("CC(=O)O")
        extra = {
            "smiles_kekule": lambda m: m.rdkit_smiles_with(kekule=True),
            "exact_mw": lambda m: m.exact_mass, "hba": lambda m: m.rdkit_hba, "hbd": lambda m: m.hbd,
            "rotatable_bonds": lambda m: m.rotatable_bonds,
            "morgan3_2048": lambda m: m.rdkit_ecfp_config(3, 2048),
            "morgan2_chiral": lambda m: m.rdkit_ecfp_config(2, 2048, include_chirality=True),
            "fp_atom_pair": lambda m: m.rdkit_atom_pair_fp(),
            "fp_torsion": lambda m: m.rdkit_torsion_fp(),
            "fp_pattern": lambda m: m.rdkit_pattern_fp(),
            "fp_layered": lambda m: m.rdkit_layered_fp(),
            "fp_rdkit": lambda m: m.rdkit_rdk_fp(),
            "cx_smarts": lambda m: m.rdkit_cx_smarts(),
            "stereoisomer_count": lambda m: m.rdkit_stereoisomer_count(),
            "chi1v": lambda m: m.chi1v, "kappa2": lambda m: m.kappa2, "labute_asa": lambda m: m.labute_asa,
            "slogp_vsa": lambda m: m.slogp_vsa(), "mqn": lambda m: m.mqn(), "num_rings": lambda m: m.num_rings,
            "add_hs_smiles": lambda m: m.add_hydrogens().rdkit_smiles,
            "distance_matrix": lambda m: m.topological_distance_matrix(),
            "tautomer_canonical": lambda m: m.canonical_tautomer().smiles,
            "tautomer_set": lambda m: sorted(t.smiles for t in m.enumerate_tautomers()),
            "rxn_amide_amine": lambda m: c.run_smirks_checked(RXN_AMIDE, [m, acid], rdkit_compat=True),
            "mmff_energy_gradient": lambda p: p[0]._rdkit_mmff_terms(p[1]),
            "uff_energy_gradient": lambda p: (p[0].rdkit_uff_energy(p[1]), p[0].rdkit_uff_gradient(p[1])),
            "__prep3d": lambda m, xyz: (m.add_hydrogens(), xyz),
        }
        return c.from_smiles, extra | {
            "canonical_smiles": lambda m: m.rdkit_smiles, "formula": lambda m: m.formula,
            "mol_wt": lambda m: m.rdkit_mw, "tpsa": lambda m: m.rdkit_tpsa,
            "logp": lambda m: m.logp, "qed": lambda m: m.qed,
            "morgan2_2048": lambda m: m.rdkit_ecfp_config(2, 2048), "maccs": lambda m: m.maccs(),
            "cip": lambda m: m.cip_stereo(mode="accurate"),
            "smarts_31": lambda m: [m.find_matches_rdkit_parity(q) for q in SMARTS_QUERIES],
            "molblock": lambda m: m.rdkit_mol_block_2d(), "inchi": lambda m: m.standard_inchi,
            "smarts_write": lambda m: m.rdkit_smarts(), "murcko": lambda m: m.rdkit_murcko_scaffold(),
            "stereoisomers": lambda m: m.rdkit_stereoisomer_smiles(),
        }
    import cosmolkit as ck
    if not ck.__version__.startswith("0.3"):
        return _cosmolkit_05_table(ck)
    qs = [ck.parse_smarts(q) for q in SMARTS_QUERIES]
    return ck.Molecule.from_smiles, {
        "canonical_smiles": lambda m: m.to_smiles(), "formula": ck.calc_mol_formula,
        "mol_wt": ck.calc_mol_wt, "tpsa": ck.calc_tpsa,
        "logp": lambda m: ck.calc_crippen_descriptors(m)[0], "qed": ck.calc_qed,
        "morgan2_2048": lambda m: m.fingerprint_morgan(radius=2, n_bits=2048),
        "maccs": lambda m: m.maccs_fingerprint(), "cip": lambda m: m.with_cip_labels(),
        "smarts_31": lambda m: [ck.get_substruct_matches(m, q, max_matches=100000) for q in qs],
        "molblock": lambda m: m.to_2d_sdf_string(), "inchi": lambda m: m.to_inchi(),
        "smarts_write": lambda m: m.to_smarts(), "murcko": lambda m: m.murcko_scaffold(),
        "stereoisomers": lambda m: sorted(x.to_smiles() for x in m.stereoisomers()),
    }


def _cosmolkit_05_table(ck):
    """COSMolKit 0.5's reworked API (generator fingerprints, *_with_params)."""
    from corpus_engines import SMARTS_QUERIES
    qs = [ck.parse_smarts(q) for q in SMARTS_QUERIES]
    gens = {key: ck.MorganFingerprintGenerator(params=ck.MorganParams(fp_size=2048, **kw))
            for key, kw in (("m2", {"radius": 2}), ("m3", {"radius": 3}),
                            ("m2c", {"radius": 2, "include_chirality": True}))}
    kek = ck.SmilesWriteParams(do_kekule=True)
    rxn = ck.parse_smirks(RXN_AMIDE)
    acid = ck.Molecule.from_smiles("CC(=O)O")
    run_params = ck.ReactionRunParams()

    def mmff(m):
        g = m.mmff_energy_gradient()
        return g.energy, g.gradient

    def uff(m):
        g = m.uff_energy_gradient()
        return g.energy, g.gradient

    return ck.Molecule.from_smiles, {
        "canonical_smiles": lambda m: m.to_smiles(),
        "smiles_kekule": lambda m: m.to_smiles_with_params(kek),
        "formula": lambda m: m.molecular_formula(), "mol_wt": lambda m: m.molecular_weight(),
        "exact_mw": lambda m: m.exact_molecular_weight(), "tpsa": lambda m: m.tpsa(),
        "logp": lambda m: m.crippen_descriptors().logp, "qed": lambda m: m.qed(),
        "hba": lambda m: m.num_hba(), "hbd": lambda m: m.num_hbd(),
        "rotatable_bonds": lambda m: m.num_rotatable_bonds(),
        "morgan2_2048": lambda m: m.morgan_fingerprint_with_generator(gens["m2"]),
        "morgan3_2048": lambda m: m.morgan_fingerprint_with_generator(gens["m3"]),
        "morgan2_chiral": lambda m: m.morgan_fingerprint_with_generator(gens["m2c"]),
        "maccs": lambda m: m.maccs_fingerprint(),
        "fp_atom_pair": lambda m: m.atom_pair_fingerprint(),
        "fp_torsion": lambda m: m.legacy_topological_torsion_fingerprint(),
        "fp_pattern": lambda m: m.pattern_fingerprint(), "fp_layered": lambda m: m.layered_fingerprint(),
        "fp_rdkit": lambda m: m.topological_fingerprint(),
        "cip": lambda m: m.with_cip_labels(),
        "smarts_31": lambda m: [m.substruct_matches(q) for q in qs],
        "molblock": lambda m: m.to_sdf_2d(), "smarts_write": lambda m: m.to_smarts(),
        "cx_smarts": lambda m: m.to_cx_smarts(),
        "stereoisomers": lambda m: sorted(x.to_smiles() for x in m.enumerate_stereoisomers()),
        "stereoisomer_count": lambda m: m.stereoisomer_count(),
        "chi1v": lambda m: m.chi_1_v(), "kappa2": lambda m: m.kappa_2(),
        "labute_asa": lambda m: m.labute_asa(), "slogp_vsa": lambda m: m.slogp_vsa(),
        "mqn": lambda m: m.mqns(), "num_rings": lambda m: m.num_rings(),
        "add_hs_smiles": lambda m: m.with_hydrogens().to_smiles(),
        "distance_matrix": lambda m: m.distance_matrix(),
        "tautomer_canonical": lambda m: m.canonical_tautomer().to_smiles(),
        "tautomer_set": lambda m: m.enumerate_tautomers().canonical_smiles(),
        "rxn_amide_amine": lambda m: rxn.run([m, acid], run_params),
        "mmff_energy_gradient": mmff, "uff_energy_gradient": uff,
        "__prep3d": lambda m, xyz: m.with_hydrogens().with_only_3d_conformer(xyz),
    }


def measure(engine: str, corpus: Path, ops: list[str]) -> dict:
    warnings.simplefilter("ignore")
    parse, fns = table(engine)
    rows = [l.split()[0] for l in corpus.read_text().splitlines() if l.strip()]

    def parsed():
        out = []
        for s in rows:
            try:
                m = parse(s)
            except Exception:  # noqa: BLE001
                m = None
            if m is not None:
                out.append(m)
        return out

    coords = None
    if THREE_D & set(ops):
        # RDKit-embedded coordinates of every input (untimed, shared by the
        # engines through the same seeded embedding).
        from corpus_engines import rdkit_embedded_h
        coords = [(e[1] if (e := rdkit_embedded_h(s)) is not None else None) for s in rows]

    def prepared():
        out = []
        for s, xyz in zip(rows, coords):
            if xyz is None:
                continue
            try:
                m = parse(s)
                out.append(fns["__prep3d"](m, xyz))
            except Exception:  # noqa: BLE001
                pass
        return out

    res = {}
    t = time.perf_counter()
    mols = parsed()
    res["parse"] = time.perf_counter() - t
    for op in ops:
        if op == "parse":
            continue
        fn = fns.get(op)
        if fn is None:
            res[op] = None
            continue
        mols = prepared() if op in THREE_D else parsed()
        errors = 0
        t = time.perf_counter()
        for m in mols:
            try:
                fn(m)
            except Exception:  # noqa: BLE001 - counted, still timed
                errors += 1
        res[op] = time.perf_counter() - t
        if errors:
            res[op + "__errors"] = errors
    res["parsed"] = len(mols)
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--corpus", type=Path, required=True)
    ap.add_argument("--output", type=Path)
    ap.add_argument("--blocks", type=int, default=5)
    ap.add_argument("--ops", default=",".join(OPS))
    ap.add_argument("--python", action="append", metavar="ENGINE=PYTHON", default=[],
                    help="interpreter per engine (default: this one)")
    ap.add_argument("--_child")
    args = ap.parse_args()
    ops = OPS_05 if args.ops == "all05" else args.ops.split(",")
    if args._child:
        print(json.dumps(measure(args._child, args.corpus, ops)))
        return 0
    pythons = {"rdkit": sys.executable, "chematic": sys.executable, "cosmolkit": sys.executable}
    pythons.update(dict(p.split("=", 1) for p in args.python))
    engines = list(pythons)
    samples: dict[str, list[dict]] = {e: [] for e in engines}
    load_start = list(os.getloadavg())
    versions = {}
    for block in range(args.blocks):
        order = engines[block % len(engines):] + engines[: block % len(engines)]
        for e in order:
            out = subprocess.run([pythons[e], __file__, "--corpus", str(args.corpus), "--ops", ",".join(ops),
                                  "--_child", e], capture_output=True, text=True, check=True)
            samples[e].append(json.loads(out.stdout))
        print(f"block {block + 1}/{args.blocks} done", file=sys.stderr)
    for e in engines:
        code = {"rdkit": "from rdkit import rdBase;print(rdBase.rdkitVersion)",
                "chematic": "import chematic;print(chematic.__version__)",
                "cosmolkit": "import cosmolkit;print(cosmolkit.__version__)"}[e]
        versions[e] = subprocess.run([pythons[e], "-c", code], capture_output=True, text=True).stdout.strip()
    report = {"schema": "cosmolkit-corpus-bench/v1", "corpus": args.corpus.name,
              "blocks": args.blocks, "host": {"platform": platform.platform(),
                                               "cpus": os.cpu_count(), "python": platform.python_version()},
              "versions": versions, "median_seconds": {}, "parsed": {}, "errors": {},
              "load_average_start": load_start, "load_average_end": list(os.getloadavg())}
    for e in engines:
        report["parsed"][e] = samples[e][0]["parsed"]
        report["median_seconds"][e] = {
            op: (None if samples[e][0].get(op) is None
                 else statistics.median(s[op] for s in samples[e])) for op in ops}
        report["errors"][e] = {k: v for k, v in samples[e][0].items() if k.endswith("__errors")}
    ck, ch, rd = (report["median_seconds"][e] for e in ("cosmolkit", "chematic", "rdkit"))
    report["chematic_speedup_vs_cosmolkit"] = {
        op: (round(ck[op] / ch[op], 2) if ck.get(op) and ch.get(op) else None) for op in ops}
    report["chematic_speedup_vs_rdkit"] = {
        op: (round(rd[op] / ch[op], 2) if rd.get(op) and ch.get(op) else None) for op in ops}
    report["cosmolkit_speedup_vs_rdkit"] = {
        op: (round(rd[op] / ck[op], 2) if rd.get(op) and ck.get(op) else None) for op in ops}
    text = json.dumps(report, indent=1, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(f"{'op':18}{'chematic ms':>13}{'COSMolKit ms':>14}{'RDKit ms':>11}{'ch/ck':>8}")
    for op in ops:
        f = lambda v: f"{v * 1000:.0f}" if v else "-"  # noqa: E731
        print(f"{op:18}{f(ch.get(op)):>13}{f(ck.get(op)):>14}{f(rd.get(op)):>11}"
              f"{str(report['chematic_speedup_vs_cosmolkit'][op] or '-') + 'x':>8}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
