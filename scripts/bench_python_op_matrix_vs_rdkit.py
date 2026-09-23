#!/usr/bin/env python3
"""Per-operation Python timing matrix: chematic vs RDKit, with output checks.

For every operation the same corpus rows are timed in both libraries,
single-threaded, in the same process. Inputs are rebuilt before every repeat
(construction is not timed) so neither library can reuse memoized results
from a previous repeat:

* chematic ``Mol`` inputs are re-parsed from SMILES;
* RDKit ``Mol`` inputs are re-parsed with ``MolFromSmiles`` (``Chem.Mol(m)``
  would copy cached computed properties such as Crippen/TPSA values).

RDKit perceives rings and aromaticity during sanitization, i.e. inside
``MolFromSmiles``; chematic perceives them lazily on first use. Prepared-input
rows therefore charge that perception to chematic only; the ``parse+...`` rows
include parsing on both sides.

Each operation carries an ``equivalence`` class that decides whether it may be
counted as a win:

``checked``        outputs compared row-by-row (untimed pass); counted only
                   when agreement is reported
``same-task``      same task, different but legitimate representation
                   (e.g. canonical SMILES spelling); not counted
``not-equivalent`` different algorithm, scope or definition; never counted

Nothing here is a release or universal speed claim; see docs/benchmark.md.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import statistics
import subprocess
import sys
import time
from datetime import datetime, timezone

import chematic
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import (
    BRICS,
    QED,
    AllChem,
    Crippen,
    Descriptors,
    MACCSkeys,
    rdCIPLabeler,
    rdDepictor,
    rdFMCS,
    rdMolDescriptors,
)
from rdkit.Chem import rdFingerprintGenerator as rfg
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def cpu_model() -> str:
    try:
        with open("/proc/cpuinfo", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        pass
    try:
        return subprocess.run(
            ["sysctl", "-n", "machdep.cpu.brand_string"],
            capture_output=True, text=True, check=False,
        ).stdout.strip() or platform.processor()
    except OSError:
        return platform.processor()


def bits_of_rdkit(fp) -> bytes:
    """RDKit ExplicitBitVect -> chematic byte layout (bit i -> byte i//8, LSB first)."""
    n = fp.GetNumBits()
    out = bytearray((n + 7) // 8)
    for i in fp.GetOnBits():
        out[i // 8] |= 1 << (i % 8)
    return bytes(out)


def close(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corpus", required=True)
    ap.add_argument("--limit", type=int, default=5000)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--only", default=None, help="substring filter on operation names")
    ap.add_argument("--output-json", required=True)
    ap.add_argument("--chematic-revision", required=True, help="git revision the chematic wheel was built from")
    ap.add_argument("--chematic-artifact", default=None, help="path of the installed wheel, hashed into the record")
    args = ap.parse_args()

    smis = [l.split()[0] for l in open(args.corpus, encoding="utf-8") if l.strip()][: args.limit]
    pairs = []
    c_fail = r_fail = 0
    for s in smis:
        try:
            cm = chematic.from_smiles(s)
        except Exception:
            c_fail += 1
            continue
        rm = Chem.MolFromSmiles(s)
        if rm is None:
            r_fail += 1
            continue
        pairs.append((s, cm, rm))
    S = [p[0] for p in pairs]
    CM = [p[1] for p in pairs]
    RM = [p[2] for p in pairs]
    n = len(pairs)

    record = {
        "schema": "chematic-python-op-matrix-vs-rdkit/v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "chematic": {
            "version": chematic.__version__,
            "revision": args.chematic_revision,
            "artifact": os.path.basename(args.chematic_artifact) if args.chematic_artifact else None,
            "artifact_sha256": sha256_file(args.chematic_artifact) if args.chematic_artifact else None,
        },
        "rdkit": {"version": Chem.rdBase.rdkitVersion},
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "machine": platform.machine(),
            "cpu_model": cpu_model(),
            "cpu_count": os.cpu_count(),
            "threads": "single-threaded calls in one process",
        },
        "corpus": {
            "path": os.path.relpath(args.corpus, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
            "sha256": sha256_file(args.corpus),
            "rows_read": len(smis),
            "common_valid_rows": n,
            "chematic_parse_failures": c_fail,
            "rdkit_parse_failures_after_chematic_success": r_fail,
        },
        "protocol": {
            "repeats": args.repeats,
            "aggregation": "median over repeats of total wall time / rows (perf_counter)",
            "cold_inputs": "inputs rebuilt before every repeat (untimed): chematic re-parses SMILES, RDKit re-parses with MolFromSmiles",
            "order": "per operation: all chematic repeats, then all RDKit repeats",
            "exceptions": "counted per repeat as failures, never timed as wins",
            "ratio": "rdkit_us / chematic_us (>1 means chematic faster)",
        },
        "operations": [],
    }

    def fresh(items):
        if items and isinstance(items[0], chematic.Mol):
            return [chematic.from_smiles(S[i]) for i in range(len(items))]
        if items and isinstance(items[0], Chem.Mol):
            return [Chem.MolFromSmiles(S[i]) for i in range(len(items))]
        return list(items)

    def timeit(fn, items):
        per_repeat = []
        fails = 0
        for _ in range(args.repeats):
            fails = 0
            batch = fresh(items)
            t0 = time.perf_counter()
            for x in batch:
                try:
                    fn(x)
                except Exception:
                    fails += 1
            per_repeat.append((time.perf_counter() - t0) / len(batch) * 1e6)
        return per_repeat, fails

    def op(name, c_fn, r_fn, *, eq, note="", c_items=None, r_items=None, sub=None,
           compare=None, per_item_divisor=1):
        if args.only and args.only not in name:
            return
        ci = c_items if c_items is not None else CM
        ri = r_items if r_items is not None else RM
        if sub:
            ci, ri = ci[:sub], ri[:sub]
        agreement = None
        if compare is not None:
            agree = compared = 0
            for x, y in zip(fresh(ci), fresh(ri)):
                try:
                    a, b = c_fn(x), r_fn(y)
                except Exception:
                    continue
                compared += 1
                agree += bool(compare(a, b))
            agreement = {"compared": compared, "agree": agree}
        c_rep, c_f = timeit(c_fn, ci)
        r_rep, r_f = timeit(r_fn, ri)
        c_rep = [t / per_item_divisor for t in c_rep]
        r_rep = [t / per_item_divisor for t in r_rep]
        ct, rt = statistics.median(c_rep), statistics.median(r_rep)
        ratio = rt / ct if ct > 0 else math.inf
        counted = eq == "checked" and agreement is not None
        row = {
            "op": name, "equivalence": eq, "note": note, "rows": len(ci),
            "chematic_us_per_row": ct, "rdkit_us_per_row": rt, "ratio": ratio,
            "chematic_repeats_us": c_rep, "rdkit_repeats_us": r_rep,
            "chematic_failures": c_f, "rdkit_failures": r_f,
            "output_agreement": agreement, "counted": counted,
        }
        record["operations"].append(row)
        agr = f"{agreement['agree']}/{agreement['compared']}" if agreement else "-"
        print(f"{name:44s} {eq:15s} chematic {ct:10.2f} us  rdkit {rt:10.2f} us  x{ratio:7.2f}  agree {agr}", flush=True)

    tol = lambda t: (lambda a, b: close(float(a), float(b), t))
    eq_int = lambda a, b: int(a) == int(b)

    # --- parsing / writing
    op("parse_smiles", chematic.from_smiles, Chem.MolFromSmiles, eq="same-task",
       note="RDKit sanitizes (rings/aromaticity) at parse; chematic defers perception", c_items=S, r_items=S)
    op("canonical_smiles(prepared)", lambda m: m.smiles, Chem.MolToSmiles, eq="same-task",
       note="different canonical spelling")
    op("parse+canonical_smiles", lambda s: chematic.from_smiles(s).smiles,
       lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s)), eq="same-task", c_items=S, r_items=S)
    op("mol_block_write", lambda m: m.to_mol_block(), Chem.MolToMolBlock, eq="not-equivalent",
       note="RDKit computes 2D coordinates; chematic writes without layout", sub=2000)
    cb = [m.to_mol_block() for m in CM[:2000]]
    rb = [Chem.MolToMolBlock(m) for m in RM[:2000]]
    op("mol_block_read", chematic.from_mol_block, Chem.MolFromMolBlock, eq="same-task",
       note="each library reads its own writer's block", c_items=cb, r_items=rb)
    op("inchi", lambda m: m.inchi, Chem.MolToInchi, eq="not-equivalent", sub=1000, compare=lambda a, b: a == b,
       note="chematic default Mol.inchi is its own non-standard layering, not IUPAC Standard InChI "
            "(standard_inchi needs the native-inchi build); agreement is recorded, never counted")
    op("inchikey", lambda m: m.inchikey, Chem.MolToInchiKey, eq="not-equivalent", sub=1000,
       compare=lambda a, b: a == b, note="see inchi")

    # --- descriptors
    op("mw", lambda m: m.rdkit_mw, Descriptors.MolWt, eq="checked", compare=tol(1e-3))
    op("exact_mass", lambda m: m.exact_mass, Descriptors.ExactMolWt, eq="checked", compare=tol(1e-3))
    op("logp", lambda m: m.logp, Crippen.MolLogP, eq="checked", compare=tol(1e-6))
    op("mr", lambda m: m.molar_refractivity, Crippen.MolMR, eq="checked", compare=tol(1e-2))
    op("tpsa", lambda m: m.tpsa, rdMolDescriptors.CalcTPSA, eq="checked", compare=tol(0.1),
       note="chematic default includes S/P (RDKit includeSandP=False default); agreement shows the gap")
    op("hbd", lambda m: m.hbd, rdMolDescriptors.CalcNumHBD, eq="checked", compare=eq_int)
    op("hba", lambda m: m.rdkit_hba, rdMolDescriptors.CalcNumHBA, eq="checked", compare=eq_int)
    op("rotatable_bonds", lambda m: m.rotatable_bonds, rdMolDescriptors.CalcNumRotatableBonds,
       eq="checked", compare=eq_int)
    op("fsp3", lambda m: m.fsp3, rdMolDescriptors.CalcFractionCSP3, eq="checked", compare=tol(1e-3))
    op("ring_count", lambda m: m.ring_count, rdMolDescriptors.CalcNumRings, eq="checked", compare=eq_int,
       note="RDKit ring info is computed inside MolFromSmiles (untimed here)")
    op("aromatic_ring_count", lambda m: m.aromatic_ring_count, rdMolDescriptors.CalcNumAromaticRings,
       eq="checked", compare=eq_int, note="RDKit aromaticity is computed inside MolFromSmiles (untimed here)")
    op("qed", lambda m: m.qed, QED.qed, eq="checked", compare=tol(1e-3))
    op("kappa2", lambda m: m.kappa2, Descriptors.Kappa2, eq="checked", compare=tol(1e-3))
    op("chi1v", lambda m: m.chi1v, Descriptors.Chi1v, eq="checked", compare=tol(1e-3))
    op("bertz_ct", lambda m: m.bertz_ct, Descriptors.BertzCT, eq="checked", compare=tol(1e-2))
    op("formula", lambda m: m.formula, rdMolDescriptors.CalcMolFormula, eq="checked", compare=lambda a, b: a == b)
    op("num_stereocenters", lambda m: m.num_stereocenters,
       lambda m: len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)),
       eq="checked", compare=eq_int)
    op("lipinski_bundle(mw,logp,hbd,hba)", lambda m: (m.rdkit_mw, m.logp, m.hbd, m.rdkit_hba),
       lambda m: (Descriptors.MolWt(m), Crippen.MolLogP(m), rdMolDescriptors.CalcNumHBD(m), rdMolDescriptors.CalcNumHBA(m)),
       eq="checked", compare=lambda a, b: close(a[0], b[0], 1e-3) and close(a[1], b[1], 1e-6) and a[2] == b[2] and a[3] == b[3])
    op("parse+logp", lambda s: chematic.from_smiles(s).logp, lambda s: Crippen.MolLogP(Chem.MolFromSmiles(s)),
       eq="checked", compare=tol(1e-6), c_items=S, r_items=S)
    op("parse+tpsa", lambda s: chematic.from_smiles(s).tpsa, lambda s: rdMolDescriptors.CalcTPSA(Chem.MolFromSmiles(s)),
       eq="checked", compare=tol(0.1), c_items=S, r_items=S)
    op("parse+ring_count", lambda s: chematic.from_smiles(s).ring_count,
       lambda s: rdMolDescriptors.CalcNumRings(Chem.MolFromSmiles(s)), eq="checked", compare=eq_int, c_items=S, r_items=S)
    op("parse+aromatic_ring_count", lambda s: chematic.from_smiles(s).aromatic_ring_count,
       lambda s: rdMolDescriptors.CalcNumAromaticRings(Chem.MolFromSmiles(s)), eq="checked", compare=eq_int,
       c_items=S, r_items=S)

    # --- fingerprints
    mg = rfg.GetMorganGenerator(radius=2, fpSize=2048)
    op("morgan_r2_2048(native ecfp4)", lambda m: m.ecfp4(), mg.GetFingerprint, eq="not-equivalent",
       note="chematic native hash differs from RDKit Morgan")
    op("morgan_r2_2048(rdkit-compatible)", lambda m: m.rdkit_ecfp4(), mg.GetFingerprint, eq="checked",
       compare=lambda a, b: bytes(a) == bits_of_rdkit(b))
    def maccs_agree(a, b):
        # chematic stores MACCS key k at bit k-1; RDKit at bit k (bit 0 unused).
        c_on = {i + 1 for i in range(166) if a[i // 8] >> (i % 8) & 1}
        return c_on == set(b.GetOnBits()) - {0}

    op("maccs", lambda m: m.maccs(), MACCSkeys.GenMACCSKeys, eq="checked", compare=maccs_agree,
       note="agreement = all 166 keys identical (chematic bit k-1 <-> RDKit bit k)")
    apg = rfg.GetAtomPairGenerator(fpSize=2048)
    op("atom_pair(rdkit-compatible)", lambda m: m.rdkit_atom_pair_fp(), apg.GetFingerprint, eq="checked",
       compare=lambda a, b: bytes(a) == bits_of_rdkit(b))
    ttg = rfg.GetTopologicalTorsionGenerator(fpSize=2048)
    op("torsion(rdkit-compatible)", lambda m: m.rdkit_torsion_fp(), ttg.GetFingerprint, eq="checked",
       compare=lambda a, b: bytes(a) == bits_of_rdkit(b))
    rdkg = rfg.GetRDKitFPGenerator(fpSize=2048)
    op("rdkit_fp(rdkit-compatible)", lambda m: m.rdkit_rdk_fp(), rdkg.GetFingerprint, eq="checked", sub=2000,
       compare=lambda a, b: bytes(a) == bits_of_rdkit(b))
    op("pattern_fp(rdkit-compatible)", lambda m: m.rdkit_pattern_fp(), Chem.PatternFingerprint, eq="checked",
       sub=2000, compare=lambda a, b: bytes(a) == bits_of_rdkit(b))

    cfps = [m.rdkit_ecfp4() for m in CM]
    rfps = [mg.GetFingerprint(m) for m in RM]
    op("tanimoto_1xN(per target, compatible Morgan)", lambda q: chematic.tanimoto_slice(q, cfps),
       lambda q: DataStructs.BulkTanimotoSimilarity(q, rfps), eq="checked",
       compare=lambda a, b: len(a) == len(b) and all(close(x, y, 1e-6) for x, y in zip(a, b)),
       c_items=cfps[:50], r_items=rfps[:50], per_item_divisor=n,
       note=f"50 queries x {n} targets; time reported per target; chematic returns f32")

    # --- substructure
    for p in ["c1ccccc1", "[OH]", "C(=O)N", "[#7;R]", "c1ccc2ccccc2c1", "[CX3](=O)[OX2H1]",
              "[NX3;H2,H1;!$(NC=O)]", "*~*~*~*~*~*"]:
        rp = Chem.MolFromSmarts(p)
        op(f"has_substructure {p}", lambda m, p=p: m.has_substructure(p),
           lambda m, rp=rp: m.HasSubstructMatch(rp), eq="checked", sub=2000, compare=lambda a, b: a == b,
           note="chematic takes the SMARTS string (parse cached); RDKit uses a pre-parsed pattern")
    rp = Chem.MolFromSmarts("[#6]~[#7]")
    op("find_matches [#6]~[#7]", lambda m: m.find_matches("[#6]~[#7]"),
       lambda m: m.GetSubstructMatches(rp, uniquify=True), eq="checked", sub=2000,
       compare=lambda a, b: sorted(tuple(sorted(x)) for x in a) == sorted(tuple(sorted(x)) for x in b))

    # --- structure transforms
    def canon(m):
        """RDKit canonical SMILES of a chematic or RDKit output (None if RDKit can't read it)."""
        if isinstance(m, chematic.Mol):
            m = Chem.MolFromSmiles(m.smiles)
        return Chem.MolToSmiles(m) if m is not None else None

    def same_structure(a, b):
        ca, cb = canon(a), canon(b)
        return ca is not None and ca == cb
    op("murcko_scaffold", lambda m: m.scaffold(), MurckoScaffold.GetScaffoldForMol, eq="checked",
       compare=same_structure, note="agreement via RDKit canonical SMILES of both outputs")
    op("add_hydrogens", lambda m: m.add_hydrogens(), Chem.AddHs, eq="same-task")
    op("remove_hydrogens", lambda m: m.remove_hydrogens(), Chem.RemoveHs, eq="same-task")
    op("largest_fragment", lambda m: m.largest_fragment(), rdMolStandardize.LargestFragmentChooser().choose,
       eq="checked", sub=2000, compare=same_structure,
       note="agreement via RDKit canonical SMILES of both outputs")
    op("neutralize", lambda m: m.neutralize(), rdMolStandardize.Uncharger().uncharge, eq="checked", sub=2000,
       compare=same_structure, note="agreement via RDKit canonical SMILES of both outputs")
    te = rdMolStandardize.TautomerEnumerator()
    op("standardize vs Cleanup", lambda m: m.standardize(), rdMolStandardize.Cleanup, eq="not-equivalent",
       sub=300, note="chematic default also neutralizes and canonicalizes the tautomer")
    op("standardize vs Cleanup+Uncharge+Canonicalize", lambda m: m.standardize(),
       lambda m: te.Canonicalize(rdMolStandardize.Uncharger().uncharge(rdMolStandardize.Cleanup(m))),
       eq="checked", sub=300, compare=same_structure,
       note="closest RDKit pipeline; agreement via RDKit canonical SMILES")
    op("canonical_tautomer", lambda m: m.canonical_tautomer(), te.Canonicalize, eq="checked", sub=300,
       compare=same_structure, note="different tautomer scoring rules expected")
    op("brics_fragments", lambda m: m.brics_fragments(), lambda m: list(BRICS.BRICSDecompose(m)),
       eq="not-equivalent", sub=500, note="chematic returns fragment Mols without RDKit's dummy-labelled SMILES set")
    op("cip_labels", lambda m: m.cip_stereo(), lambda m: rdCIPLabeler.AssignCIPLabels(Chem.Mol(m)),
       eq="not-equivalent", sub=1000, note="chematic default LegacyFast vs RDKit new CIP labeler")
    op("sssr_rings", lambda m: m.ring_membership(), lambda m: Chem.GetSymmSSSR(m), eq="not-equivalent",
       note="SSSR vs symmetrized SSSR")

    # --- depiction / 3D / MCS
    def rd2d(m):
        m = Chem.Mol(m)
        rdDepictor.Compute2DCoords(m)
        return m

    def rsvg(m):
        d = rdMolDraw2D.MolDraw2DSVG(300, 300)
        rdMolDraw2D.PrepareAndDrawMolecule(d, m)
        d.FinishDrawing()
        return d.GetDrawingText()

    op("svg_depiction", lambda m: m.svg(), rsvg, eq="not-equivalent", sub=300, note="different layout and renderer")
    op("2d_layout", lambda m: m.depict_data(), rd2d, eq="not-equivalent", sub=500, note="different layout algorithms")

    def rd3d(m):
        mh = Chem.AddHs(m)
        AllChem.EmbedMolecule(mh, randomSeed=42)
        return mh

    def rdmmff(m):
        mh = rd3d(m)
        AllChem.MMFFOptimizeMolecule(mh)
        return mh

    op("embed_3d", lambda m: m.generate_3d_etkdg(), rd3d, eq="not-equivalent", sub=100,
       note="no geometry-quality check; see the dedicated 3D records")
    op("embed+minimize_mmff94", lambda m: m.minimize_mmff94(m.generate_3d_etkdg()), rdmmff,
       eq="not-equivalent", sub=10, note="legacy numeric-gradient chematic path; no quality check")
    mcs_c = [[CM[i], CM[i + 1]] for i in range(0, 100, 2)]
    mcs_r = [[RM[i], RM[i + 1]] for i in range(0, 100, 2)]
    op("mcs_pair", lambda ms: chematic.find_mcs(ms, timeout_ms=2000), lambda ms: rdFMCS.FindMCS(ms, timeout=2),
       eq="not-equivalent", c_items=mcs_c, r_items=mcs_r, note="different defaults/timeouts; inputs not rebuilt")

    ops = record["operations"]
    counted = [o for o in ops if o["counted"]]
    record["summary"] = {
        "operations": len(ops),
        "counted_operations": len(counted),
        "counted_faster": sum(o["ratio"] > 1.0 for o in counted),
        "counted_faster_with_full_agreement": sum(
            o["ratio"] > 1.0 and o["output_agreement"]["agree"] == o["output_agreement"]["compared"] for o in counted
        ),
        "same_task_operations": sum(o["equivalence"] == "same-task" for o in ops),
        "not_equivalent_operations": sum(o["equivalence"] == "not-equivalent" for o in ops),
    }
    with open(args.output_json, "w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=1)
        fh.write("\n")
    print(json.dumps(record["summary"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
