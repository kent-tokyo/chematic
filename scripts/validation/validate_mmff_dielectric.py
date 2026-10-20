#!/usr/bin/env python3
"""Live MMFF dielectric comparisons plus an original Open Babel #2811 probe."""
import argparse
import gzip
import hashlib
import json
import math
import platform
import subprocess
from pathlib import Path

import chematic
import chematic.chematic as extension
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / 'validation/openbabel_issues/mmff-dielectric/cases.json'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def reference_field(mol, variant, epsilon):
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=variant)
    props.SetMMFFDielectricConstant(epsilon)
    props.SetMMFFDielectricModel(1)  # constant dielectric
    field = AllChem.MMFFGetMoleculeForceField(mol, props, nonBondedThresh=100., ignoreInterfragInteractions=False)
    field.Initialize()
    return field


def compare(actual, expected, absolute=1e-8, relative=1e-10):
    assert len(actual) == len(expected)
    errors = []
    for a, b in zip(actual, expected):
        assert math.isfinite(a) and math.isfinite(b)
        error = abs(a-b)
        assert error <= absolute + relative*abs(b), (a,b,error)
        errors.append(error)
    return max(errors, default=0.)


def set_coords(mol, coords):
    for i, xyz in enumerate(coords):
        mol.GetConformer().SetAtomPosition(i, xyz)


def diagnostic_error(actual, expected):
    assert len(actual) == len(expected)
    assert all(math.isfinite(x) for x in actual + expected)
    return max((abs(a-b) for a,b in zip(actual,expected)), default=0.)


def cases(data):
    for row in data['cases']:
        for seed in (7, 19, 71):
            mol = Chem.AddHs(Chem.MolFromSmiles(row['smiles']))
            if '.' in row['smiles']:
                conformer = Chem.Conformer(mol.GetNumAtoms())
                for i in range(mol.GetNumAtoms()):
                    conformer.SetAtomPosition(i, ((2 + seed % 7) * i, .2 * i, -.1 * i))
                mol.AddConformer(conformer)
            else:
                assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0, (row['smiles'], seed)
            coords = mol.GetConformer().GetPositions().tolist()
            yield row['smiles'], seed, mol, coords


def check_goldens(data):
    count = 0
    for row in data['cases']:
        mol = Chem.AddHs(Chem.MolFromSmiles(row['smiles']))
        conformer = Chem.Conformer(mol.GetNumAtoms())
        for i, xyz in enumerate(row['coords']):
            conformer.SetAtomPosition(i, xyz)
        mol.AddConformer(conformer)
        for expected in row['references']:
            field = reference_field(mol, expected['variant'], expected['epsilon'])
            compare([field.CalcEnergy()], [expected['energy']], absolute=1e-12, relative=1e-12)
            compare(list(field.CalcGrad()), expected['gradient'], absolute=1e-12, relative=1e-12)
            count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--openbabel-probe', required=True, type=Path)
    parser.add_argument('--obabel', default='/opt/homebrew/bin/obabel')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    assert rdBase.rdkitVersion == '2026.09.1'
    data = json.loads(FIXTURE.read_text(encoding='utf-8'))
    goldens = check_goldens(data)
    raw = []
    max_energy = max_gradient = max_opt_energy = max_opt_coord = 0.
    evaluations = optimizations = 0
    minimization_residuals = []
    for smiles, seed, reference, coords in cases(data):
        native = chematic.from_smiles(smiles).add_hydrogens()
        assert native.atom_count == reference.GetNumAtoms()
        assert native.formal_charge_per_atom() == [a.GetFormalCharge() for a in reference.GetAtoms()]
        for variant in ('MMFF94', 'MMFF94s'):
            field = native.mmff_force_field(coords, variant, ignore_interfrag_interactions=False)
            initial = field.energy_gradient()
            for epsilon in (1., 4., 78., .5, 4., 1.):
                field.set_dielectric_constant_(epsilon)
                oracle = reference_field(reference, variant, epsilon)
                energy, nested_gradient = field.energy_gradient()
                gradient = [x for p in nested_gradient for x in p]
                want_energy, want_gradient = oracle.CalcEnergy(), list(oracle.CalcGrad())
                max_energy = max(max_energy, compare([energy], [want_energy]))
                max_gradient = max(max_gradient, compare(gradient, want_gradient))
                assert field.positions == coords and field.dielectric_constant == epsilon
                raw.append(dict(lane='energy_gradient', smiles=smiles, seed=seed, variant=variant, epsilon=epsilon, coords=coords, native_energy=energy, rdkit_energy=want_energy, native_gradient=gradient, rdkit_gradient=want_gradient))
                evaluations += 1
            assert field.energy_gradient() == initial
            for epsilon in (1., 4., 78.):
                # The native object is reused after previous evaluations and
                # minimizations. Reset only coordinates and screening. RDKit
                # builds an independent oracle with the same starting points.
                field.set_positions_(coords).set_dielectric_constant_(epsilon)
                set_coords(reference, coords)
                oracle = reference_field(reference, variant, epsilon)
                expected_status = oracle.Minimize(maxIts=10)
                expected_energy = oracle.CalcEnergy()
                status, energy = field.minimize_(10)
                fresh_native = native.mmff_force_field(coords, variant, ignore_interfrag_interactions=False)
                fresh_native.set_dielectric_constant_(epsilon)
                assert fresh_native.minimize_(10) == (status, energy)
                assert fresh_native.positions == field.positions
                actual_positions = field.positions
                expected_positions = reference.GetConformer().GetPositions().tolist()
                energy_error = diagnostic_error([energy], [expected_energy])
                coord_error = diagnostic_error([x for p in actual_positions for x in p], [x for p in expected_positions for x in p])
                max_opt_energy = max(max_opt_energy, energy_error)
                max_opt_coord = max(max_opt_coord, coord_error)
                if status != expected_status or energy_error > 1e-6 + 1e-9*abs(expected_energy) or coord_error > 5e-6:
                    minimization_residuals.append(dict(smiles=smiles, seed=seed, variant=variant, epsilon=epsilon, native_status=status, rdkit_status=expected_status, energy_abs_error=energy_error, coordinate_max_abs_error=coord_error))
                raw.append(dict(lane='minimize_10_steps', smiles=smiles, seed=seed, variant=variant, epsilon=epsilon, start=coords, native_status=status, rdkit_status=expected_status, native_energy=energy, rdkit_energy=expected_energy, native_coords=actual_positions, rdkit_coords=expected_positions))
                optimizations += 1
            set_coords(reference, coords)
    probe = subprocess.run([str(args.openbabel_probe.resolve())], text=True, capture_output=True, check=True)
    ob_rows = json.loads(probe.stdout)
    for mode in ('live', 'setup_again'):
        rows = [r for r in ob_rows if r['mode'] == mode]
        assert len(rows) == 4 and len({r['electrostatic'] for r in rows}) == 1
        assert all(r['epsilon'] == r['reported_epsilon'] for r in rows)
    fresh = [r for r in ob_rows if r['mode'] == 'fresh']
    assert len(fresh) == 4
    for row in fresh:
        compare([row['electrostatic']], [fresh[0]['electrostatic'] / row['epsilon']], absolute=1e-12)
    args.output.mkdir(parents=True, exist_ok=True)
    packet = gzip.compress(json.dumps({'comparisons':raw, 'openbabel':ob_rows, 'openbabel_stdout':probe.stdout, 'openbabel_stderr':probe.stderr},sort_keys=True).encode(),mtime=0)
    (args.output/'comparisons.json.gz').write_bytes(packet)
    sources = ['crates/chematic-ff/src/rdkit_mmff.rs', 'crates/chematic-ff/src/rdkit_mmff/dielectric.rs', 'crates/chematic-py/src/force_field.rs', 'scripts/validation/validate_mmff_dielectric.py', 'validation/openbabel_issues/mmff-dielectric/cases.json', 'validation/openbabel_issues/mmff-dielectric/openbabel_probe.cpp']
    result = dict(schema='chematic.mmff-dielectric-validation.v1', source_base=subprocess.check_output(['git','rev-parse','HEAD'],text=True,cwd=ROOT).strip(), source_sha256={s:sha((ROOT/s).read_bytes()) for s in sources}, runtime=dict(rdkit=rdBase.rdkitVersion, python=platform.python_version(), platform=platform.platform(), python_extension_sha256=sha(Path(extension.__file__).read_bytes()), openbabel=subprocess.check_output([args.obabel,'-V'],text=True).strip(), openbabel_probe_sha256=sha(args.openbabel_probe.read_bytes())), results=dict(golden_energy_gradient_checks=goldens, live_energy_gradient_checks=evaluations, minimization_fresh_native_matches=optimizations, minimization_rdkit_diagnostics=optimizations, minimization_rdkit_status_matches=optimizations-sum(r['native_status'] != r['rdkit_status'] for r in minimization_residuals), minimization_diagnostic_tolerances=dict(energy_absolute=1e-6, energy_relative=1e-9, coordinate_absolute=5e-6), minimization_rdkit_residuals=minimization_residuals, max_energy_abs_error=max_energy, max_gradient_abs_error=max_gradient, max_minimization_energy_abs_error=max_opt_energy, max_minimization_coordinate_abs_error=max_opt_coord, energy_gradient_failures=0, fresh_native_minimization_failures=0), openbabel_issue=dict(url='https://github.com/openbabel/openbabel/issues/2811', observed='live updates and Setup on the same object keep stale energy; a fresh instance respects epsilon', fixture='authored Na+/Cl- pair at 4 Angstrom'), packet=dict(path='comparisons.json.gz',sha256=sha(packet)), scope='local MMFF94/MMFF94s constant screening; minimizer oracle residuals are reported separately; no general solvation, distance-dependent dielectric, UFF/GAFF, speed or package release claim')
    (args.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
