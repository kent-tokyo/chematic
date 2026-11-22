"""Conservative 3D checks for every RDKit-recognized declared stereo element.

The existing recovery scorer consumes ``atoms`` and ``bonds`` as failures.
Unevaluable elements are included there so they cannot count as recovery,
while separate fields distinguish uncertainty from a contradictory geometry.
"""
import math

from rdkit import Chem


def _raw_tetrahedral_parity(mol, index, coords, ranks):
    """Compare a four-neighbor raw tag without requiring CIP assignment."""
    atom = mol.GetAtomWithIdx(index)
    tetra = (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
    neighbors = [a.GetIdx() for a in atom.GetNeighbors()]
    if (atom.GetChiralTag() not in tetra or len(neighbors) != 4
            or len({ranks[i] for i in neighbors}) != 4):
        return None
    p = [coords[i] for i in neighbors]
    vectors = [[p[i][k] - p[3][k] for k in range(3)] for i in range(3)]
    a, b, c = vectors
    volume = (a[0] * (b[1] * c[2] - b[2] * c[1])
              - a[1] * (b[0] * c[2] - b[2] * c[0])
              + a[2] * (b[0] * c[1] - b[1] * c[0]))
    scale = math.prod(math.sqrt(sum(v * v for v in vector)) for vector in vectors)
    if scale == 0.0 or abs(volume) <= 1e-6 * scale:
        return None
    return (volume > 0.0) == (atom.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW)


def stereo_violations_strict(mol, conformers):
    declared_atoms = {a.GetIdx() for a in mol.GetAtoms()
                      if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED}
    declared_bonds = {b.GetIdx() for b in mol.GetBonds()
                      if b.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ,
                                           Chem.BondStereo.STEREOCIS, Chem.BondStereo.STEREOTRANS)}
    expected = Chem.Mol(mol)
    Chem.AssignStereochemistry(expected, cleanIt=True, force=True)
    constraints = [s for s in Chem.FindPotentialStereo(expected)
                   if s.specified == Chem.StereoSpecified.Specified]
    recognized_atoms = {s.centeredOn for s in constraints if str(s.type).startswith('Atom_')}
    recognized_bonds = {s.centeredOn for s in constraints if str(s.type).startswith('Bond_')}
    ranks = list(Chem.CanonicalRankAtoms(expected, breakTies=False, includeChirality=False))
    raw_tetra_ranks = (list(Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=True))
                       if declared_atoms - recognized_atoms else None)
    results = []
    for coords in conformers:
        if len(coords) != mol.GetNumAtoms() or any(
                len(point) != 3 or not all(math.isfinite(value) for value in point)
                for point in coords):
            raise ValueError('nonfinite or invalid conformer coordinates')
        probe = Chem.Mol(mol)
        probe.RemoveAllConformers()
        Chem.RemoveStereochemistry(probe)
        for atom in probe.GetAtoms():
            if atom.HasProp('_chiralPermutation'):
                atom.ClearProp('_chiralPermutation')
        conf = Chem.Conformer(probe.GetNumAtoms())
        for i, point in enumerate(coords):
            conf.SetAtomPosition(i, point)
        probe.AddConformer(conf)
        Chem.AssignStereochemistryFrom3D(probe, replaceExistingTags=True)
        Chem.AssignStereochemistry(probe, cleanIt=True, force=True)
        violated_atoms, unknown_atoms = set(), declared_atoms - recognized_atoms
        violated_bonds, unknown_bonds = set(), declared_bonds - recognized_bonds
        for index in tuple(unknown_atoms):
            status = _raw_tetrahedral_parity(mol, index, coords, raw_tetra_ranks)
            if status is not None:
                unknown_atoms.remove(index)
                if not status:
                    violated_atoms.add(index)
        for stereo in constraints:
            kind, index = str(stereo.type), stereo.centeredOn
            if kind == 'Atom_Tetrahedral':
                want = expected.GetAtomWithIdx(index).GetChiralTag()
                actual = probe.GetAtomWithIdx(index).GetChiralTag()
                tetra = (Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                         Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
                if want not in tetra or actual not in tetra:
                    unknown_atoms.add(index)
                elif want != actual:
                    violated_atoms.add(index)
            elif kind in ('Atom_SquarePlanar', 'Atom_TrigonalBipyramidal', 'Atom_Octahedral'):
                atom = expected.GetAtomWithIdx(index)
                actual = probe.GetAtomWithIdx(index)
                neighbor_ranks = [ranks[a.GetIdx()] for a in atom.GetNeighbors()]
                # Equal-ranked ligands may allow equivalent permutations. Until
                # equivalence is resolved, do not guess whether they agree.
                if (len(set(neighbor_ranks)) != len(neighbor_ranks)
                        or not atom.HasProp('_chiralPermutation')
                        or actual.GetChiralTag() != atom.GetChiralTag()
                        or not actual.HasProp('_chiralPermutation')):
                    unknown_atoms.add(index)
                elif atom.GetIntProp('_chiralPermutation') != actual.GetIntProp('_chiralPermutation'):
                    violated_atoms.add(index)
            elif kind == 'Bond_Double':
                want = expected.GetBondWithIdx(index).GetStereo()
                actual = probe.GetBondWithIdx(index).GetStereo()
                ez = (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ)
                if want not in ez or actual not in ez:
                    unknown_bonds.add(index)
                elif want != actual:
                    violated_bonds.add(index)
            elif kind.startswith('Atom_'):
                unknown_atoms.add(index)
            else:
                unknown_bonds.add(index)
        results.append({
            'atoms': sorted(violated_atoms | unknown_atoms),
            'bonds': sorted(violated_bonds | unknown_bonds),
            'violated_atoms': sorted(violated_atoms),
            'violated_bonds': sorted(violated_bonds),
            'unevaluable_atoms': sorted(unknown_atoms),
            'unevaluable_bonds': sorted(unknown_bonds),
        })
    return results


def annotate_result(path, wrapper_path):
    """Mark only a complete run as having used this strict verifier."""
    import hashlib
    import json
    from pathlib import Path

    path = Path(path)
    result = json.loads(path.read_text())
    result['stereo_policy'] = 'rdkit-declared-stereo-strict-v1'
    result['stereo_verifier_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['stereo_wrapper_sha256'] = hashlib.sha256(Path(wrapper_path).read_bytes()).hexdigest()
    temporary = path.with_suffix(path.suffix + '.stereo.tmp')
    temporary.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)
