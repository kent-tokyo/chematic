"""Independently record RMSD, declared stereo and geometry for each candidate."""
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from platinum_development_probe import rms_values
from platinum_geometry import check_geometry
from platinum_stereo_strict import stereo_violations_strict


def annotate(out, mol, reference, checks):
    candidates = out.pop('candidates')
    coords = [point for _, point in candidates]
    out['statuses'] = [status for status, _ in candidates]
    out['accepted'] = len(coords)
    out['rmsd'] = []
    out['rmsd_errors'] = []
    for index, points in enumerate(coords):
        try:
            rms = rms_values(mol, reference, [points])[0]
        except RuntimeError:
            try:
                probe = Chem.Mol(mol)
                probe.RemoveAllConformers()
                conf = Chem.Conformer(probe.GetNumAtoms())
                for i, point in enumerate(points):
                    conf.SetAtomPosition(i, point)
                probe.AddConformer(conf)
                heavy = Chem.RemoveHs(probe)
                rms = rdMolAlign.GetBestRMS(Chem.Mol(reference), heavy, maxMatches=100000)
            except Exception as exc:
                rms = None
                out['rmsd_errors'].append({'candidate': index, 'error': repr(exc)})
        except Exception as exc:
            rms = None
            out['rmsd_errors'].append({'candidate': index, 'error': repr(exc)})
        out['rmsd'].append(rms)
    out['stereo'] = stereo_violations_strict(mol, coords)
    out['geometry'] = [check_geometry(c, checks) for c in coords]
    return out
