#!/usr/bin/env python3
"""Run the strict development probe with the checked compatibility gradient."""
import sys

import platinum_development_probe as probe
from platinum_stereo_strict import annotate_result, stereo_violations_strict


class CompatibleMolecule:
    def __init__(self, molecule):
        self.molecule = molecule

    def add_hydrogens(self):
        return CompatibleMolecule(self.molecule.add_hydrogens())

    def __getattr__(self, name):
        if name == "etkdg_embed_multiple_confs_mmff_corrected":
            return self.molecule.etkdg_embed_multiple_confs_mmff_compatible_checked
        return getattr(self.molecule, name)


if __name__ == "__main__":
    original_from_smiles = probe.chematic.from_smiles
    probe.chematic.from_smiles = lambda smiles: CompatibleMolecule(original_from_smiles(smiles))
    probe.stereo_violations = stereo_violations_strict
    probe.main()
    annotate_result(sys.argv[sys.argv.index("--output") + 1], __file__)
