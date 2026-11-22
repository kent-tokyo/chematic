#!/usr/bin/env python3
"""Run the strict development probe with the checked mixed-gradient profile."""
import json
import sys
from pathlib import Path

import platinum_development_probe as probe
from platinum_stereo_strict import annotate_result, stereo_violations_strict


class MixedMolecule:
    def __init__(self, molecule):
        self.molecule = molecule

    def add_hydrogens(self):
        return MixedMolecule(self.molecule.add_hydrogens())

    def __getattr__(self, name):
        if name == "etkdg_embed_multiple_confs_mmff_corrected":
            return self.molecule.etkdg_embed_multiple_confs_mmff_mixed_checked
        return getattr(self.molecule, name)


if __name__ == "__main__":
    original_from_smiles = probe.chematic.from_smiles
    probe.chematic.from_smiles = lambda smiles: MixedMolecule(original_from_smiles(smiles))
    probe.stereo_violations = stereo_violations_strict
    probe.main()
    output = Path(sys.argv[sys.argv.index("--output") + 1])
    annotate_result(output, __file__)
    result = json.loads(output.read_text())
    result["profile"] = "mixed_checked"
    temporary = output.with_suffix(output.suffix + ".profile.tmp")
    temporary.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    temporary.replace(output)
