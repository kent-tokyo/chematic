#!/usr/bin/env python3
"""Time one-attempt mixed-gradient ensembles with strict stereo checks."""
import json
import sys
from pathlib import Path

import platinum_time_budget as probe
from platinum_stereo_strict import annotate_result, stereo_violations_strict


class MixedAttempts:
    def __init__(self, molecule, compatible_attempts):
        self.molecule = molecule
        self.compatible_attempts = compatible_attempts
        self.attempt = 0

    def etkdg_embed_multiple_confs_mmff_corrected(self, *args, **kwargs):
        if self.attempt < self.compatible_attempts:
            method = self.molecule.etkdg_embed_multiple_confs_mmff_compatible_checked
        else:
            method = self.molecule.etkdg_embed_multiple_confs_mmff_corrected
        self.attempt += 1
        return method(*args, **kwargs)


if __name__ == '__main__':
    original_attempts = probe.attempts

    def mixed_attempts(engine, mol, other, seconds, max_candidates, repair_iters=0):
        if engine == 'chematic':
            other = MixedAttempts(other, max_candidates // 2)
        return original_attempts(engine, mol, other, seconds, max_candidates, repair_iters)

    probe.attempts = mixed_attempts
    probe.stereo_violations = stereo_violations_strict
    probe.main()
    output = Path(sys.argv[sys.argv.index('--output') + 1])
    annotate_result(output, __file__)
    result = json.loads(output.read_text())
    result['profile'] = 'mixed_checked'
    result['gradient_schedule'] = 'first floor(max_candidates / 2) compatible, remainder corrected'
    temporary = output.with_suffix(output.suffix + '.profile.tmp')
    temporary.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    temporary.replace(output)
