#!/usr/bin/env python3
"""Run the existing fixed-count scorer with conservative declared-stereo checks."""
import sys

import platinum_development_probe as probe
from platinum_stereo_strict import annotate_result, stereo_violations_strict


if __name__ == '__main__':
    probe.stereo_violations = stereo_violations_strict
    probe.main()
    annotate_result(sys.argv[sys.argv.index('--output') + 1], __file__)
