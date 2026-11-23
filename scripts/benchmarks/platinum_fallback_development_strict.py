#!/usr/bin/env python3
"""Evaluate the fallback profile only on exposed development records."""
import sys
import platinum_fallback_development_probe as probe
from platinum_stereo_strict import annotate_result, stereo_violations_strict

if __name__ == '__main__':
    probe.stereo_violations = stereo_violations_strict
    probe.main()
    annotate_result(sys.argv[sys.argv.index('--output') + 1], __file__)
