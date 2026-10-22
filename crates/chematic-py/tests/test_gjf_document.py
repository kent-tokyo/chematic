"""Explicit all-job retention; no implicit first/last Cartesian job selection."""
import json
import struct
from pathlib import Path

import pytest
import chematic

ROOT = Path(__file__).parents[3]
CASES = json.loads((ROOT / 'validation/gjf-document-2026-10-22/cases.json').read_text(encoding="utf-8"))['cases']
SOURCE = (ROOT / 'validation/gjf-document-2026-10-22/two-jobs.gjf').read_bytes().decode()


def assert_jobs(source, actual, expected):
    assert len(actual) == len(expected)
    raw = source.encode('utf8')
    for index, (job, control) in enumerate(zip(actual, expected)):
        assert job['index'] == index
        for key in ['source_start', 'source_end', 'elements', 'freeze_codes', 'charge', 'multiplicity']:
            assert job[key] == control[key]
        for got, wanted in zip(job['coords'][0], control['coords'][0]):
            assert struct.pack('>d', got) == struct.pack('>d', wanted)
        piece = raw[job['source_start']:job['source_end']].decode('utf8')
        assert '--Link1--' not in piece


@pytest.mark.parametrize('case', CASES, ids=lambda c: c['id'])
def test_all_jobs_or_whole_document_refusal(case):
    source = case['source']
    if case['expected'] == 'refuse':
        for call in [chematic.roundtrip_gjf_document, chematic.parse_gjf_document_json]:
            with pytest.raises(ValueError):
                call(source)
        return
    assert chematic.roundtrip_gjf_document(source).encode() == source.encode()
    report = json.loads(chematic.parse_gjf_document_json(source))
    assert report['schema'] == 'chematic.gjf-document.v1'
    assert_jobs(source, report['jobs'], case['jobs'])


@pytest.mark.parametrize('limits', [dict(max_jobs=0), dict(max_jobs=1), dict(max_atoms=0), dict(max_atoms=1), dict(max_input_bytes=len(SOURCE.encode())-1)])
def test_cumulative_limits(limits):
    for call in [chematic.roundtrip_gjf_document, chematic.parse_gjf_document_json]:
        with pytest.raises(ValueError):
            call(SOURCE, **limits)


def test_exact_limits_and_bom_budget():
    source = '\ufeff'+SOURCE
    assert chematic.roundtrip_gjf_document(source, max_jobs=2, max_atoms=2, max_input_bytes=len(source.encode())) == source
    with pytest.raises(ValueError):
        chematic.roundtrip_gjf_document(source, max_input_bytes=len(source.encode())-1)


def test_single_job_routes_keep_explicit_refusal():
    for call in [chematic.roundtrip_gjf_record, chematic.parse_gjf_geometry,
                 lambda s: chematic.convert_format(s, 'gjf', 'com')]:
        with pytest.raises(ValueError):
            call(SOURCE)


def test_job_count_boundary_and_late_invalid_job():
    one = '# HF/STO-3G\n\nHe\n\n0 1\nHe 0 0 0\n\n'
    source = '--Link1--\n'.join([one]*128)
    assert chematic.roundtrip_gjf_document(source) == source
    assert len(json.loads(chematic.parse_gjf_document_json(source))['jobs']) == 128
    for tail in [one, 'invalid']:
        with pytest.raises(ValueError):
            chematic.roundtrip_gjf_document(source+'--Link1--\n'+tail)


def test_atom_budget_counts_across_jobs():
    prefix = '# HF/STO-3G\n\nHe\n\n0 1\n'
    half = prefix+'He 0 0 0\n'*50_000+'\n'
    source = half+'--Link1--\n'+half
    assert chematic.roundtrip_gjf_document(source) == source
    with pytest.raises(ValueError):
        chematic.roundtrip_gjf_document(source.rstrip('\n')+'\nHe 0 0 0\n\n')


@pytest.mark.parametrize('case', json.loads((ROOT / 'validation/gjf-cartesian-2026-10-22/source-cases.json').read_text(encoding="utf-8"))['cases'], ids=lambda c: c['id'])
def test_existing_source_packet_under_explicit_document_profile(case):
    supported = case['expected'] == 'retain' or case['id'].startswith('multiple-jobs-')
    if supported:
        assert chematic.roundtrip_gjf_document(case['source']) == case['source']
    else:
        with pytest.raises(ValueError):
            chematic.roundtrip_gjf_document(case['source'])
