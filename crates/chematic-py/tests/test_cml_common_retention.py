"""Normal CML self-conversion retains source instead of projecting a graph."""
import json
from pathlib import Path
import pytest
import chematic

PACKET = json.loads((Path(__file__).resolve().parents[3]/'validation/cml_source_preservation/cases.json').read_text(encoding="utf-8"))['cases']
NS = 'http://www.xml-cml.org/schema'
OPAQUE = [f'<molecule xmlns="{NS}"/>',
          f'<molecule xmlns="{NS}"><bondArray><bond atomRefs2="missing other" order="uninterpreted"/></bondArray></molecule>',
          f'<molecule xmlns="{NS}" xmlns:v="urn:vendor"><v:molecule stereo="opaque"/></molecule>',
          f'<molecule xmlns="{NS}"><atomArray><atom id="a" elementType="C" x3="uninterpreted"/></atomArray></molecule>']
INVALID = ['<molecule/>', f'<cml xmlns="{NS}"><molecule/></cml>',
           f'<molecule xmlns="{NS}"><molecule/></molecule>',
           f'<!DOCTYPE molecule><molecule xmlns="{NS}"/>',
           f'<?xml version="1.1"?><molecule xmlns="{NS}"/>',
           f'<?xml version="1.0" encoding="ISO-8859-1"?><molecule xmlns="{NS}"/>',
           f'<molecule xmlns="{NS}"><atom id="a" id="b"/></molecule>',
           f'<molecule xmlns="{NS}"><atom/></molecule><bad/>',
           f'<molecule xmlns="{NS}">&#0;</molecule>']


@pytest.mark.parametrize('source', [c['source'] for c in PACKET] + OPAQUE)
def test_direct_and_common_return_exact_source(source):
    for value in [source, source.encode('utf8')]:
        assert chematic.roundtrip_cml_source(value) == source
    assert chematic.convert_format(source, '.CML', 'cml') == source
    # The existing common wrapper accepts str only; the direct API also takes bytes.
    with pytest.raises(ValueError, match='text must be a string'):
        chematic.convert_format(source.encode('utf8'), 'cml', 'cml')
    report = chematic.convert_with_report(source, 'cml', 'cml', policy='preserve_source')
    assert report['output'] == source
    assert report['candidate_graph_verified'] is False


@pytest.mark.parametrize('source', INVALID)
def test_entire_xml_document_is_checked_before_output(source):
    with pytest.raises(ValueError):
        chematic.roundtrip_cml_source(source)
    with pytest.raises(ValueError):
        chematic.convert_format(source, 'cml', 'cml')


@pytest.mark.parametrize('override', [dict(coords=[]), dict(charges=[]), dict(name='new'), dict(comment='new')])
def test_source_retention_never_silently_accepts_an_edit(override):
    with pytest.raises(ValueError, match='does not accept overrides'):
        chematic.convert_format(PACKET[0]['source'], 'cml', 'cml', **override)


def test_utf8_and_byte_limits():
    with pytest.raises(UnicodeDecodeError):
        chematic.roundtrip_cml_source(b'\xff')
    for call in [chematic.roundtrip_cml_source,
                 lambda text: chematic.convert_format(text, 'cml', 'cml')]:
        with pytest.raises(ValueError):
            call(' ' * (16 * 1024 * 1024 + 1))
