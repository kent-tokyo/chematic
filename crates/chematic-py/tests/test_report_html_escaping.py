"""Untrusted report labels must remain text in Jupyter and saved HTML."""

import chematic


def test_report_escapes_names_and_title_in_saved_html(tmp_path):
    name = '<img src=x onerror=alert(1)>'
    title = '</title><script>alert(1)</script>'
    path = tmp_path / 'report.html'

    report = chematic.report(
        [chematic.from_smiles('CCO')], names=[name], title=title, output=path
    )
    html = report._repr_html_()

    assert path.read_text(encoding='utf-8') == html
    assert '<img' not in html and '<script' not in html
    assert '&lt;img src=x onerror=alert(1)&gt;' in html
    assert html.count('&lt;/title&gt;&lt;script&gt;alert(1)&lt;/script&gt;') == 2


def test_compare_escapes_names_and_explicit_or_derived_title():
    mol = chematic.from_smiles('CCO')
    name = '<img src=x onerror=alert(1)>'
    title = '</title><script>alert(1)</script>'

    for report in (
        chematic.compare(mol, mol, names=(name, 'Other')),
        chematic.compare(mol, mol, names=(name, 'Other'), title=title),
    ):
        html = report._repr_html_()
        assert '<img' not in html and '<script' not in html
        assert '&lt;img src=x onerror=alert(1)&gt;' in html
    assert '&lt;/title&gt;&lt;script&gt;alert(1)&lt;/script&gt;' in html
