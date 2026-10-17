"""Run the published Parent examples against the real extension."""

from pathlib import Path
import re


def test_parent_guide_examples():
    guide = Path(__file__).parents[3] / "docs/standardization-parents.md"
    blocks = re.findall(r"```python\n(.*?)\n```", guide.read_text(), re.S)
    assert blocks, "the published Parent guide must contain executable examples"
    namespace = {"__name__": "parent_api_examples"}
    for index, source in enumerate(blocks):
        exec(compile(source, f"{guide}:python-block-{index + 1}", "exec"), namespace)
