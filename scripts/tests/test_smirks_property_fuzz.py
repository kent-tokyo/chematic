"""The SMIRKS property fuzz lane's generator (#754); chematic is not needed."""

from __future__ import annotations

import importlib.util
import random
import sys
from pathlib import Path

import pytest

pytest.importorskip("rdkit")
from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
SPEC = importlib.util.spec_from_file_location("smirks_property_fuzz", SCRIPTS / "smirks_property_fuzz.py")
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_generated_templates_parse_and_mostly_match_their_source():
    rng = random.Random(754)
    pool = module.load_pool(module.DEFAULT_REACTANTS)
    made = matched = 0
    while made < 200:
        source = rng.choice(pool)
        template = module.generate(Chem.MolFromSmiles(source), rng)
        if template is None:
            continue
        made += 1
        rxn = AllChem.ReactionFromSmarts(template["smirks"])
        assert rxn.GetNumReactantTemplates() == 1
        query = rxn.GetReactantTemplate(0)
        matched += Chem.MolFromSmiles(source).HasSubstructMatch(query)
    # About one template in seven carries a primitive its atom fails.
    assert matched >= 150


def test_the_generator_is_reproducible():
    pool = module.load_pool(module.DEFAULT_REACTANTS)

    def run(seed):
        rng = random.Random(seed)
        out = []
        while len(out) < 30:
            t = module.generate(Chem.MolFromSmiles(rng.choice(pool)), rng)
            if t is not None:
                out.append(t["smirks"])
        return out

    assert run(1) == run(1)
    assert run(1) != run(2)


def test_ring_closures_are_written_for_cycles():
    mol = Chem.MolFromSmiles("C1CC1")
    labels = {i: f"[#6:{i + 1}]" for i in range(3)}
    bonds = {(0, 1): "-", (1, 2): "-", (0, 2): "-"}
    smarts = module.write_graph([0, 1, 2], labels, bonds)
    assert Chem.MolFromSmarts(smarts).GetNumBonds() == 3
    assert mol.HasSubstructMatch(Chem.MolFromSmarts(smarts))


def test_simplifications_drop_one_term_at_a_time():
    out = module.simplifications("[#6;H1:1]-;!@[#7&D2:2]>>[*:1].[*:2]")
    assert "[#6:1]-;!@[#7&D2:2]>>[*:1].[*:2]" in out
    assert "[#6;H1:1]-[#7&D2:2]>>[*:1].[*:2]" in out
