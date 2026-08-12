"""
Opting into the accurate CIP engine (Milestone 5A/5B).

`Mol.cip_stereo()` (no args, or `mode="legacy"`) is unchanged from every prior
release -- ~96.3% agreement with RDKit's modern rdCIPLabeler oracle, always
returns *an* answer, never "I don't know".

`Mol.cip_stereo(mode="accurate")` uses a hierarchical-digraph engine instead
for tetrahedral R/S (~99.6% oracle-stable agreement -- see
docs/rfcs/cip_accurate_rfc.md), merged with legacy's E/Z and allene answers (the
accurate engine doesn't compute either). It is strictly opt-in: nothing about
your existing code changes unless you pass mode="accurate" explicitly.

The accurate engine is honest about what it doesn't know: a genuine CIP-rule
tie or a computation-budget overrun surfaces via `cip_stereo_unresolved()`
instead of a silently-guessed label. As of this writing that affects ~0.4% of
stereocenters on a 5,000-molecule benchmark (see docs/rfcs/cip_accurate_rfc.md's
Milestone 5B entry) -- almost entirely a known, already-documented family of
symmetric cage/adamantane-like stereocenters (Milestone 4A-2, needs
symmetry/automorphism detection) plus 2 cyclophosphazene rows where no RDKit
oracle has a stable answer either (Milestone 4C-1).
"""

import chematic

CHIRAL_AMINO_ACID = "C[C@H](N)C(=O)O"  # L-alanine
EZ_AND_RS = "C/C=C/[C@H](N)C(=O)O"  # both a double bond and a stereocenter
# A cyclophosphazene where the accurate engine ties (Milestone 4C-1) but legacy
# still produces an answer -- legacy's answer isn't wrong here, it's simply
# unverifiable: no RDKit oracle has a stable label for this molecule either
# (see docs/rfcs/cip_accurate_rfc.md), which is exactly the kind of case the
# accurate engine is designed to flag rather than silently resolve.
TIED_WITH_LEGACY_FALLBACK = (
    "CNP1(NC)=N[P@](NC)(N2CC2)=NP(NC)(NC)=N[P@@](NC)(N2CC2)=N1"
)

# ── 1. Default behavior is unchanged ────────────────────────────────────────
mol = chematic.from_smiles(CHIRAL_AMINO_ACID)
assert mol.cip_stereo() == mol.cip_stereo(mode="legacy")
print(f"legacy (default): {mol.cip_stereo()}")

# ── 2. Opt in per call ───────────────────────────────────────────────────────
print(f"accurate:         {mol.cip_stereo(mode='accurate')}")

# ── 3. Accurate mode merges tetrahedral R/S with legacy E/Z (it computes no
#      bond stereo itself) -- you get both, from one call ──────────────────
mol_ez = chematic.from_smiles(EZ_AND_RS)
stereo = mol_ez.cip_stereo(mode="accurate")
print(f"R/S + E/Z merge:  {stereo}")

# ── 4. Ties are reported, never guessed -- and never silently backfilled
#      with legacy's (less rigorous) answer, even though legacy has one ────
mol_tied = chematic.from_smiles(TIED_WITH_LEGACY_FALLBACK)
print(f"legacy:     {mol_tied.cip_stereo(mode='legacy')}")
print(f"accurate:   {mol_tied.cip_stereo(mode='accurate')}  (empty -- no guess)")
print(f"unresolved: {mol_tied.cip_stereo_unresolved()}")

# ── 5. A drop-in helper for "give me the best answer, tell me if you can't" ─
def cip_stereo_best_effort(mol):
    """Accurate R/S/E/Z where resolvable, falls back to legacy only for atoms
    the accurate engine explicitly couldn't resolve. Marks the fallback so
    callers can distinguish "rigorous" from "best guess" per atom."""
    accurate = {d["atom_idx"]: d["descriptor"] for d in mol.cip_stereo(mode="accurate")}
    legacy = {d["atom_idx"]: d["descriptor"] for d in mol.cip_stereo(mode="legacy")}
    unresolved_idx = {d["atom_idx"] for d in mol.cip_stereo_unresolved()}
    result = []
    for idx, descriptor in accurate.items():
        result.append({"atom_idx": idx, "descriptor": descriptor, "source": "accurate"})
    for idx in unresolved_idx:
        if idx in legacy:
            result.append({"atom_idx": idx, "descriptor": legacy[idx], "source": "legacy_fallback"})
    return result


print(f"\nbest-effort (tied): {cip_stereo_best_effort(mol_tied)}")
