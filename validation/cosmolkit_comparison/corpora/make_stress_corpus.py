#!/usr/bin/env python3
"""Build ``stress-v1.smi``: odd but RDKit-valid SMILES for the corpus harness.

Every row is a SMILES that RDKit 2026.03.1's ``Chem.MolFromSmiles`` accepts.
Rows come from hand-written families (antiaromatic rings and annulenes,
charged and zwitterionic aromatics, radicals, isotopes, explicit ``[H]``,
hypervalent main-group atoms, metals and dative bonds, non-tetrahedral
stereo, allenes and cumulenes, ring E/Z, macrocycles, large fused PAHs,
bridgehead and quaternary stereo, atom maps, tiny molecules, dummy atoms)
and from seeded random mutations of ChEMBL corpus rows (element, charge,
isotope, radical, bond-order, H-count and chirality edits, atom insertion
and deletion, random atom orders). The output is deterministic for a given
RDKit version and seed.

    python make_stress_corpus.py --chembl ../../../scripts/chembl_accuracy_corpus_4999.smi \
        --output stress-v1.smi
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase

RDLogger.DisableLog("rdApp.*")

SEED = 20261009


def annulenes() -> list[str]:
    out = []
    for n in range(3, 31):
        out.append("c1" + "c" * (n - 1) + "1")
        if n % 2 == 0:
            out.append("C1=" + "CC=" * (n // 2 - 1) + "C1")
        out.append("C1" + "C=C" * (n // 2) + ("C" if n % 2 else "") + "1")
    for n in range(5, 20):
        out.append("[cH-]1" + "c" * (n - 1) + "1")
        out.append("[cH+]1" + "c" * (n - 1) + "1")
        out.append("[nH]1" + "c" * (n - 1) + "1")
        out.append("o1" + "c" * (n - 1) + "1")
        out.append("s1" + "c" * (n - 1) + "1")
        out.append("[n-]1" + "c" * (n - 1) + "1")
        out.append("[o+]1" + "c" * (n - 1) + "1")
        out.append("b1" + "c" * (n - 1) + "1")
        out.append("[se]1" + "c" * (n - 1) + "1")
    out += [
        "C1=CC=C1", "C1=CC=CC=CC=C1", "c1ccc2ccc2c1", "c1cc2cccc2c1", "C1=CC2=CC=CC2=C1",
        "c1ccc2c(c1)C=C2", "c1ccc2c(c1)-c1ccccc1-2", "c1cc2ccc3cccc4ccc(c1)c2c34",
        "C1=CC=CC=CC=CC=C1", "c1cc2cc3cccc3cc2c1", "O=C1C=CC1=O", "O=c1ccc1=O",
        "c1ccc1C", "Cc1ccc1C", "c1cc1", "C1=C1", "[CH+]1C=C1", "c1cc[cH-]1", "C1=CC=CC1",
        "c1ccccc1c1ccc1", "B1C=CC=C1", "c1cc[b-]cc1", "c1ccc[c-]1", "[c-]1cccc1",
        "c1=cc=cc=c1", "c1ccc2c(c1)c1ccccc1c1ccccc21", "C1=CC=C2C=CC=CC=C2C=C1",
    ]
    return out


def charged_aromatics() -> list[str]:
    return [
        "c1cc[n+](C)cc1", "c1cc[n+]([O-])cc1", "[O-][n+]1ccccc1", "C[n+]1ccn(C)c1",
        "Cn1cc[n+](C)c1", "c1ccc2[n+](c1)cc[n-]2", "[NH3+]c1ccccc1C(=O)[O-]",
        "O=C([O-])c1cccc[n+]1C", "c1cc[o+]cc1", "c1cc[s+]cc1", "C[o+]1ccccc1", "[S-]c1ccccc1",
        "[O-]c1ccccc1", "[O-]c1cc[n+](C)cc1", "C[N+](C)(C)c1ccccc1", "[NH+]1=CC=CC=C1",
        "c1cc[nH+]cc1", "[nH+]1ccccc1", "c1c[nH+]c[nH]1", "c1c[nH]c[nH+]1",
        "c1ccc2[nH+]c3ccccc3cc2c1", "C[n+]1c2ccccc2c(N)c2ccccc21", "Cc1cc(C)[o+]c(C)c1",
        "[O-][N+](=O)c1ccc([O-])cc1", "c1ccc(-[n+]2ccccc2)cc1", "[n-]1cccc1", "[n-]1nnnc1C",
        "c1ccc2c(c1)[n-]c1ccccc12", "C[n+]1ccccc1-c1cccc[n+]1C", "[Cl-].c1cc[n+](C)cc1",
        "c1cc[n+](cc1)CC(=O)[O-]", "C[n+]1c(-c2ccccc2)[nH]c2ccccc21", "c1ccc[cH-]1.[Fe+2].c1ccc[cH-]1",
        "[cH-]1cccc1.[cH-]1cccc1.[Fe+2]", "c1cc[c-](c1)C", "[O-]c1cccc[n+]1[O-]",
        "C[N+]1=CC=CC=C1", "C[N+]1=CC=C(C=C1)[O-]", "[NH2+]=c1ccccc1", "O=c1cc[nH+]cc1",
        "c1nc[n-]n1", "[O-]S(=O)(=O)c1ccc(cc1)[n+]1ccccc1", "Cc1c[n+](C)c[nH]1",
        "C[n+]1cccc2ccccc21", "[N-]=[N+]=Nc1ccccc1", "[N-]=[N+]=[N-]", "C=[N+]=[N-]",
        "[C-]#[N+]c1ccccc1", "[C-]#[O+]", "[NH3+]CC(=O)[O-]", "C[N+](C)(C)CC(=O)[O-]",
        "[NH3+][C@@H](Cc1c[nH]c[nH+]1)C(=O)[O-]", "O=C([O-])C[n+]1ccccc1",
        "C[S+](C)c1ccccc1", "C[Se+](C)C", "c1ccc([I+]c2ccccc2)cc1", "[O-][Cl+3]([O-])([O-])[O-]",
        "[O-][I+2]([O-])c1ccccc1",
    ]


def radicals() -> list[str]:
    return [
        "[CH3]", "[CH2]", "[CH]", "[C]", "[NH2]", "[NH]", "[N]", "[OH]", "[O]", "[O][O]",
        "C[CH2]", "C[CH]C", "CC(C)(C)[O]", "[O]N1C(C)(C)CCCC1(C)C", "c1cc[c]cc1",
        "[CH2]c1ccccc1", "C[N]C", "O=[N]", "[O][N+](=O)[O-]", "[N]=O", "C[C](C)C",
        "[CH2]C=C", "C=C[CH2]", "[c]1ccccc1", "[n]1cccc1", "[Si]", "[SiH3]", "[B]", "[BH2]",
        "[S]", "[SH]", "C[S]", "[P]", "[PH2]", "Cl[C](Cl)Cl", "[CH2][CH2]", "[O-][O]",
        "c1ccc2c(c1)[N]c1ccccc12", "[C]=O", "[C-]#[O+]", "[H]", "[H][H]", "[He]",
        "O=C1C=CC(=O)C=C1[O]", "[O]c1ccccc1", "C1=C[CH]C=C1", "[CH]1C=CC=C1",
        "C[C]=O", "[CH]=O", "N#[C]", "[C]#N", "CC[N]C(C)=O", "[O]S(=O)(=O)C",
    ]


def isotopes() -> list[str]:
    out = [
        "[2H]C([2H])([2H])[2H]", "[13CH4]", "[2H]O[2H]", "[3H]C", "[14C]", "[18OH2]",
        "C[13C](=O)O", "[2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H]", "[15NH3]", "[13cH]1ccccc1",
        "[11CH3]Oc1ccccc1", "[18F]c1ccccc1", "[125I]c1ccccc1", "[2H][C@](C)(O)F",
        "[2H]C(C)(O)F", "[13CH3][C@H](N)C(=O)O", "C[C@H]([2H])O", "[2H]/C=C/[2H]",
        "[1H]C", "[1H][1H]", "[0CH4]", "[2H+]", "[3H]O[3H]", "[2H]N([2H])[2H]",
        "[99Tc]", "[235U]", "[12C]", "[2H]C1=CC=CC=C1", "OC[13C@H]1O[13C@@H](O)[13C@H](O)[13C@@H](O)[13C@@H]1O",
        "[2H]C([2H])([2H])C(=O)O", "[19F]C(F)(F)F", "[32P](=O)(O)(O)O", "[35S]C",
        "c1cc[15n]cc1", "[15n]1cccc1", "[2H]n1cccc1", "[2H]c1cc[nH]c1",
    ]
    return out


def explicit_h() -> list[str]:
    return [
        "[H]C([H])([H])[H]", "[H]O[H]", "[H]N([H])[H]", "[H][C@@](C)(N)C(=O)O",
        "[H][C@](F)(Cl)Br", "C([H])([H])([H])C", "[H]c1ccccc1", "[H]/C(C)=C(/[H])C",
        "[H]C(=O)O", "[H][H]", "[H+]", "[H-]", "[H]Cl", "[H][N+]([H])([H])[H]",
        "[2H][H]", "[H]C#C[H]", "[H]N1C=CC=C1", "[H]n1cccc1", "[H]OC(=O)C([H])([H])N([H])[H]",
        "[H]/N=C(/N)N", "[H]C1([H])CC1", "[H][C@]12CC[C@@H](C1)C2", "C[C@]([H])(O)CC",
        "[H][Si]([H])([H])[H]", "[H]B([H])[H]", "[H]P([H])[H]", "[H]S[H]", "[H][Se][H]",
        "[H]Oc1ccccc1O[H]", "[H]C([H])=C([H])[H]", "[CH4]", "[NH4+]", "[OH-]", "[OH3+]",
        "[CH2]=[CH2]", "[H]c1ccc([H])cc1", "[H]C([H])([H])[C@@]([H])(N)C(=O)O",
        "[H][C]([H])[H]", "[H]N([H])c1ccccc1",
    ]


def hypervalent() -> list[str]:
    return [
        "CS(=O)(=O)C", "CS(C)=O", "O=S(=O)(O)O", "FS(F)(F)(F)(F)F", "F[S](F)(F)(F)(F)F",
        "CS(=O)(=O)N", "C[S+]([O-])C", "O=P(O)(O)O", "CP(C)(C)=O", "C=P(c1ccccc1)(c1ccccc1)c1ccccc1",
        "FP(F)(F)(F)F", "F[P-](F)(F)(F)(F)F", "ClP(Cl)(Cl)(Cl)Cl", "FI(F)F", "FI(F)(F)(F)F",
        "O=I(=O)c1ccccc1", "CC(=O)OI(OC(C)=O)c1ccccc1", "O=[Cl](=O)(=O)O", "O=Cl(=O)(=O)[O-]",
        "[O-][Cl+3]([O-])([O-])[O-]", "OCl(=O)=O", "OCl=O", "O=[Br](=O)O", "FBr(F)F",
        "F[Xe]F", "F[Xe](F)(F)F", "O=[Xe](=O)(=O)=O", "O=[Se](=O)(O)O", "[Se]=C=[Se]",
        "C[Te](C)(Cl)Cl", "O=[N+]([O-])c1ccccc1", "CN(=O)=O", "O=N(=O)c1ccccc1",
        "C[N](=O)=O", "ON(=O)=O", "O=N(O)=O", "CS(=O)(=N)c1ccccc1", "CS(C)(=O)=NC",
        "N#S(F)(F)F", "O=S1(=O)CCCC1", "O=S1(=O)c2ccccc2-c2ccccc21", "c1ccc2c(c1)S(=O)(=O)N2",
        "S(=O)(=O)=O", "O=S=O", "O=[S]=O", "S(F)(F)F[S](F)(F)F", "C[S@](=O)c1ccccc1",
        "C[S@@+]([O-])c1ccccc1", "C[P@](=O)(O)c1ccccc1", "O=P1(O)OCCO1", "c1ccc2c(c1)[I](OC2=O)O",
        "C[As](C)(C)=O", "O=[As](O)(O)O", "[O-][As](=O)([O-])[O-]", "F[B-](F)(F)F", "[BH4-]",
        "C[Si](C)(C)C", "F[Si-2](F)(F)(F)(F)F", "C[Sn](C)(C)C", "O=[Si]=O", "[SiH4]",
        "Fc1ccccc1[S](F)(F)(F)(F)F", "CS(=O)(=O)CS(=O)(=O)C", "OS(=O)(=O)c1ccc2ccccc2c1",
        "C[S+](C)C", "CS(C)(C)(C)C", "C=S(C)C", "CN=S(=O)(C)C", "O=[S+]c1ccccc1",
    ]


def metals() -> list[str]:
    return [
        "[Fe]", "[Fe+2]", "[Fe+3]", "[Na+].[Cl-]", "[Na]Cl", "[Li]C", "C[Mg]Br", "CC[Zn]CC",
        "C[Pt](Cl)(Cl)N", "N[Pt](N)(Cl)Cl", "[NH3][Pt]([NH3])(Cl)Cl", "N->[Cu]<-N",
        "[NH3]->[Co+3](<-[NH3])(<-[NH3])(<-[NH3])(<-[NH3])<-[NH3]", "O->[Fe]", "C1=CC=CC1[Fe]C1C=CC=C1",
        "c1ccc2c(c1)O[Cu]Oc1ccccc1-2", "[K+].[O-]C(=O)c1ccccc1", "[Ca+2].[O-]C(=O)[O-]",
        "[Mg+2].[Cl-].[Cl-]", "Cl[Ti](Cl)(Cl)Cl", "O=[V](=O)O", "O=[Mn](=O)(=O)[O-].[K+]",
        "[Cr](=O)(=O)(Cl)Cl", "C[Hg]Cl", "C[Hg]C", "[Ag+].[NH3]", "Cl[Au]", "Cl[Pd]Cl",
        "Cl[Pd](Cl)(<-P(c1ccccc1)(c1ccccc1)c1ccccc1)<-P(c1ccccc1)(c1ccccc1)c1ccccc1",
        "[Rh]", "C[Al](C)C", "CCCC[Li]", "C[Si](C)(C)[Li]", "[Gd+3]", "[U]", "[Cu+]",
        "[Zn+2].[O-]S(=O)(=O)[O-]", "O=C=[Fe](=C=O)(=C=O)(=C=O)=C=O", "[C-]#[O+]->[Ni](<-[C-]#[O+])(<-[C-]#[O+])<-[C-]#[O+]",
        "CC(=O)O[Pb](OC(C)=O)(OC(C)=O)OC(C)=O", "[Sn](Cl)(Cl)(Cl)Cl", "[Ga+3]", "[Bi]",
        "c1ccc(cc1)[Sn](c1ccccc1)(c1ccccc1)Cl", "C1CC[Ru]CC1", "[Ir]", "[Os](=O)(=O)(=O)=O",
        "N1C=CC=C1->[Zn]", "c1cc2ccc3cccc4ccc(c1)c2c34.[Na]", "O=C(O[Na])c1ccccc1", "C(=O)([O-])[O-].[Na+].[Na+]",
        "[Na+].[Na+].[O-]S([O-])(=O)=O", "[Cs+].[F-]", "[Ba+2]", "[Sr]", "[Be]", "[Mo]", "[W]",
        "[Co]", "[Ni+2]", "[Pt+2]", "Cl[Hg]Cl", "[Tl+]", "[La]", "[Ce+4]", "[Eu]",
        "C(#N)[Cu]", "c1ccc(-c2cccc[n]2->[Ru+2](<-n2ccccc2-c2ccccn2)<-n2ccccc2)nc1",
        "O->[Mg](<-O)(<-O)<-O", "[Fe]1(Cl)(Cl)Cl.Cl1",
    ]


def nontet() -> list[str]:
    return [
        "F[Pt@SP1](Cl)(Br)I", "F[Pt@SP2](Cl)(Br)I", "F[Pt@SP3](Cl)(Br)I",
        "C[Pt@SP1](Cl)(N)N", "N[Pt@SP2](N)(Cl)Cl", "Cl[Pt@SP1](Cl)([NH3])[NH3]",
        "F[As@TB1](Cl)(Br)(I)C", "F[As@TB5](Cl)(Br)(I)C", "F[As@TB10](Cl)(Br)(I)C",
        "F[P@TB1](F)(F)(F)F", "C[P@TB3](F)(F)(F)F", "F[S@OH1](F)(F)(F)(F)F",
        "F[S@OH2](Cl)(F)(F)(F)F", "C[Co@OH1](N)(N)(N)(N)Cl", "C[Co@OH25](N)(O)(F)(Cl)Br",
        "N[Co@OH8](N)(N)(Cl)(Cl)Cl", "F[Pt@SP1](Cl)Br", "[Pt@SP1](F)(Cl)(Br)I",
        "C[C@SP1](F)(Cl)Br", "F[Pt@](Cl)(Br)I", "C[Fe@OH1](C)(C)(C)(C)C", "C[Fe@TB2](F)(Cl)(Br)I",
        "O[Sn@TB7](C)(C)(C)Cl", "Cl[Pd@SP3](Cl)(N)N", "F[Pt@SP2](F)(Cl)Cl",
        "[NH3][Pt@SP1]([NH3])(Cl)Cl", "C[Ru@OH12](O)(N)(Cl)(Br)I", "F[Cr@OH30](Cl)(Br)(I)(N)O",
    ]


def cumulenes() -> list[str]:
    return [
        "C=C=C", "CC=C=CC", "C[C@H]=C=C[C@H]C", "CC=C=C(C)C", "C=C=C=C", "CC=C=C=CC",
        "C/C=C=C=C/C", "C/C=C=C=C\\C", "C=C=C=C=C", "C=C=C=C=C=C", "O=C=O", "O=C=C=C=O",
        "N=C=N", "C=C=O", "CC(C)=C=O", "S=C=S", "C=C=CBr", "ClC=C=CCl", "F/C=C=C=C/F",
        "C1=C=CCCCCC1", "C1CCC=C=CCC1", "C=C=C1CC1", "CC(=C=C(C)C)C", "OC=C=CO",
        "N#CC=C=CC#N", "C(=C=C=C=C=C=C)", "c1ccccc1C=C=Cc1ccccc1", "C=[N+]=[N-]", "C=N=N",
        "[N-]=[N+]=NC", "C(=C=C=C)C", "CC=C=C=C=CC", "C/C(F)=C=C=C(/F)C",
    ]


def ring_ez() -> list[str]:
    out = []
    for n in range(8, 21):
        out.append(f"C1/C=C/{'C' * (n - 3)}1")
        out.append(f"C1/C=C\\{'C' * (n - 3)}1")
        out.append(f"C/1=C/{'C' * (n - 2)}1")
    out += [
        "C1CCCCCCC/C=C/1", "C1CCCCCC/C=C\\1", "C1=C\\CCCCCC/1", "F/C1=C/CCCCCCC1",
        "C/C1=C(\\C)CCCCCC1", "C1CC/C=C/CC/C=C/1", "C1CCCC/C=C/C=C/CCC1", "C1=CC=CC=CC=C1",
        "C1CCC/C=C/CCC1", "O=C1/C=C/CCCCCCC1", "C1CCCC/C=C/CCCC/C=C/CCC1", "C/C=C1/CCCCC1",
        "C/C=C1\\CCC1", "O/N=C1/CCCCC1", "C/N=C1/CCCC(C)C1", "CC1=CCCCC1", "C1=C/CCCCCC\\1",
        "C1C/C=C\\C/C=C\\C1", "c1ccc2c(c1)/C=C\\c1ccccc1/C=C\\2",
    ]
    return out


def macrocycles() -> list[str]:
    out = []
    for n in range(12, 61, 4):
        out.append("C1" + "C" * (n - 1) + "1")
        out.append("O1" + "CCO" * (n // 3 - 1) + "CC1")
    out += [
        "C1CCOCCOCCOCCOCCOCCO1", "C1CNCCNCCNCCN1", "c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3",
        "C1=CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2",
        "CC[C@H](C)[C@@H]1NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@@H](CC(C)C)NC(=O)[C@H](C)NC1=O",
        "OC[C@H]1O[C@@H]2O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O[C@H]1[C@H](O)[C@H]2O",
        "C1CCCCCCCCCCC2CCCCCCCCCCC12", "c1cc2ccc1CCc1ccc(cc1)CC2", "c1cc2ccc1-c1ccc(cc1)-c1ccc-2cc1",
        "C1CCC2(CC1)CCCCC2", "O=C1CCCCCCCCCCCCCCN1", "C1=CC=CC=CC=CC=CC=CC=CC=C1",
        "C1CC2CCC1CC2", "c1ccc2cc3ccc4cc5ccccc5cc4c3cc2c1",
    ]
    return out


def pahs() -> list[str]:
    out = [
        "c1ccc2ccccc2c1", "c1ccc2cc3ccccc3cc2c1", "c1ccc2cc3cc4ccccc4cc3cc2c1",
        "c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1", "c1ccc2cc3cc4cc5cc6ccccc6cc5cc4cc3cc2c1",
        "c1ccc2c(c1)ccc1ccccc12", "c1ccc2c(c1)ccc1c2ccc2ccccc21", "c1cc2ccc3cccc4ccc(c1)c2c34",
        "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61", "c1ccc2c(c1)c1cccc3cccc2c31",
        "c1ccc2c(c1)-c1cccc3cccc-2c13", "c1cc2cccc3c4cccc5cccc(c(c1)c23)c54",
        "c1ccc2c(c1)c1ccccc1c1ccccc21", "c1cc2ccc3ccc4ccc5cccc6c(c1)c2c3c4c56",
        "c1ccc2cc3c(ccc4ccccc43)cc2c1", "c1ccc2c(c1)ccc1c3ccccc3ccc21",
        "c1cc2c3c(c1)ccc1cccc(c13)C2", "c1ccc2c(c1)-c1cccc3c1c-2cc1ccccc13",
        "C1=Cc2cccc3cccc1c23", "c1ccc2c(c1)C=Cc1ccccc1-2",
        "c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c1c6c6c%11c7c2c2c7c8c3c3c8c9c4c4c9c%10c5c5c1c1c6c6c%11c2c2c7c3c3c8c4c4c9c5c1c1c6c2c3c41",
        "c1cc2ccc3ccc4ccc5ccc6ccc7ccc1c1c2c3c4c5c6c71",
    ]
    return out


def bridgehead_quaternary() -> list[str]:
    return [
        "C1C[C@H]2CC[C@@H]1C2", "C1C[C@@H]2CC[C@H]1C2", "C[C@]12CC[C@H](C1)C2(C)C",
        "C[C@@]12CC[C@@H](CC1=O)C2(C)C", "CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2", "C1C2CC3CC1CC(C2)C3",
        "C1[C@@H]2C[C@H]3C[C@@H]1C[C@H](C2)C3", "C[C@]1(O)CC[C@@](C)(N)CC1", "C[C@@](F)(Cl)CC",
        "CC[C@](C)(O)c1ccccc1", "C[C@](N)(C(=O)O)c1ccccc1", "C[C@@]1(c2ccccc2)CCCN1",
        "C1CC2CCC1C2", "N12CCC(CC1)CC2", "C1CN2CCN1CC2", "[C@@H]12CC[C@@H](CC1)O2",
        "O[C@@H]1C[C@H]2CC[C@@H]1C2", "C[C@]12CCCC[C@H]1CCCC2", "C[C@@]12CCCC[C@@H]1CCCC2",
        "C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc43)[C@@H]1CC[C@@H]2O",
        "C[C@H]1C[C@@]2(CCCC2)C[C@@H]1C", "C1C[C@]23CC[C@@]1(CC2)C3", "C[C@]1(F)C[C@](C)(Cl)C1",
        "C[C@@H]1C[C@H](C)C1", "C[C@H]1C[C@@H](C)C[C@H](C)C1", "C[C@@]1(O)C[C@@](C)(O)C1",
        "N[C@@]12CC[C@@H](CC1)CC2", "C[N@+]12CCC(CC1)CC2", "C[N+]1(C)CCCC1",
        "C[N@@+](CC)(CCC)CCCC", "C[P@+](CC)(CCC)c1ccccc1", "C[Si@](F)(Cl)c1ccccc1",
        "C1=C[C@@H]2C[C@H]1C=C2", "C1[C@H]2CC[C@@H]1[C@H]1C=CC[C@@H]21",
        "C12C3C4C1C5C2C3C45", "C12C3C1C23", "C1C2CC12", "C12CC1C2", "C1C23CC12C3",
        "[C@@]1234C5C1C2C3C45", "CC1(C)C2CCC1(C)C(=O)C2",
    ]


def labelled() -> list[str]:
    return [
        "[CH3:1][OH:2]", "[C:1]([H:2])([H:3])([H:4])[H:5]", "[NH2:1][CH2:2][C:3](=[O:4])[OH:5]",
        "[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1", "[*:1]C", "[*:1]c1ccccc1[*:2]",
        "C[C@H:7](N)O", "[CH3:3][C@@H:2]([NH2:1])[C:4](=O)O", "[OH:10]C(=O)[CH2:20]Cl",
        "[13CH3:1]O", "[2H:1]C", "[Na+:1].[Cl-:2]", "[N:1]#[N:2]", "[O:1]=[C:2]=[O:3]",
        "[CH2:1]=[CH2:2]", "[C:1]1[C:2][C:3]1", "c1cc[n:1]cc1", "[nH:3]1cccc1", "C[S:9](=O)(=O)C",
        "[CH3:999]C", "[CH3:1]C.[CH3:1]C", "Cl[C@@H:5](F)Br", "[C@@H:1]1(O)CCCCC1",
    ]


def dummies_tiny() -> list[str]:
    return [
        "*", "[*]", "**", "*C*", "*c1ccccc1", "[*]C(=O)[*]", "*N(*)*", "[1*]C", "[2*]CC[3*]",
        "c1cc(*)ccc1*", "*1CCCC1", "*1*****1", "c1cc*cc1", "C(*)(*)(*)*", "[*-]", "[*+]",
        "[*H]", "[*H2]", "*=C", "*#N", "*/C=C/*", "[*]C[C@H](N)C(=O)O", "*C(=O)O*", "*.*",
        "[1*]c1ccc([2*])cc1", "*:c", "c1ccc2c(c1)**2", "C", "N", "O", "F", "Cl", "Br", "I", "S", "P", "B",
        "[C]", "[N]", "[He]", "[Ne]", "[Ar]", "[Kr]", "[Xe]", "[Rn]", "[H+]", "[H-]", "[Li]", "[Li+]",
        "C.C", "O.O", "C.N.O", "[H]O[H].[H]O[H]", "CC", "C=C", "C#C", "C1CC1", "C1C2C12",
        "[C-]", "[C+]", "[C-4]", "[N+4]", "[O+2]", "[Cl+]", "[I-]", "[Se]", "[Te]", "[Po]",
        "[Og]", "[Ts]", "[Lv]", "[Fl]", "[Mc]", "[Nh]", "[Cn]", "[Rg]", "[Ds]", "[Mt]",
        "[Ac]", "[Th]", "[Pa]", "[Np]", "[Pu]", "[Am]", "[Cm]", "[Bk]", "[Cf]", "[Es]",
        "[Fm]", "[Md]", "[No]", "[Lr]", "[Rf]", "[Db]", "[Sg]", "[Bh]", "[Hs]",
        "C.[Na+].[Cl-]", "[NH4+].[NH4+].[O-]S(=O)(=O)[O-]", "Cl.Cl.Cl.N", "C1CC1.C1CC1",
    ]


def odd_valence_charge() -> list[str]:
    return [
        "C[N+](C)(C)C", "[NH4+]", "C[O+](C)C", "[O-][N+]#N", "C[C-](C)C", "C[CH+]C",
        "[CH3-]", "[CH3+]", "C=[N+]=C", "N=[N+]=[N-]", "O=[O+][O-]", "[O-][O+]=O", "C[S-]",
        "C[Se-]", "[B-](C)(C)(C)C", "C[N-]C", "[N-]1C=CC=C1", "C1=CC=C[CH-]1", "[CH-]1C=CC=C1",
        "C1=CC=CC=C[CH+]1", "C1=CC=CC=C[CH-]1", "c1cccccc1", "[cH+]1cccccc1", "[cH-]1ccccccc1",
        "O=C1C=CC(=O)C=C1", "C[N+]([O-])=O", "C[NH+]=O", "[NH2-]", "[NH-]C", "[N-2]", "[O-2]",
        "[S-2]", "[P-3]", "[C+4]", "[Fe-4]", "[Al-](C)(C)(C)C", "[Si+]", "Cl[I-]Cl",
        "[I-].[I-]", "I[I-]I", "[Br-].[Br][Br]", "O=C([O-])[O-]",
        "C[N+](=O)[O-]", "OB(O)O", "O[B-](O)(O)O", "C1=CN=C1",
        "[N+](=O)([O-])[O-]", "N#[N+][O-]", "C#[N+][O-]", "C#[O+]", "[C-]#[C-]", "[CH2-][N+](C)(C)C",
        "C[N+]1=CC=CC=C1", "C=[N+]1C=CC=C1", "c1ccc2[nH+]ccc2c1",
    ]


def heteroaromatics() -> list[str]:
    return [
        "c1cnccn1", "c1ncncn1", "c1cnnnc1", "c1nncnn1", "c1nnn[nH]1", "c1nn[nH]n1",
        "c1cc2[nH]ccc2[nH]1", "c1cc2occc2o1", "c1ccc2[se]ccc2c1", "c1ccc2[te]ccc2c1", "c1cc[se]c1",
        "c1cc[te]c1", "c1ccc2c(c1)[nH]c1ccccc12", "c1ccc2c(c1)oc1ccccc12", "c1ccc2c(c1)sc1ccccc12",
        "c1cc[pH]c1", "c1ccpcc1", "c1cc[as]cc1", "c1cc[sb]cc1", "c1cocn1", "c1cscn1",
        "c1conc1", "c1csnc1", "n1ccsc1", "c1nonc1", "c1nsnc1", "c1ccc2nonc2c1", "O=c1cc[nH]cc1",
        "O=c1ccoc2ccccc12", "O=c1ccc2ccccc2o1", "S=c1cc[nH]cc1", "Nc1ncnc2[nH]cnc12",
        "O=c1[nH]c(=O)c2[nH]cnc2[nH]1", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "c1ccc2c(c1)nc1ccccc1n2",
        "c1ccc2nc3ccccc3cc2c1", "B1OB(O1)", "b1ccccc1", "c1ccc2c(c1)[bH]c1ccccc12", "c1cc2ccc[n+]3ccc(c1)c23",
        "O=c1occc1", "O=c1[nH]cccc1=O", "c1ccc[n+]2ccccc12", "c12ccccc1cc[n+]2", "c1c[nH]c2nccc2c1",
        "c1ccc2c(c1)nc1ccccc12", "[nH]1c2ccccc2c2ncccc12",
        "C1=Cc2ccccc2N1", "c1cc[nH]c1-c1ccc[nH]1", "c1ccc(-c2ncccn2)nc1",
    ]


def mutate_rows(chembl: list[str], rng: random.Random, n_target: int) -> list[str]:
    elems = [5, 6, 7, 8, 9, 14, 15, 16, 17, 34, 35, 53, 0, 1, 3, 11, 26, 29, 30, 33, 50, 78]
    out: list[str] = []
    seen: set[str] = set()
    tries = 0
    while len(out) < n_target and tries < n_target * 40:
        tries += 1
        base = Chem.MolFromSmiles(rng.choice(chembl))
        if base is None or base.GetNumAtoms() == 0:
            continue
        rw = Chem.RWMol(base)
        if rng.random() < 0.6:
            try:
                Chem.Kekulize(rw, clearAromaticFlags=True)
            except Exception:  # noqa: BLE001
                continue
        n_edits = rng.choice([1, 1, 1, 2, 2, 3])
        ok = True
        for _ in range(n_edits):
            kind = rng.choice(["elem", "charge", "iso", "rad", "bond", "hs", "chiral",
                               "add", "del", "arom", "map", "frag"])
            na = rw.GetNumAtoms()
            if na == 0:
                ok = False
                break
            a = rw.GetAtomWithIdx(rng.randrange(na))
            if kind == "elem":
                a.SetAtomicNum(rng.choice(elems))
            elif kind == "charge":
                a.SetFormalCharge(rng.choice([-2, -1, -1, 1, 1, 2, 3]))
            elif kind == "iso":
                a.SetIsotope(rng.choice([1, 2, 3, 13, 14, 15, 18, 99, 200]))
            elif kind == "rad":
                a.SetNumRadicalElectrons(rng.choice([1, 1, 2]))
                a.SetNoImplicit(True)
            elif kind == "hs":
                a.SetNumExplicitHs(rng.choice([0, 1, 2, 3]))
                a.SetNoImplicit(rng.random() < 0.5)
            elif kind == "chiral":
                a.SetChiralTag(rng.choice([Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                                           Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
                                           Chem.ChiralType.CHI_UNSPECIFIED]))
            elif kind == "bond" and rw.GetNumBonds():
                b = rw.GetBondWithIdx(rng.randrange(rw.GetNumBonds()))
                b.SetBondType(rng.choice([Chem.BondType.SINGLE, Chem.BondType.DOUBLE,
                                          Chem.BondType.TRIPLE, Chem.BondType.DATIVE]))
                b.SetIsAromatic(False)
            elif kind == "add":
                new = rw.AddAtom(Chem.Atom(rng.choice(elems)))
                rw.AddBond(a.GetIdx(), new, rng.choice([Chem.BondType.SINGLE, Chem.BondType.DOUBLE]))
            elif kind == "del" and na > 1:
                rw.RemoveAtom(a.GetIdx())
            elif kind == "arom":
                a.SetIsAromatic(not a.GetIsAromatic())
            elif kind == "map":
                a.SetAtomMapNum(rng.randrange(1, 50))
            elif kind == "frag":
                rw.AddAtom(Chem.Atom(rng.choice([11, 17, 19, 35, 8, 26])))
        if not ok:
            continue
        try:
            rw.UpdatePropertyCache(strict=False)
            if rng.random() < 0.3:
                smi = Chem.MolToSmiles(rw, canonical=False, doRandom=True)
            else:
                smi = Chem.MolToSmiles(rw, canonical=rng.random() < 0.5)
        except Exception:  # noqa: BLE001
            continue
        if not smi or smi in seen or len(smi) > 400:
            continue
        if Chem.MolFromSmiles(smi) is None:
            continue
        seen.add(smi)
        out.append(smi)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--chembl", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--mutations", type=int, default=3000)
    args = ap.parse_args()
    assert rdBase.rdkitVersion == "2026.03.1", rdBase.rdkitVersion
    rng = random.Random(SEED)
    families = [annulenes, charged_aromatics, radicals, isotopes, explicit_h, hypervalent, metals,
                nontet, cumulenes, ring_ez, macrocycles, pahs, bridgehead_quaternary, labelled,
                dummies_tiny, odd_valence_charge, heteroaromatics]
    rows: list[str] = []
    seen: set[str] = set()

    def add(smi: str) -> None:
        if smi and smi not in seen and Chem.MolFromSmiles(smi) is not None:
            seen.add(smi)
            rows.append(smi)

    for fam in families:
        for smi in fam():
            add(smi)
    # Random atom orders of the hand-written rows.
    for smi in list(rows):
        m = Chem.MolFromSmiles(smi)
        if m.GetNumAtoms() > 1 and rng.random() < 0.35:
            for r in Chem.MolToRandomSmilesVect(m, 1, randomSeed=rng.randrange(1 << 30)):
                add(r)
    chembl = [l.split()[0] for l in args.chembl.read_text().splitlines() if l.strip()]
    for smi in mutate_rows(chembl, rng, args.mutations):
        add(smi)
    args.output.write_text("".join(s + "\n" for s in rows))
    print(f"{len(rows)} rows -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
