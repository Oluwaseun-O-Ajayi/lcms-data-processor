"""
Glycan target library: m/z values for extracted-ion chromatograms of released N-glycans.

Compositions are written as residue counts, either long form ("HexNAc4Hex5dHex1NeuAc2")
or short form ("N4H5F1S2"). Monoisotopic masses are computed from residue masses
(glycan = sum of residues + H2O), the reducing-end label shift and the adduct.

Label shifts assume reductive amination (label + glycan - H2O + 2H):

    2-AB (2-aminobenzamide, C7H8N2O)        +120.068748
    2-AA (2-aminobenzoic acid, C7H7NO2)     +121.052764
    procainamide (C13H21N3O)                +219.173547
    reduced (alditol, +2H)                  +2.015650
    none (free reducing end)                 0

These are candidate m/z values. Assignments still need retention time and MS/MS evidence.
"""

import re

import pandas as pd

# Monoisotopic element masses (IUPAC/AME)
ELEMENT = {'C': 12.0, 'H': 1.00782503207, 'N': 14.0030740048, 'O': 15.99491461956}
PROTON = 1.007276466812

# Residue (dehydrated monosaccharide) formulas
RESIDUE_FORMULA = {
    'Hex': {'C': 6, 'H': 10, 'O': 5},
    'HexNAc': {'C': 8, 'H': 13, 'N': 1, 'O': 5},
    'dHex': {'C': 6, 'H': 10, 'O': 4},
    'NeuAc': {'C': 11, 'H': 17, 'N': 1, 'O': 8},
    'NeuGc': {'C': 11, 'H': 17, 'N': 1, 'O': 9},
}
SHORT = {'H': 'Hex', 'N': 'HexNAc', 'F': 'dHex', 'S': 'NeuAc', 'G': 'NeuGc'}

LABEL_FORMULA = {
    '2-AB': {'C': 7, 'H': 8, 'N': 2, 'O': 1},
    '2-AA': {'C': 7, 'H': 7, 'N': 1, 'O': 2},
    'procainamide': {'C': 13, 'H': 21, 'N': 3, 'O': 1},
}
ADDUCT = {'H': PROTON, 'Na': 22.98922070, 'NH4': 18.03382555}


def formula_mass(formula):
    """Monoisotopic mass of an element-count dict."""
    return sum(ELEMENT[el] * n for el, n in formula.items())


RESIDUE = {name: formula_mass(f) for name, f in RESIDUE_FORMULA.items()}
WATER = formula_mass({'H': 2, 'O': 1})
H2 = formula_mass({'H': 2})


def label_shift(label):
    """Mass added to the free glycan by the reducing-end label."""
    if label in (None, '', 'none'):
        return 0.0
    if label == 'reduced':
        return H2
    if label not in LABEL_FORMULA:
        raise ValueError(f"Unknown label '{label}'. Options: none, reduced, {', '.join(LABEL_FORMULA)}")
    return formula_mass(LABEL_FORMULA[label]) - WATER + H2


def parse_composition(text):
    """
    Parse a composition string into residue counts.

    >>> parse_composition('N4H5F1S2')
    {'HexNAc': 4, 'Hex': 5, 'dHex': 1, 'NeuAc': 2}
    """
    text = str(text).strip()
    long_form = re.findall(r'(HexNAc|Hex|dHex|NeuAc|NeuGc)(\d+)', text)
    if long_form and ''.join(f'{a}{b}' for a, b in long_form) == text:
        tokens = long_form
    else:
        short = re.findall(r'([HNFSG])(\d+)', text)
        if not short or ''.join(f'{a}{b}' for a, b in short) != text:
            raise ValueError(f"Cannot parse composition '{text}'. Use e.g. N4H5F1S2 or HexNAc4Hex5dHex1NeuAc2")
        tokens = [(SHORT[a], b) for a, b in short]
    counts = {}
    for name, n in tokens:
        if name in counts:
            raise ValueError(f"Residue {name} appears twice in '{text}'")
        counts[name] = int(n)
    return counts


def glycan_mass(composition, label='none'):
    """Neutral monoisotopic mass of a (labelled) glycan."""
    counts = parse_composition(composition) if isinstance(composition, str) else composition
    return sum(RESIDUE[k] * n for k, n in counts.items()) + WATER + label_shift(label)


def glycan_mz(composition, charge=1, label='none', adduct='H'):
    """m/z of [M + z adduct]^z+."""
    if charge < 1:
        raise ValueError("charge must be a positive integer")
    return (glycan_mass(composition, label) + charge * ADDUCT[adduct]) / charge


def target_table(compositions, label='none', charges=(1, 2), adduct='H'):
    """
    Build an XIC target table.

    Args:
        compositions: mapping of glycan name to composition string
        label: reducing-end label (see module docstring)
        charges: charge states to extract
    Returns:
        DataFrame with columns glycan, composition, charge, mz
    """
    rows = [{'glycan': name, 'composition': comp, 'charge': z,
             'mz': glycan_mz(comp, z, label, adduct)}
            for name, comp in compositions.items() for z in charges]
    return pd.DataFrame(rows)


# Common candidate compositions (names use Oxford-style shorthand).
IGG_CANDIDATES = {
    'G0': 'N4H3', 'G0F': 'N4H3F1', 'G1F': 'N4H4F1', 'G2F': 'N4H5F1',
    'G1': 'N4H4', 'G2': 'N4H5', 'G0FB': 'N5H3F1', 'Man5': 'N2H5',
    'G1FS1': 'N4H4F1S1', 'G2FS1': 'N4H5F1S1', 'G2FS2': 'N4H5F1S2',
}
FETUIN_CANDIDATES = {
    'A2G2S1': 'N4H5S1', 'A2G2S2': 'N4H5S2',
    'A3G3S2': 'N5H6S2', 'A3G3S3': 'N5H6S3', 'A3G3S4': 'N5H6S4',
}
