"""Residues for glycan classes beyond N-glycans: masses checked against elemental formulas."""

import pytest

from glycan_targets import glycan_mass, glycan_mz, parse_composition

MASS = {'C': 12.0, 'H': 1.00782503207, 'N': 14.0030740048, 'O': 15.99491461956, 'S': 31.97207100}


def mono(**counts):
    return sum(MASS[el] * n for el, n in counts.items())


def test_free_glucuronic_acid():
    # GlcA, C6H10O7
    assert glycan_mass('HexA1') == pytest.approx(mono(C=6, H=10, O=7), abs=1e-6)


def test_chondroitin_disaccharides():
    # Delta-Di-0S (DeltaHexA-GalNAc), C14H21NO11; Delta-Di-4S adds SO3, C14H21NO14S
    assert glycan_mass('DeltaHexA1HexNAc1') == pytest.approx(mono(C=14, H=21, N=1, O=11), abs=1e-6)
    assert glycan_mass('U1N1Su1') == pytest.approx(mono(C=14, H=21, N=1, O=14, S=1), abs=1e-6)


def test_heparin_trisulfated_disaccharide():
    # DeltaUA2S-GlcNS6S, C12H19NO19S3
    assert glycan_mass('U1Hn1Su3') == pytest.approx(mono(C=12, H=19, N=1, O=19, S=3), abs=1e-6)


def test_negative_ion_and_2amac_label():
    m = mono(C=14, H=21, N=1, O=14, S=1)
    proton = 1.007276466812
    assert glycan_mz('U1N1Su1', charge=1, adduct='-H') == pytest.approx(m - proton, abs=1e-6)
    # 2-AMAC by reductive amination: + C13H10N2O - H2O + 2H
    shift = mono(C=13, H=10, N=2, O=1) - mono(H=2, O=1) + mono(H=2)
    assert glycan_mass('U1N1', label='2-AMAC') == pytest.approx(mono(C=14, H=21, N=1, O=11) + shift, abs=1e-6)


def test_o_glycan_and_pentose():
    # Core 1 sialyl-T (NeuAc-Gal-GalNAc), C25H42N2O19; xylose C5H10O5
    assert glycan_mass('N1H1S1') == pytest.approx(mono(C=25, H=42, N=2, O=19), abs=1e-6)
    assert glycan_mass('P1') == pytest.approx(mono(C=5, H=10, O=5), abs=1e-6)


def test_long_and_short_forms_agree():
    assert parse_composition('DeltaHexA1HexN1Sulfo3') == parse_composition('U1Hn1Su3')
    assert parse_composition('HexNAc4Hex5dHex1NeuAc2') == {'HexNAc': 4, 'Hex': 5, 'dHex': 1, 'NeuAc': 2}
    with pytest.raises(ValueError):
        parse_composition('X1')
