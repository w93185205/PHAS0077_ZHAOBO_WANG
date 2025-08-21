import sys
import os
import pytest

# 让 pytest 能找到上一级目录的代码文件
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import mass_and_radius as mnr


def test_planet_dicts_have_expected_keys():
    """Check that planet dictionaries have expected keys"""
    assert "sun" in mnr.planet_mass
    assert "jupiter barycenter" in mnr.planet_mass
    assert "earth" in mnr.planet_radius


def test_masses_and_radii_positive():
    """Check all masses and radii are positive"""
    for m in mnr.planet_mass.values():
        assert m > 0
    for r in mnr.planet_radius.values():
        assert r > 0


def test_name_consistency_subset():
    """Check that planets appearing in both dicts are consistent"""
    for name in set(mnr.planet_mass.keys()) & set(mnr.planet_radius.keys()):
        assert isinstance(mnr.planet_mass[name], (int, float))
        assert isinstance(mnr.planet_radius[name], (int, float))
