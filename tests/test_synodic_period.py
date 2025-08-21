# tests/test_synodic_period.py
import importlib
import math

def test_inferior_and_superior_examples():
    mod = importlib.import_module("synodic_period")

    T_earth = mod.planet_periods["Earth"]
    T_mercury = mod.planet_periods["Mercury"]   # inferior
    T_venus = mod.planet_periods["Venus"]       # inferior
    T_mars = mod.planet_periods["Mars"]         # superior

    S_mercury = mod.synodic_period_inferior(T_mercury, T_earth)
    S_venus = mod.synodic_period_inferior(T_venus, T_earth)
    S_mars = mod.synodic_period_superior(T_mars, T_earth)

    assert math.isclose(S_mercury, 115.88, rel_tol=3e-3)   # ~0.3%
    assert math.isclose(S_venus,   583.92, rel_tol=3e-3)
    assert math.isclose(S_mars,    779.94, rel_tol=4e-3)
