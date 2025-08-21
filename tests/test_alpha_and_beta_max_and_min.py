# tests/test_alpha_and_beta_max_and_min.py
import importlib
import numpy as np

def test_alpha_simplified_vs_normal_at_90deg():
    mod = importlib.import_module("alpha_and_beta_max_and_min")
    # Monkeypatch module-level beta_rad to fixed 90 deg
    beta_rad_backup = mod.beta_rad.copy()
    try:
        mod.beta_rad = np.deg2rad(np.array([90.0]))
        M, r = 5.97237e24, 1.0 * mod.AU
        a_s, a_n, a_d = mod.compute_deflection_angles(M, r)
        # Each is an array length 1
        assert np.isclose(a_s[0], 2.0 * a_n[0], rtol=1e-6)
        assert a_d[0] >= -1e-12
    finally:
        mod.beta_rad = beta_rad_backup

def test_alpha_simplified_always_ge_normal_for_all_beta():
    mod = importlib.import_module("alpha_and_beta_max_and_min")
    M, r = 5.97237e24, 1.0 * mod.AU
    a_s, a_n, a_d = mod.compute_deflection_angles(M, r)
    assert np.all(a_s >= a_n - 1e-12)
    assert np.all(a_d >= -1e-12)
