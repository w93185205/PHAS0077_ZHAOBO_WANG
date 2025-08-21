# tests/test_maximum_deflection_alphamax.py
import math
import importlib

def test_alpha_formula_matches_definition(phys_constants, sun_params, tolerances):
    mod = importlib.import_module("maximum_deflection_alphamax")
    G, c = phys_constants["G"], phys_constants["c"]

    M = sun_params["mass"]
    b_m = sun_params["radius_km"] * 1000.0

    alpha = mod.calculate_light_deflection_simple(M, b_m)
    expected = 4 * G * M / (c**2 * b_m)
    assert math.isclose(alpha, expected, rel_tol=tolerances["rel"])

def test_alpha_monotonic_mass_and_b(phys_constants):
    mod = importlib.import_module("maximum_deflection_alphamax")

    M1, M2 = 1.0e25, 2.0e25
    b = 1.0e7
    assert mod.calculate_light_deflection_simple(M2, b) > mod.calculate_light_deflection_simple(M1, b)

    M = 1.0e25
    b1, b2 = 1.0e7, 2.0e7
    assert mod.calculate_light_deflection_simple(M, b1) > mod.calculate_light_deflection_simple(M, b2)

def test_solar_limb_is_about_1p75_arcsec(phys_constants, sun_params):
    mod = importlib.import_module("maximum_deflection_alphamax")
    rad_to_uas = phys_constants["rad_to_uas"]

    alpha_rad = mod.calculate_light_deflection_simple(
        sun_params["mass"], sun_params["radius_km"] * 1000.0
    )
    alpha_uas = alpha_rad * rad_to_uas
    # ~1.75 arcsec = 1.75e6 μas, accept 2% tolerance
    assert abs(alpha_uas - 1.75e6) / 1.75e6 < 0.02
