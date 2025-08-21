# tests/test_mean_and_std.py
import importlib
import numpy as np

def test_compute_alpha_stats_mean_matches_central_value():
    mod = importlib.import_module("mean_and_std")
    G0 = mod.G0
    c = mod.c
    rad_to_uas = 206264.806 * 1e6

    # pick Earth as a stable case
    mass = mod.all_bodies["earth"]["mass"]
    radius_km = mod.all_bodies["earth"]["radius"]
    b_m = radius_km * 1000.0

    mean, std, _, _ = mod.compute_alpha_stats(mass, radius_km, samples=5000)
    expected_uas = (4 * G0 * mass / (c**2 * b_m)) * rad_to_uas

    # 3% tolerance due to finite sampling and rounding
    assert abs(mean - expected_uas) / expected_uas < 0.03
    assert std > 0.0

def test_sigma_scaling_for_std_like_behavior():
    mod = importlib.import_module("mean_and_std")
    mass = mod.all_bodies["earth"]["mass"]
    radius_km = mod.all_bodies["earth"]["radius"]

    # Compare std at two different sigmas
    mean1, std1, *_ = mod.compute_alpha_stats(mass, radius_km, sigma=mod.G0*1e-5, samples=4000)
    mean2, std2, *_ = mod.compute_alpha_stats(mass, radius_km, sigma=mod.G0*5e-5, samples=4000)
    assert std2 > std1
    assert std2 / std1 > 4.0  # roughly proportional, allow noise
