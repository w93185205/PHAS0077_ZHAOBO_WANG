# -*- coding: utf-8 -*-
"""
Perturbation Duration — Dual Mode
1) mode="fixed_source": Fixed radio source direction, directly compute τ10, τ1, τ0.1 from α(t) in time domain
2) mode="beta_sweep":   Sweep β (no fixed source), estimate the “upper-limit” type duration
"""

import math
import numpy as np
import pandas as pd
import datetime as dt

from skyfield.api import load, Star
from skyfield.units import Angle

# -----------------------------
# Configuration
# -----------------------------
mode = "beta_sweep"          # "fixed_source" or "beta_sweep"

START_DATE = '2025-01-01'    # Start date (UTC)
STEP_HOURS = 6               # Sampling step in hours for fixed_source mode
WINDOW_SCALE = 1.0           # Time window = synodic period * this factor
OBSERVER = 'earth'           # 'earth' or 'ssb' (approximate barycenter)

# Target source direction in fixed_source mode
SOURCE_RA_DEG  = 180.0
SOURCE_DEC_DEG = 30.0

# Angle grid for beta_sweep mode
BETA_DEG_MIN, BETA_DEG_MAX, BETA_STEPS = 0.01, 179.99, 20000

# Thresholds (μas)
THRESHOLDS = [10.0, 1.0, 0.1]

# “…” criterion: coverage ≥ 90%
ELLIPSIS_COVERAGE = 0.90

# Synodic periods (days)
PLANET_SYNODIC = {
    'mercury': 115.88,
    'venus':   583.92,
    'mars':    779.94,
    'jupiter': 398.88,
    'saturn':  378.09,
    'uranus':  369.66,
    'neptune': 367.49,
    'pluto':   366.73,
    'moon':     27.32,
}

# Physical constants
G = 6.67430e-11
c = 299792458.0
rad2uas = 206264.806e6

# Skyfield initialization
ts = load.timescale()
eph = load('de440s.bsp')

if OBSERVER.lower() == 'earth':
    OBS = eph['earth']
elif OBSERVER.lower() == 'ssb':
    # Construct a dummy observer at Solar System barycenter
    class _SSB:
        def at(self, t):
            return eph['sun'].at(t) - eph['sun'].at(t)
    OBS = _SSB()
else:
    raise ValueError("OBSERVER must be 'earth' or 'ssb'")

star = Star(ra=Angle(degrees=SOURCE_RA_DEG), dec=Angle(degrees=SOURCE_DEC_DEG))

# Target objects (in de440s some require barycenter)
TARGETS = [
    ('Mercury', 'mercury'),
    ('Venus',   'venus'),
    ('Mars',    'mars barycenter'),     # de440s requires barycenter
    ('Jupiter', 'jupiter barycenter'),
    ('Saturn',  'saturn barycenter'),
    ('Uranus',  'uranus barycenter'),
    ('Neptune', 'neptune barycenter'),
    ('Pluto',   'pluto barycenter'),
    ('Moon',    'moon'),
]

def get_mass_kg(sf_key):
    # Mass of each body (kg)
    masses = {
        'mercury':            3.3011e23,
        'venus':              4.8675e24,
        'earth':              5.97219e24,
        'mars barycenter':    6.4171e23,
        'jupiter barycenter': 1.8982e27,
        'saturn barycenter':  5.6834e26,
        'uranus barycenter':  8.6810e25,
        'neptune barycenter': 1.02413e26,
        'pluto barycenter':   1.303e22,
        'moon':               7.342e22,
    }
    return masses[sf_key]

def parse_ymd(date_str):
    # Parse date string "YYYY-MM-DD"
    y, m, d = map(int, date_str.split('-'))
    return y, m, d

def make_time_grid(start_date_str, days, step_hours):
    # Generate a time grid
    y, m, d = parse_ymd(start_date_str)
    start_dt = dt.datetime(y, m, d, 0, 0, 0)
    step = dt.timedelta(hours=step_hours)
    n = int(math.ceil(days * 24.0 / step_hours)) + 1
    dts = [start_dt + i * step for i in range(n)]
    years   = [t.year   for t in dts]
    months  = [t.month  for t in dts]
    days_   = [t.day    for t in dts]
    hours   = [t.hour   for t in dts]
    minutes = [t.minute for t in dts]
    seconds = [t.second for t in dts]
    return ts.utc(years, months, days_, hours, minutes, seconds)

def observe_vec(obs, target_key, times):
    # Observation vector (from observer to target)
    if target_key == 'star':
        return obs.at(times).observe(star).apparent()
    else:
        return obs.at(times).observe(eph[target_key]).apparent()

def compute_beta_r(planet_vec, star_vec):
    # Compute β (angle) and r (distance from observer origin)
    p_vec = planet_vec.position.km * 1000.0
    s_vec = star_vec.position.km * 1000.0
    r = np.linalg.norm(p_vec, axis=0)
    p_hat = p_vec / r
    s_hat = s_vec / np.linalg.norm(s_vec, axis=0)
    dot = np.clip(np.sum(p_hat * s_hat, axis=0), -1.0, 1.0)
    beta = np.arccos(dot)
    return beta, r

def alpha_uas(beta, r, M):
    # Compute deflection angle (μas) using formula
    eps = 1e-15
    sinb = np.sin(beta)
    sinb = np.where(np.abs(sinb) < eps, np.sign(sinb) * eps + (sinb==0)*eps, sinb)
    a_rad = (2.0 * G * M / (c**2 * r * sinb)) * (1.0 + np.cos(beta))
    return a_rad * rad2uas

def layered_tau_days(alpha_arr, dt_days, thresholds, cover=ELLIPSIS_COVERAGE):
    # Layered duration statistics τ10, τ1, τ0.1
    total_days = len(alpha_arr) * dt_days
    m10 = alpha_arr > thresholds[0]
    m1  = alpha_arr > thresholds[1]
    m01 = alpha_arr > thresholds[2]
    d10 = m10.sum() * dt_days
    d1  = m1.sum()  * dt_days
    d01 = m01.sum() * dt_days
    t10 = d10
    t1  = max(d1  - d10, 0.0)
    t01 = max(d01 - d1 , 0.0)
    out10 = '...' if (d10 / total_days) >= cover else round(t10, 3)
    out1  = '...' if (d1  / total_days) >= cover else round(t1 , 3)
    out01 = '...' if (d01 / total_days) >= cover else round(t01, 3)
    return out10, out1, out01

def rmin_rmax_over_window(sf_key, window_days, step_hours):
    """Get min and max r(t) from time samples (used in beta_sweep mode)"""
    times = make_time_grid(START_DATE, window_days, step_hours)
    vec_planet = observe_vec(OBS, sf_key, times)
    p_vec = vec_planet.position.km * 1000.0
    r = np.linalg.norm(p_vec, axis=0)
    return float(np.min(r)), float(np.max(r))

def tau_via_beta_sweep(M, T_syn, rmin, rmax, thresholds):
    # Estimate durations by sweeping β
    beta_deg = np.linspace(BETA_DEG_MIN, BETA_DEG_MAX, BETA_STEPS)
    beta_rad = np.deg2rad(beta_deg)

    a_min = alpha_uas(beta_rad, rmin, M)
    # a_max = alpha_uas(beta_rad, rmax, M)  # no longer used

    cover_fracs = []
    for thr in thresholds:
        mask = (a_min > thr)  # <-- only use rmin (optimistic upper limit)
        if not np.any(mask):
            cover_fracs.append(0.0)
        else:
            half_span = beta_deg[mask][-1] - beta_deg[mask][0]
            full_span = min(2.0 * half_span, 360.0)
            cover_fracs.append(full_span / 360.0)

    f10, f1, f01 = cover_fracs
    f10_only = max(f10, 0.0)
    f1_only  = max(f1  - f10, 0.0)
    f01_only = max(f01 - f1 , 0.0)

    def fmt(f):
        if f >= ELLIPSIS_COVERAGE:
            return '...'
        days = f * T_syn
        return 0.0 if days == 0.0 else round(float(days), 3)

    return fmt(f10_only), fmt(f1_only), fmt(f01_only)

# -----------------------------
# Main process
# -----------------------------
rows = []
for nice_name, sf_key in TARGETS:
    base_key = sf_key.split()[0] if sf_key != 'moon' else 'moon'
    T_syn = PLANET_SYNODIC[base_key]
    window_days = T_syn * WINDOW_SCALE
    M = get_mass_kg(sf_key)

    if mode == "fixed_source":
        times = make_time_grid(START_DATE, window_days, STEP_HOURS)
        dt_days = STEP_HOURS / 24.0
        vec_planet = observe_vec(OBS, sf_key, times)
        vec_star   = observe_vec(OBS, 'star', times)
        beta, r = compute_beta_r(vec_planet, vec_star)
        alpha = alpha_uas(beta, r, M)
        out10, out1, out01 = layered_tau_days(alpha, dt_days, THRESHOLDS, ELLIPSIS_COVERAGE)

    elif mode == "beta_sweep":
        # First compute r_min / r_max in time domain (observer dependent), then sweep β
        rmin, rmax = rmin_rmax_over_window(sf_key, window_days, STEP_HOURS)
        out10, out1, out01 = tau_via_beta_sweep(M, T_syn, rmin, rmax, THRESHOLDS)

    else:
        raise ValueError("mode must be 'fixed_source' or 'beta_sweep'")

    rows.append([nice_name, out10, out1, out01])

df = pd.DataFrame(rows, columns=['Object', 'τ₁₀ (days)', 'τ₁.₀ (days)', 'τ₀.₁ (days)'])
print(df.to_string(index=False))
