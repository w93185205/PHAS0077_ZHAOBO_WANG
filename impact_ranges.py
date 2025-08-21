import numpy as np
import pandas as pd
from scipy.optimize import root_scalar

# === Physical constants ===
G = 6.67430e-11  # m^3 / kg / s^2
c = 299792458  # m/s
AU = 1.496e11  # m
deg = 180 / np.pi  # rad -> deg
arcsec = deg * 3600  # rad -> arcsec

# === Extended celestial body masses (kg) ===
planet_mass = {
    'sun': 1.9885e30,
    'mercury': 3.3011e23,
    'venus': 4.8675e24,
    'moon': 7.342e22,
    'mars': 6.4171e23,
    'jupiter': 1.8982e27,
    'saturn': 5.6834e26,
    'uranus': 8.6810e25,
    'neptune': 1.02413e26,
    'pluto': 1.303e22,
    # Newly added small bodies and satellites
    'ceres': 9.3835e20,
    'pallas': 2.11e20,
    'vesta': 2.59076e20,
    'ganymede': 1.4819e23,
    'callisto': 1.0759e23,
    'io': 8.9319e22,
    'europa': 4.7998e22,
    'titan': 1.3452e23,
    'eris': 1.66e22
}

# === Average distances of celestial bodies from Earth (AU or meters) ===
planet_distance_au = {
    'sun': 1e-5,
    'mercury': 0.39,
    'venus': 0.72,
    'moon': 1.00,  # special case: 384400 km
    'mars': 1.52,
    'jupiter': 5.20,
    'saturn': 9.58,
    'uranus': 19.2,
    'neptune': 30.05,
    'pluto': 39.48,
    # Estimated average geocentric distances for newly added bodies (AU)
    'ceres': 2.77,
    'pallas': 2.77,
    'vesta': 2.36,
    'ganymede': 5.20,
    'callisto': 5.20,
    'io': 5.20,
    'europa': 5.20,
    'titan': 9.58,
    'eris': 67.67
}

# === Deflection angle thresholds (microarcseconds) and in radians ===
alpha_thresholds_uas = [10, 1, 0.1]
alpha_thresholds_rad = [a * 1e-6 / 206264.806 for a in alpha_thresholds_uas]

# === Deflection angle function ===
def deflection_angle(beta_rad, M, r):
    sin_beta = np.sin(beta_rad)
    cos_beta = np.cos(beta_rad)
    return (2 * G * M / (c ** 2 * r * sin_beta)) * (cos_beta + 1)

# === Function to solve for β ===
def solve_beta(M, r, alpha_thresh_rad):
    def f(beta):
        return deflection_angle(beta, M, r) - alpha_thresh_rad
    try:
        sol = root_scalar(f, bracket=[1e-8, np.pi - 1e-8], method='brentq')
        return sol.root if sol.converged else np.nan
    except:
        return np.nan

# === Calculate β ranges ===
results = []
for name in planet_mass:
    clean_name = name.capitalize()
    M = planet_mass[name]
    R = planet_distance_au[name]

    if name == 'sun':
        r_min = r_max = AU
    elif name == 'moon':
        r_min = r_max = 384400e3
    else:
        r_min = abs(R - 1.0) * AU
        r_max = (R + 1.0) * AU

    row = {'body': clean_name}
    for alpha_uas, alpha_rad in zip(alpha_thresholds_uas, alpha_thresholds_rad):
        beta_max_rad = solve_beta(M, r_min, alpha_rad)
        beta_min_rad = solve_beta(M, r_max, alpha_rad)

        row[f'beta{alpha_uas}_max_deg'] = round(beta_max_rad * deg, 4) if beta_max_rad else np.nan
        row[f'beta{alpha_uas}_min_deg'] = round(beta_min_rad * deg, 4) if beta_min_rad else np.nan
        row[f'beta{alpha_uas}_max_arcsec'] = round(beta_max_rad * arcsec, 1) if beta_max_rad else np.nan
        row[f'beta{alpha_uas}_min_arcsec'] = round(beta_min_rad * arcsec, 1) if beta_min_rad else np.nan

    results.append(row)

# === Convert to DataFrame and output ===
df = pd.DataFrame(results)
pd.set_option('display.max_columns', None)
print(df.to_string(index=False))

