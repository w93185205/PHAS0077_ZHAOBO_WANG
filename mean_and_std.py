import numpy as np
import scipy.stats as stats
import pandas as pd

# Constants
c = 299792458           # Speed of light (m/s)
G0 = 6.67430e-11        # Gravitational constant (central value)
sigma = 1e-4 * G0       # 10^-4 relative uncertainty

# All body data: mass (kg) and radius (km)
all_bodies = {
    'sun': {'mass': 1.9885e30, 'radius': 695700},
    'mercury': {'mass': 3.3010e23, 'radius': 2439.7},
    'venus': {'mass': 4.8673e24, 'radius': 6051.8},
    'earth': {'mass': 5.9722e24, 'radius': 6378.137},
    'moon': {'mass': 7.346e22, 'radius': 1738.1},
    'mars': {'mass': 6.4169e23, 'radius': 3389.5},
    'jupiter': {'mass': 1.89813e27, 'radius': 69911},
    'saturn': {'mass': 5.6832e26, 'radius': 58232},
    'uranus': {'mass': 8.6811e25, 'radius': 25362},
    'neptune': {'mass': 1.02409e26, 'radius': 24622},
    'pluto': {'mass': 1.303e22, 'radius': 1188},
}

# Function to compute deflection statistics under Gaussian G
def compute_alpha_stats(mass, radius_km, G0=G0, sigma=sigma, samples=10000):
    b_m = radius_km * 1000
    G_samples = np.random.normal(loc=G0, scale=sigma, size=samples)

    alpha_rad = 4 * G_samples * mass / (c ** 2 * b_m)
    alpha_microarcsec = alpha_rad * 206264.806 * 1e6

    mean = np.mean(alpha_microarcsec)
    std = np.std(alpha_microarcsec)
    ci_68 = np.percentile(alpha_microarcsec, [16, 84])
    ci_95 = np.percentile(alpha_microarcsec, [2.5, 97.5])

    return mean, std, ci_68, ci_95

# Run for all bodies
results = []

for name, info in all_bodies.items():
    mean, std, ci_68, ci_95 = compute_alpha_stats(info['mass'], info['radius'])
    results.append({
        'Name': name,
        'α_mean (μas)': round(mean, 2),
        'α_std (μas)': round(std, 2),
        '68% CI (μas)': f"[{ci_68[0]:.2f}, {ci_68[1]:.2f}]",
        '95% CI (μas)': f"[{ci_95[0]:.2f}, {ci_95[1]:.2f}]"
    })

# Output as table
df = pd.DataFrame(results)
df = df.sort_values(by='α_mean (μas)', ascending=False)
print(df.to_string(index=False))
df.to_csv("deflection_stats.csv", index=False)
