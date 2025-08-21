def calculate_light_deflection_simple(M, b):
    """
    Calculate the gravitational light deflection angle (in radians)
    using the general relativity formula: alpha = 4GM / (c^2 * b)

    Parameters:
    - M: mass of the lensing body (kg)
    - b: impact parameter (m), typically the radius of the body
    Returns:
    - alpha in radians
    """
    G = 6.67430e-11  # gravitational constant, m^3 / kg / s^2
    c = 299792458  # speed of light, m / s
    alpha = 4 * G * M / (c ** 2 * b)
    return alpha


# Data for 59 celestial objects: mass (kg) and radius (km)
# Source: Geng's dissertation, Table 5.2 and 5.3
all_bodies = {
    # Sun & Planets
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

    # Jupiter’s moons
    'io': {'mass': 8.9319e22, 'radius': 1821.6},
    'europa': {'mass': 4.7998e22, 'radius': 1560.8},
    'ganymede': {'mass': 1.4819e23, 'radius': 2631.2},
    'callisto': {'mass': 1.0759e23, 'radius': 2410.3},

    # Saturn’s moons
    'titan': {'mass': 1.3455e23, 'radius': 2574.7},
    'rhea': {'mass': 2.31e21, 'radius': 764.3},
    'iapetus': {'mass': 1.81e21, 'radius': 735.6},
    'dione': {'mass': 1.1e21, 'radius': 561.7},
    'tethys': {'mass': 6.18e20, 'radius': 533.0},
    'enceladus': {'mass': 1.08e20, 'radius': 252.1},
    'mimas': {'mass': 3.79e19, 'radius': 198.2},
    'phoebe': {'mass': 8.3e18, 'radius': 106.5},
    'hyperion': {'mass': 5.6e18, 'radius': 135.0},
    'janus': {'mass': 1.9e18, 'radius': 89.5},
    'epimetheus': {'mass': 5.3e17, 'radius': 58.1},
    'prometheus': {'mass': 1.6e17, 'radius': 43.1},
    'pandora': {'mass': 1.4e17, 'radius': 40.7},

    # Uranus’s moons
    'titania': {'mass': 3.42e21, 'radius': 788.9},
    'oberon': {'mass': 2.8834e21, 'radius': 761.4},
    'ariel': {'mass': 1.2948e21, 'radius': 578.9},
    'umbriel': {'mass': 1.2214e21, 'radius': 584.7},
    'miranda': {'mass': 6.5941e19, 'radius': 235.8},

    # Neptune’s moons
    'triton': {'mass': 2.1395e22, 'radius': 1353.4},
    'proteus': {'mass': 5.0355e19, 'radius': 210.0},
    'nereid': {'mass': 3.0873e19, 'radius': 170.0},
    'larissa': {'mass': 4.9456e18, 'radius': 97.0},
    'galatea': {'mass': 3.7467e18, 'radius': 88.0},
    'despina': {'mass': 2.0981e18, 'radius': 75.0},
    'thalassa': {'mass': 3.7467e17, 'radius': 41.0},
    'naiad': {'mass': 1.9483e17, 'radius': 33.0},
    'halimede': {'mass': 8.992e16, 'radius': 31.0},
    'neso': {'mass': 1.6485e17, 'radius': 30.0},
    'sao': {'mass': 8.992e16, 'radius': 22.0},

    # Pluto’s moons
    'charon': {'mass': 1.5466e21, 'radius': 603.6},
    'hydra': {'mass': 9.8912e17, 'radius': 36.0},

    # Asteroids
    'ceres': {'mass': 9.47e20, 'radius': 476.2},
    'pallas': {'mass': 2.14e20, 'radius': 272.5},
    'vesta': {'mass': 2.59e20, 'radius': 262.7},
    'hygiea': {'mass': 1.05e20, 'radius': 203.6},
    'interamnia': {'mass': 7.49e19, 'radius': 153.2},
    'psyche': {'mass': 2.29e19, 'radius': 113.0},
    'kalliope': {'mass': 7.36e18, 'radius': 83.8},
    'camilla': {'mass': 1.12e18, 'radius': 105.2}
}

# Output
print("Maximum gravitational light deflection for 59 solar system bodies:\n")

for name, info in all_bodies.items():
    M = info['mass']
    b_km = info['radius']
    b_m = b_km * 1000

    alpha_rad = calculate_light_deflection_simple(M, b_m)
    alpha_microarcsec = alpha_rad * 206264.806 * 1e6

    print(f"{name:<12}: {alpha_microarcsec:10.2f} μas")
