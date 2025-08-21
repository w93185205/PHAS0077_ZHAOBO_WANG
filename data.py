from skyfield.api import load

# 1. Load the JPL DE440s ephemeris file
planets = load(r"C:\Users\lenovo\Desktop\de440s.bsp")
ts = load.timescale()
t = ts.utc(2025, 1, 1, 0, 0, 0)

# 2. List of main bodies supported by de440s.bsp
body_names = [
    'sun', 'mercury', 'venus', 'earth', 'moon',
    'mars barycenter', 'jupiter barycenter', 'saturn barycenter',
    'uranus barycenter', 'neptune barycenter', 'pluto barycenter'
]

# 3. NASA/IAU physical parameters (mass in kg, radius in km, orbital period in days)
planet_mass = {
    'sun': 1.9885e30,
    'mercury': 3.3011e23,
    'venus': 4.8675e24,
    'earth': 5.97237e24,
    'moon': 7.342e22,
    'mars barycenter': 6.4171e23,
    'jupiter barycenter': 1.8982e27,
    'saturn barycenter': 5.6834e26,
    'uranus barycenter': 8.6810e25,
    'neptune barycenter': 1.02413e26,
    'pluto barycenter': 1.303e22
}
planet_radius = {
    'sun': 695700,
    'mercury': 2439.7,
    'venus': 6051.8,
    'earth': 6371.0,
    'moon': 1737.4,
    'mars barycenter': 3389.5,
    'jupiter barycenter': 69911,
    'saturn barycenter': 58232,
    'uranus barycenter': 25362,
    'neptune barycenter': 24622,
    'pluto barycenter': 1188.3
}
planet_period = {
    'mercury': 87.969,        # days
    'venus': 224.701,
    'earth': 365.256,
    'moon': 27.322,
    'mars barycenter': 686.980,
    'jupiter barycenter': 4332.589,
    'saturn barycenter': 10759.22,
    'uranus barycenter': 30685.4,
    'neptune barycenter': 60189,
    'pluto barycenter': 90560
}

# 4. Retrieve and print all parameters for each body
all_bodies = {}
for name in body_names:
    position = planets[name].at(t).position.km  # (x, y, z) in km
    mass = planet_mass.get(name, None)
    radius = planet_radius.get(name, None)
    period = planet_period.get(name, None)
    all_bodies[name] = {
        'position': position,
        'mass': mass,
        'radius': radius,
        'period': period
    }
    print(f"{name}:")
    print(f"  position (km): {position}")
    print(f"  mass (kg): {mass}")
    print(f"  radius (km): {radius}")
    print(f"  period (days): {period}")

