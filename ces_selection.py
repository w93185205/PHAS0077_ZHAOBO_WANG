from skyfield.api import load
import numpy as np

def get_planet_position(planets, planet_name, t):
    """
    Get the position of the specified body (with respect to the Solar System Barycenter) at time t.
    Returns: numpy array, shape (3,), units km
    """
    return planets[planet_name].at(t).position.km

def compute_CES_position_from_planet(planet_position, multiplier=1e6):
    """
    Given a body's position vector (with respect to the SSB), return the position vector
    of the CES in the same direction, at a distance 'multiplier' times farther.
    """
    return planet_position * multiplier

if __name__ == "__main__":
    # Load ephemeris and time
    planets = load(r"C:\Users\lenovo\Desktop\de440s.bsp")
    ts = load.timescale()
    t = ts.utc(2025, 1, 1, 12, 0, 0)

    # Supported names in de440s.bsp
    body_names = [
        'sun', 'mercury', 'venus', 'earth', 'moon',
        'mars barycenter', 'jupiter barycenter', 'saturn barycenter',
        'uranus barycenter', 'neptune barycenter', 'pluto barycenter'
    ]

    # Map internal names to display names without the word "barycenter"
    name_map = {
        'mars barycenter': 'mars',
        'jupiter barycenter': 'jupiter',
        'saturn barycenter': 'saturn',
        'uranus barycenter': 'uranus',
        'neptune barycenter': 'neptune',
        'pluto barycenter': 'pluto'
    }

    for name in body_names:
        # Use the internal name for lookup, but a clean name for display
        display_name = name_map.get(name, name)
        planet_position = get_planet_position(planets, name, t)
        CES_position = compute_CES_position_from_planet(planet_position, multiplier=1e6)

        print(f"{display_name}:")
        print(f"  SSB position (km): {planet_position}")
        print(f"  CES position (km): {CES_position}\n")
