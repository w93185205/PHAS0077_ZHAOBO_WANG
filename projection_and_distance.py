import numpy as np

def distance_between(a, b):
    """
    Calculate the Euclidean distance between two points a and b in 3D space.
    """
    return np.linalg.norm(a - b)

def project_point_on_line(E, C, P):
    """
    Project point P onto the line defined by points E and C in 3D space.
    Returns: The projected point P' as a numpy array
    """
    EC = C - E
    EP = P - E
    t = np.dot(EC, EP) / np.dot(EC, EC)
    P_prime = E + t * EC
    return P_prime

all_bodies = {
    'sun': {'ssb': np.array([-1300684.80446007,  -215295.30096281,   -58332.18929034])},
    'mercury': {'ssb': np.array([-9886075.93833861, 39856428.3702694,  22237843.66320414])},
    'venus': {'ssb': np.array([-34596104.40665917, -95233416.0950379,  -40706214.01096347])},
    'earth': {'ssb': np.array([ 2.27521239e+07, -1.37998606e+08, -5.97859794e+07])},
    'moon': {'ssb': np.array([ 2.26378647e+07, -1.38312140e+08, -5.99449984e+07])},
    'mars barycenter': {'ssb': np.array([-2.45888805e+08,  3.67176291e+07,  2.34812462e+07])},
    'jupiter barycenter': {'ssb': np.array([6.48508646e+08, 3.34096342e+08, 1.27419348e+08])},
    'saturn barycenter': {'ssb': np.array([ 1.28746527e+09, -6.21450035e+08, -3.12143424e+08])},
    'uranus barycenter': {'ssb': np.array([1.91832564e+09, 2.04749955e+09, 8.69616983e+08])},
    'neptune barycenter': {'ssb': np.array([ 4.45781956e+09, -2.87584928e+08, -2.28694152e+08])},
    'pluto barycenter': {'ssb': np.array([ 2.49588934e+09, -4.09264531e+09, -2.02919358e+09])},
}

multiplier = 1e6
earth_name = 'earth'
E = all_bodies[earth_name]['ssb']

for name, info in all_bodies.items():
    if name == earth_name:
        continue
    P = info['ssb']
    C = P * multiplier  # CES position as per the paper's definition

    P_prime = project_point_on_line(E, C, P)
    b = distance_between(P, P_prime)

    print(f"{name}:")
    print(f"  Projection of planet on line EC (P'): {P_prime}")
    print(f"  Impact parameter b (km): {b}\n")
