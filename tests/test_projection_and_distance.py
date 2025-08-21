# tests/test_projection_and_distance.py
import numpy as np
import importlib

def test_distance_between():
    mod = importlib.import_module("projection_and_distance")
    a = np.array([0.0, 0.0, 0.0])
    b = np.array([3.0, 4.0, 0.0])
    assert mod.distance_between(a, b) == 5.0

def test_project_point_on_line_perpendicularity():
    mod = importlib.import_module("projection_and_distance")
    E = np.array([0.0, 0.0, 0.0])
    C = np.array([1.0, 0.0, 0.0])
    P = np.array([0.5, 3.0, 0.0])

    Pp = mod.project_point_on_line(E, C, P)
    # vector from P' to P should be perpendicular to EC
    EC = C - E
    PPp = P - Pp
    dot = np.dot(EC, PPp)
    assert abs(dot) < 1e-12

def test_impact_parameter_is_shortest_distance():
    mod = importlib.import_module("projection_and_distance")
    E = np.array([0.0, 0.0, 0.0])
    C = np.array([0.0, 10.0, 0.0])
    P = np.array([3.0, 5.0, 0.0])
    Pp = mod.project_point_on_line(E, C, P)
    b = mod.distance_between(P, Pp)
    # brute force: any other point on line farther than b
    t_vals = np.linspace(-5, 15, 101)
    distances = [np.linalg.norm(P - (E + t*(C-E))) for t in t_vals]
    assert b == min(round(d, 12) for d in distances)
