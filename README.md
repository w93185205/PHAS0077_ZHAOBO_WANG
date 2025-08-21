# PHAS0077_ZHAOBO_WANG


# Gravitational Light Deflection in the Solar System

This repository provides a Python-based framework to study gravitational light deflection caused by Solar System objects.  
It implements analytic formulas, numerical calculations, and statistical simulations for light bending effects and impact parameters, and includes automated tests for validation.

## Project Structure

├── alpha_and_beta_max_and_min.py # Compute simplified vs. normal deflection angles α and β  
├── ces_selection.py # Select objects/events for case studies (CES = Candidate Event Selection)  
├── data.py # Utility data and constants (e.g., astronomical units)  
├── dual_deflection.py # Model dual-body deflection scenarios  
├── impact_ranges.py # Evaluate ranges of impact parameters  
├── maximum_deflection_alphamax.py # Compute maximum deflection angle (α_max)  
├── mass_and_radius.py # Planetary masses and radii dictionary  
├── mean_and_std.py # Monte Carlo statistics (mean/std) for deflection values  
├── projection_and_distance.py # Geometry utilities: distance, projection, impact parameter  
│  
├── tests/ # Unit tests (pytest-based)  
│ ├── test_alpha_and_beta_max_and_min.py  
│ ├── test_mass_and_radius.py  
│ ├── test_maximum_deflection_alphamax.py  
│ ├── test_mean_and_std.py  
│ └── test_projection_and_distance.py  


## Requirements

- Python 3.9+  
- NumPy  
- PyTest (for testing)  

Install dependencies:

```bash
pip install -r requirements.txt




## Usage
1: Light Deflection at the Solar Limb：

import mass_and_radius as mnr
from maximum_deflection_alphamax import calculate_light_deflection_simple

M = mnr.planet_mass["sun"]
R = mnr.planet_radius["sun"] * 1000.0  # convert km → m

alpha_rad = calculate_light_deflection_simple(M, R)
print("Deflection at solar limb [rad]:", alpha_rad)



2：Monte Carlo Statistics for Earth

from mean_and_std import compute_alpha_stats
import mass_and_radius as mnr

mean, std, samples, alphas = compute_alpha_stats(
    mass=mnr.planet_mass["earth"],
    radius_km=mnr.planet_radius["earth"],
    sigma=1e-5
)

print("Mean deflection [μas]:", mean)
print("Standard deviation [μas]:", std)



3：Impact Parameter Geometry

import numpy as np
from projection_and_distance import distance_between, project_point_on_line

E = np.array([0.0, 0.0, 0.0])  # Earth
C = np.array([0.0, 10.0, 0.0]) # Celestial object
P = np.array([3.0, 5.0, 0.0])  # Photon path point

Pp = project_point_on_line(E, C, P)
b = distance_between(P, Pp)

print("Impact parameter (shortest distance) =", b)





## RunningTest

To validate correctness, run the unit tests with pytest:
pytest tests/


This runs all consistency checks, including:

Simplified vs. normal α comparison,Maximum deflection scaling with mass and distance, planet mass/radius dictionary validation

Monte Carlo statistics consistency

Projection and geometry correctness
