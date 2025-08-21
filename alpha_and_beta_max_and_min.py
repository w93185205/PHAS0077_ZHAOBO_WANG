import numpy as np
import matplotlib.pyplot as plt

# --- Physical constants ---
G = 6.67430e-11  # m^3 / kg / s^2
c = 299792458    # m/s
AU = 1.496e11    # m
rad2uas = 206264.806 * 1e6  # Radian to microarcseconds

# --- Masses of the eight planets (kg) and their average heliocentric distances (AU) ---
planet_data = {
    'Mercury': (3.3011e23, 0.39),
    'Venus': (4.8675e24, 0.72),
    'Earth': (5.97237e24, 1.00),
    'Mars': (6.4171e23, 1.52),
    'Jupiter': (1.8982e27, 5.20),
    'Saturn': (5.6834e26, 9.58),
    'Uranus': (8.6810e25, 19.20),
    'Neptune': (1.02413e26, 30.05)
}

# --- β range (degrees & radians) ---
beta_deg = np.linspace(0.01, 179.99, 1000)
beta_rad = np.deg2rad(beta_deg)

# --- Function: calculate deflection angles ---
def compute_deflection_angles(M, r):
    alpha_normal = (2 * G * M / (c**2 * r * np.sin(beta_rad))) * (1 + np.cos(beta_rad)) * rad2uas
    alpha_simplified = (4 * G * M / (c**2 * r * np.sin(beta_rad))) * rad2uas
    alpha_difference = alpha_simplified - alpha_normal
    return alpha_simplified, alpha_normal, alpha_difference

# --- Loop over the eight planets and plot ---
for planet, (M, a) in planet_data.items():
    # Set maximum/minimum Earth-planet distances
    r_min = abs(a - 1.0) * AU  # Minimum distance (conjunction)
    r_max = (a + 1.0) * AU     # Maximum distance (opposition)

    # Calculate deflection angles
    alpha_simplified_min, alpha_normal_min, alpha_diff_min = compute_deflection_angles(M, r_min)
    alpha_simplified_max, alpha_normal_max, alpha_diff_max = compute_deflection_angles(M, r_max)

    # Create subplots
    fig, axs = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

    for ax, alpha_simplified, alpha_normal, alpha_diff, title in zip(
        axs,
        [alpha_simplified_max, alpha_simplified_min],
        [alpha_normal_max, alpha_normal_min],
        [alpha_diff_max, alpha_diff_min],
        ['Maximum Distance', 'Minimum Distance']
    ):
        # Background heat zones
        ax.fill_between(beta_deg, 1e-5, 1e4, where=(alpha_simplified > 10), color='red', alpha=0.1)
        ax.fill_between(beta_deg, 1e-5, 1e4, where=((alpha_simplified <= 10) & (alpha_simplified > 1)), color='orange', alpha=0.1)
        ax.fill_between(beta_deg, 1e-5, 1e4, where=((alpha_simplified <= 1) & (alpha_simplified > 0.1)), color='yellow', alpha=0.1)

        # Plot curves
        ax.plot(beta_deg, alpha_simplified, 'r', label='Simplified')
        ax.plot(beta_deg, alpha_normal, 'g', label='Normal')
        ax.plot(beta_deg, alpha_diff, 'b', label='Difference')

        # Threshold lines
        ax.axhline(y=10.0, color='gray', linestyle=':', linewidth=1, label=r'$10\ \mu as$')
        ax.axhline(y=1.0, color='gray', linestyle='-', linewidth=1, label=r'$1.0\ \mu as$')
        ax.axhline(y=0.1, color='black', linestyle='--', linewidth=1, label=r'$0.1\ \mu as$')

        # Plot settings
        ax.set_yscale('log')
        ax.set_ylim(1e-5, 1e4)
        ax.set_xlabel(r'Angle between Planet and CES Seen from Earth, $\beta$ (°)')
        ax.set_title(f'{planet} - {title}')
        ax.grid(True, which="both", ls=':', linewidth=0.5)

    # y-axis label and legend
    axs[0].set_ylabel(r'Deflection Angle, $\alpha$ ($\mu$as)')
    axs[1].legend(loc='upper center')

    # Adjust layout and display
    plt.tight_layout()
    plt.show()
