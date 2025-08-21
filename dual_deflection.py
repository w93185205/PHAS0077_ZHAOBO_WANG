from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Iterable, Tuple

# ------------------------
# Physical constants and units
# ------------------------
G = 6.67430e-11            # m^3 kg^-1 s^-2
c = 299792458.0            # m s^-1
AU = 1.495978707e11        # m
RAD_TO_UAS = (180.0/np.pi)*3600.0*1e6  # rad -> microarcsec

# ------------------------
# Planetary masses (kg)
# ------------------------
PLANET_MASS = {
    "Jupiter": 1.89813e27,
    "Saturn": 5.68319e26,
    "Uranus": 8.6810e25,
    "Neptune": 1.02413e26,
}

@dataclass
class BodyGeom:
    name: str
    mass: float
    r_obs: float       # [m]
    beta: float        # [rad]
    n_hat: np.ndarray  # unit vector (3D) in tangent plane

# ------------------------
# Utility functions
# ------------------------
def unit(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    return v if n == 0.0 else v / n

def tangent_direction(l_hat: np.ndarray, p_hat: np.ndarray) -> np.ndarray:
    """
    Tangent-plane direction:
      t = p_hat - (p_hat·l_hat) l_hat
      n_hat = t / ||t||
    """
    t = p_hat - np.dot(p_hat, l_hat) * l_hat
    norm = np.linalg.norm(t)
    if norm == 0.0:
        # Degenerate case: construct any direction orthogonal to l_hat
        v = np.array([1.0, 0.0, 0.0], dtype=float)
        if abs(np.dot(v, l_hat)) > 0.9:
            v = np.array([0.0, 1.0, 0.0], dtype=float)
        t = v - np.dot(v, l_hat) * l_hat
        t = t / np.linalg.norm(t)
        return t
    return t / norm

def alpha_general(M: float, r: float, beta: float) -> float:
    """
    General single-body deflection:
        α(β) = [2 G M / (c^2 r sinβ)] (1 + cosβ)
    Return: rad
    """
    sb = np.sin(beta)
    if sb == 0.0:
        return np.inf
    return (2.0 * G * M / (c**2 * r * sb)) * (1.0 + np.cos(beta))

def make_body_by_rb(name: str, mass: float, r_au: float, beta_deg: float,
                    l_hat: np.ndarray | None = None,
                    p_hat: np.ndarray | None = None) -> BodyGeom:
    """
    Construct BodyGeom using (r[AU], beta[deg]). If l_hat/p_hat not given, set l_hat=+z,
    and p_hat lying in x–z plane forming angle beta with l_hat.
    """
    r = r_au * AU
    beta = np.deg2rad(beta_deg)
    if (l_hat is None) or (p_hat is None):
        l_hat = np.array([0.0, 0.0, 1.0], dtype=float)
        p_hat = unit(np.array([np.sin(beta), 0.0, np.cos(beta)], dtype=float))
    n_hat = tangent_direction(l_hat, p_hat)
    return BodyGeom(name=name, mass=mass, r_obs=r, beta=beta, n_hat=n_hat)

def dual_deflection_vector(
    M1: float, r1: float, beta1: float, n1: np.ndarray,
    M2: float, r2: float, beta2: float, n2: np.ndarray
) -> tuple[np.ndarray, float, float, float]:
    """
    Dual-body vector superposition:
      α_vec_tot = α1 n1 + α2 n2
    Return: (α_vec_tot[3], α1, α2, Δφ)
    """
    a1 = alpha_general(M1, r1, beta1)
    a2 = alpha_general(M2, r2, beta2)
    cos_dphi = float(np.clip(np.dot(n1, n2), -1.0, 1.0))
    dphi = float(np.arccos(cos_dphi))
    a_vec = a1 * n1 + a2 * n2
    return a_vec, a1, a2, dphi

# ------------------------
# Sweep and output (tidy format without NaNs)
# ------------------------
def sweep_pair_clean(
    pair_name: str,
    body1_name: str, mass1: float, r1_au: float,
    body2_name: str, mass2: float, r2_au: float,
    beta_deg_list: Iterable[float],
    delta_phi_deg: float = 30.0
) -> pd.DataFrame:
    """
    Return a tidy table without NaNs:
    Columns: pair, beta_deg, body1, body2, alpha_1_uas, alpha_2_uas, alpha_total_uas, delta_phi_deg
    """
    rows = []
    l_hat = np.array([0.0, 0.0, 1.0], dtype=float)
    rot = np.deg2rad(delta_phi_deg)

    for beta_deg in beta_deg_list:
        b1 = make_body_by_rb(body1_name, mass1, r1_au, beta_deg, l_hat=l_hat)
        beta = b1.beta
        # Rotate second body’s direction in tangent plane around l_hat by rot to form Δφ
        p2 = unit(np.array([np.sin(beta)*np.cos(rot), np.sin(beta)*np.sin(rot), np.cos(beta)]))
        n2 = tangent_direction(l_hat, p2)
        b2 = BodyGeom(name=body2_name, mass=mass2, r_obs=r2_au*AU, beta=beta, n_hat=n2)

        a_vec, a1, a2, dphi = dual_deflection_vector(
            b1.mass, b1.r_obs, b1.beta, b1.n_hat,
            b2.mass, b2.r_obs, b2.beta, b2.n_hat
        )

        rows.append({
            "pair": pair_name,
            "beta_deg": float(beta_deg),
            "body1": body1_name,
            "body2": body2_name,
            "alpha_1_uas": a1 * RAD_TO_UAS,
            "alpha_2_uas": a2 * RAD_TO_UAS,
            "alpha_total_uas": np.linalg.norm(a_vec) * RAD_TO_UAS,
            "delta_phi_deg": np.rad2deg(dphi),
        })
    return pd.DataFrame(rows)

def save_pair_wide(df_clean: pd.DataFrame, out_csv: str) -> pd.DataFrame:
    """
    Convert tidy table into a “wide” table (with dedicated body1/body2 columns), only for a single pair.
    Return wide table and save CSV.
    """
    assert df_clean["pair"].nunique() == 1, "df_clean should only contain a single pair"
    b1 = df_clean["body1"].iloc[0]
    b2 = df_clean["body2"].iloc[0]
    wide = df_clean[["beta_deg", "alpha_1_uas", "alpha_2_uas", "alpha_total_uas"]].copy()
    wide = wide.rename(columns={
        "alpha_1_uas": f"alpha_{b1}_uas",
        "alpha_2_uas": f"alpha_{b2}_uas",
    })
    wide.to_csv(out_csv, index=False)
    return wide

def plot_total_vs_beta(df_all_clean: pd.DataFrame, out_png: str) -> None:
    """
    Plot total deflection vs beta for different pairs (log-log scale)
    """
    plt.figure()
    for pair_name, g in df_all_clean.groupby("pair"):
        g = g.sort_values("beta_deg")
        plt.plot(g["beta_deg"], g["alpha_total_uas"], marker="o", label=pair_name)
    plt.xscale("log"); plt.yscale("log")
    plt.xlabel(r"$\beta$ (deg)")
    plt.ylabel(r"$\alpha_{\rm total}$ (µas)")
    plt.title("Dual-body total deflection vs beta")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=220)
    plt.close()

def plot_components(df_clean: pd.DataFrame, out_png: str) -> None:
    """
    For a single pair, plot body1/body2/total comparison curves (log-log scale)
    """
    assert df_clean["pair"].nunique() == 1, "df_clean should only contain a single pair"
    b1 = df_clean["body1"].iloc[0]
    b2 = df_clean["body2"].iloc[0]
    g = df_clean.sort_values("beta_deg")
    plt.figure()
    plt.plot(g["beta_deg"], g["alpha_1_uas"], marker="o", linestyle="--", label=b1)
    plt.plot(g["beta_deg"], g["alpha_2_uas"], marker="o", linestyle="--", label=b2)
    plt.plot(g["beta_deg"], g["alpha_total_uas"], marker="o", label="Total")
    plt.xscale("log"); plt.yscale("log")
    plt.xlabel(r"$\beta$ (deg)")
    plt.ylabel(r"$\alpha$ (µas)")
    plt.title(f"{b1}–{b2} components vs total")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=220)
    plt.close()

# ------------------------
# Main program
# ------------------------
if __name__ == "__main__":
    # Sweep angles (avoid 0° singularity)
    beta_list = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 45, 89]

    # Typical Earth–planet distances (AU) (can be replaced with r_min / r_max)
    r_au = {
        "Jupiter": 5.2,
        "Saturn": 9.5,
        "Uranus": 19.2,
        "Neptune": 30.1,
    }

    # Jupiter–Saturn (tidy table without NaNs)
    df_js_clean = sweep_pair_clean(
        pair_name="Jupiter–Saturn",
        body1_name="Jupiter", mass1=PLANET_MASS["Jupiter"], r1_au=r_au["Jupiter"],
        body2_name="Saturn",  mass2=PLANET_MASS["Saturn"],  r2_au=r_au["Saturn"],
        beta_deg_list=beta_list,
        delta_phi_deg=30.0
    )

    # Uranus–Neptune (tidy table without NaNs)
    df_un_clean = sweep_pair_clean(
        pair_name="Uranus–Neptune",
        body1_name="Uranus",  mass1=PLANET_MASS["Uranus"],  r1_au=r_au["Uranus"],
        body2_name="Neptune", mass2=PLANET_MASS["Neptune"], r2_au=r_au["Neptune"],
        beta_deg_list=beta_list,
        delta_phi_deg=30.0
    )

    # Combined tidy table (no NaNs)
    df_all_clean = pd.concat([df_js_clean, df_un_clean], ignore_index=True)
    df_all_clean.to_csv("dual_deflection_all.csv", index=False)
    print("Saved: dual_deflection_all.csv (tidy table without NaNs)")

    # Wide tables (dedicated columns for each pair)
    js_wide = save_pair_wide(df_js_clean, "dual_deflection_JS.csv")
    un_wide = save_pair_wide(df_un_clean, "dual_deflection_UN.csv")
    print("Saved: dual_deflection_JS.csv, dual_deflection_UN.csv (wide tables)")

    plot_total_vs_beta(df_all_clean, "Dual-body_total_deflection_vs_beta.png")
    plot_components(df_js_clean, "Jupiter–Saturn_components_vs_total.png")
    plot_components(df_un_clean, "Uranus–Neptune_components_vs_total.png")
    print("Saved: Dual-body_total_deflection_vs_beta.png, Jupiter–Saturn_components_vs_total.png, Uranus–Neptune_components_vs_total.png")

    with pd.option_context("display.max_columns", None, "display.width", 140):
        print(df_all_clean.head(12))
        print("\nJupiter–Saturn wide table:\n", js_wide.head(10))
        print("\nUranus–Neptune wide table:\n", un_wide.head(10))
