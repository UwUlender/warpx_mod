#!/usr/bin/env python3
"""
Analysis script for dusty plasma simulations in WarpX

This script analyzes the output from dusty plasma simulations and produces:
1. Dust charge evolution over time
2. Dust velocity distribution
3. Plasma-dust interaction statistics
4. Comparison with OML theory predictions

Usage:
    python analysis_dusty_plasma.py <path_to_plotfiles>
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from scipy.constants import e, m_e, m_p, epsilon_0, k as k_B
import sys
import os

try:
    import openpmd_api as io
except ImportError:
    print("Warning: openpmd_api not found. Trying legacy reader...")
    import yt
    use_yt = True
else:
    use_yt = False


def calculate_oml_equilibrium_charge(a, n_e, T_e, n_i, T_i, m_i):
    """
    Calculate equilibrium dust charge using OML theory.

    Parameters:
    -----------
    a : float
        Dust radius (m)
    n_e, n_i : float
        Electron and ion densities (m^-3)
    T_e, T_i : float
        Electron and ion temperatures (eV)
    m_i : float
        Ion mass (kg)

    Returns:
    --------
    Q_dust : float
        Equilibrium charge (in units of elementary charge)
    phi_dust : float
        Surface potential (V)
    """
    # Convert temperatures to Joules
    T_e_J = T_e * e
    T_i_J = T_i * e

    # Iteratively solve for equilibrium potential
    # Initial guess: floating potential
    phi_float = -(T_e / 2.0) * np.log(2 * np.pi * m_e / m_i * (T_e / T_i))

    phi_dust = phi_float
    for _ in range(100):  # Newton iteration
        # Electron current
        I_e = -e * n_e * np.pi * a**2 * np.sqrt(8 * k_B * T_e_J / (np.pi * m_e)) * \
              np.exp(e * phi_dust / (k_B * T_e_J))

        # Ion current
        I_i = e * n_i * np.pi * a**2 * np.sqrt(8 * k_B * T_i_J / (np.pi * m_i)) * \
              (1 - e * phi_dust / (k_B * T_i_J))

        # Check convergence
        if abs(I_e + I_i) < 1e-30:
            break

        # Update (simplified Newton step)
        dphi = -0.01 * (I_e + I_i) / abs(I_e)
        phi_dust += dphi

    # Calculate charge
    Q_dust = 4 * np.pi * epsilon_0 * a * phi_dust / e

    return Q_dust, phi_dust


def analyze_dust_charging(plotfile_dir):
    """
    Analyze dust charging evolution from simulation data.
    """
    print(f"Analyzing dust charging from: {plotfile_dir}")

    # Get list of plotfiles
    plotfiles = sorted([f for f in os.listdir(plotfile_dir) if f.startswith('diag')])

    if len(plotfiles) == 0:
        print("No plotfiles found!")
        return

    times = []
    avg_dust_charge = []
    std_dust_charge = []
    max_dust_charge = []
    min_dust_charge = []

    for pf in plotfiles:
        pf_path = os.path.join(plotfile_dir, pf)

        try:
            if use_yt:
                ds = yt.load(pf_path)
                ad = ds.all_data()
                dust_charges = ad['dust', 'dust_charge'].v
            else:
                # Use openPMD API
                series = io.Series(pf_path, io.Access.read_only)
                it = list(series.iterations)[0]
                dust = series.iterations[it].particles["dust"]
                dust_charges = dust["dust_charge"][io.Mesh_Record_Component.SCALAR].load_chunk()
                series.flush()

            times.append(it * 1e-10)  # Assuming 0.1 ns time step
            avg_dust_charge.append(np.mean(dust_charges))
            std_dust_charge.append(np.std(dust_charges))
            max_dust_charge.append(np.max(dust_charges))
            min_dust_charge.append(np.min(dust_charges))

        except Exception as e:
            print(f"Error reading {pf}: {e}")
            continue

    times = np.array(times)
    avg_dust_charge = np.array(avg_dust_charge)
    std_dust_charge = np.array(std_dust_charge)
    max_dust_charge = np.array(max_dust_charge)
    min_dust_charge = np.array(min_dust_charge)

    # Plot charging evolution
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Average charge
    ax1.plot(times * 1e9, avg_dust_charge, 'b-', linewidth=2, label='Average')
    ax1.fill_between(times * 1e9,
                      avg_dust_charge - std_dust_charge,
                      avg_dust_charge + std_dust_charge,
                      alpha=0.3, color='b', label='±1 std dev')

    # OML theory prediction
    # Typical parameters: n_e = 1e16 m^-3, T_e = 3 eV, T_i = 0.03 eV, a = 1 micron
    Q_oml, phi_oml = calculate_oml_equilibrium_charge(
        a=1e-6, n_e=1e16, T_e=3.0, n_i=1e16, T_i=0.03, m_i=40*m_p
    )
    ax1.axhline(Q_oml, color='r', linestyle='--', linewidth=2,
                label=f'OML equilibrium: {Q_oml:.1f} e')

    ax1.set_xlabel('Time (ns)', fontsize=12)
    ax1.set_ylabel('Dust Charge (elementary charges)', fontsize=12)
    ax1.set_title('Dust Particle Charging Evolution', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    # Min/max range
    ax2.plot(times * 1e9, max_dust_charge, 'g-', linewidth=1.5, label='Maximum')
    ax2.plot(times * 1e9, min_dust_charge, 'orange', linewidth=1.5, label='Minimum')
    ax2.plot(times * 1e9, avg_dust_charge, 'b-', linewidth=2, label='Average')
    ax2.axhline(Q_oml, color='r', linestyle='--', linewidth=2, label='OML prediction')

    ax2.set_xlabel('Time (ns)', fontsize=12)
    ax2.set_ylabel('Dust Charge (elementary charges)', fontsize=12)
    ax2.set_title('Dust Charge Distribution Range', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('dust_charging_evolution.png', dpi=300, bbox_inches='tight')
    print("Saved: dust_charging_evolution.png")
    plt.close()

    # Summary statistics
    print("\n" + "="*60)
    print("DUST CHARGING ANALYSIS SUMMARY")
    print("="*60)
    print(f"Final average dust charge: {avg_dust_charge[-1]:.2f} e")
    print(f"Final std deviation:       {std_dust_charge[-1]:.2f} e")
    print(f"OML theory prediction:     {Q_oml:.2f} e")
    print(f"Surface potential (OML):   {phi_oml:.2f} V")
    print(f"Relative difference:       {abs(avg_dust_charge[-1] - Q_oml)/abs(Q_oml)*100:.1f}%")
    print("="*60)


def analyze_dust_dynamics(plotfile_dir):
    """
    Analyze dust particle dynamics (velocities, trajectories).
    """
    print(f"\nAnalyzing dust dynamics from: {plotfile_dir}")

    plotfiles = sorted([f for f in os.listdir(plotfile_dir) if f.startswith('diag')])

    if len(plotfiles) == 0:
        return

    # Analyze final state
    final_pf = os.path.join(plotfile_dir, plotfiles[-1])

    try:
        if use_yt:
            ds = yt.load(final_pf)
            ad = ds.all_data()
            dust_ux = ad['dust', 'ux'].v
            dust_uy = ad['dust', 'uy'].v
            dust_uz = ad['dust', 'uz'].v
        else:
            series = io.Series(final_pf, io.Access.read_only)
            it = list(series.iterations)[0]
            dust = series.iterations[it].particles["dust"]
            dust_ux = dust["momentum"]["x"][io.Mesh_Record_Component.SCALAR].load_chunk()
            dust_uy = dust["momentum"]["y"][io.Mesh_Record_Component.SCALAR].load_chunk()
            dust_uz = dust["momentum"]["z"][io.Mesh_Record_Component.SCALAR].load_chunk()
            series.flush()

        # Calculate velocity magnitudes (assuming normalized momentum)
        dust_v = np.sqrt(dust_ux**2 + dust_uy**2 + dust_uz**2)

        # Plot velocity distribution
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.hist(dust_v, bins=50, edgecolor='black', alpha=0.7)
        ax.set_xlabel('Dust Velocity (normalized)', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title('Dust Velocity Distribution (Final State)', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig('dust_velocity_distribution.png', dpi=300, bbox_inches='tight')
        print("Saved: dust_velocity_distribution.png")
        plt.close()

        print(f"Average dust speed: {np.mean(dust_v):.2e}")
        print(f"Max dust speed:     {np.max(dust_v):.2e}")

    except Exception as e:
        print(f"Error analyzing dust dynamics: {e}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python analysis_dusty_plasma.py <plotfile_directory>")
        print("Example: python analysis_dusty_plasma.py ./diags/")
        sys.exit(1)

    plotfile_dir = sys.argv[1]

    if not os.path.exists(plotfile_dir):
        print(f"Error: Directory {plotfile_dir} not found!")
        sys.exit(1)

    # Perform analyses
    analyze_dust_charging(plotfile_dir)
    analyze_dust_dynamics(plotfile_dir)

    print("\nAnalysis complete!")
