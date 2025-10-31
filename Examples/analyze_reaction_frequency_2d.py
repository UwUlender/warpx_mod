#!/usr/bin/env python3
"""
Analysis Script for 2D Reaction Frequency

Analyzes the 2D ionization simulation to show where reactions
(ionization events) occur most frequently in space.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os

try:
    import openpmd_api as io
except ImportError:
    print("Error: openpmd_api not installed")
    print("Install with: pip install openPMD-api")
    exit(1)

# ##############################################################################
# Configuration
# ##############################################################################

DIAG_DIR = 'diags_2d/particle_diag'
ION_SPECIES = 'he_ions'  # Product species from ionization

# Spatial binning
NBINS_X = 128
NBINS_Z = 128

# Domain bounds (should match simulation)
XMIN, XMAX = -10.0e-6, 10.0e-6
ZMIN, ZMAX = -10.0e-6, 10.0e-6

# ##############################################################################
# Helper Functions
# ##############################################################################

def find_iterations(diag_dir):
    """Find all available iterations."""
    series = io.Series(f"{diag_dir}/openpmd_%T.h5", io.Access.read_only)
    iterations = sorted(list(series.iterations))
    series.close()
    return iterations


def load_particle_data_2d(diag_dir, iteration, species):
    """Load 2D particle data (x, z, weight)."""
    series = io.Series(f"{diag_dir}/openpmd_%T.h5", io.Access.read_only)

    if iteration not in series.iterations:
        series.close()
        raise ValueError(f"Iteration {iteration} not found")

    it = series.iterations[iteration]

    try:
        particles = it.particles[species]
    except KeyError:
        series.close()
        return None  # Species doesn't exist yet

    # Load data
    x = particles['position']['x'].load_chunk()
    z = particles['position']['z'].load_chunk()
    w = particles['weighting'][io.Mesh_Record_Component.SCALAR].load_chunk()

    series.flush()
    series.close()

    return {
        'x': np.array(x),
        'z': np.array(z),
        'w': np.array(w)
    }


def compute_2d_histogram(data, nbins_x, nbins_z, xrange, zrange):
    """Compute 2D weighted histogram."""
    if data is None or len(data['x']) == 0:
        return np.zeros((nbins_x, nbins_z)), None

    hist, x_edges, z_edges = np.histogram2d(
        data['x'], data['z'],
        bins=[nbins_x, nbins_z],
        range=[xrange, zrange],
        weights=data['w']
    )

    return hist, (x_edges, z_edges)


def plot_2d_map(hist, edges, title, output_file=None, vmax=None):
    """Create 2D visualization."""
    if edges is None:
        print(f"Skipping plot: {title} (no data)")
        return

    x_edges, z_edges = edges
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])

    fig, ax = plt.subplots(figsize=(10, 9))

    # Handle log scale
    hist_plot = hist.copy()
    hist_plot[hist_plot == 0] = np.nan  # Don't show zeros in log scale

    vmin = np.nanmin(hist_plot) if np.any(~np.isnan(hist_plot)) else 1
    if vmax is None:
        vmax = np.nanmax(hist_plot) if np.any(~np.isnan(hist_plot)) else 1

    if vmax / vmin > 10:
        norm = LogNorm(vmin=max(vmin, vmax/1e6), vmax=vmax)
    else:
        norm = None

    im = ax.pcolormesh(x_centers * 1e6, z_centers * 1e6, hist_plot.T,
                       shading='auto', cmap='hot', norm=norm)

    cbar = plt.colorbar(im, ax=ax, label='Reaction Count (weighted)')

    ax.set_xlabel('X Position (μm)')
    ax.set_ylabel('Z Position (μm)')
    ax.set_title(title)
    ax.set_aspect('equal')

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"  Saved: {output_file}")

    return fig, ax


# ##############################################################################
# Main Analysis
# ##############################################################################

def main():
    print("="*70)
    print("2D Reaction Frequency Analysis")
    print("="*70)

    if not os.path.exists(DIAG_DIR):
        print(f"\nError: Directory not found: {DIAG_DIR}")
        print("Run simulation first: python reaction_frequency_2d_simple.py")
        return

    # Find all iterations
    print(f"\nSearching for data in: {DIAG_DIR}")
    iterations = find_iterations(DIAG_DIR)

    if not iterations:
        print("Error: No data found")
        return

    print(f"Found {len(iterations)} iterations: {iterations}")

    # Create output directory
    os.makedirs('reaction_analysis_2d', exist_ok=True)

    # Analyze each iteration
    print("\nAnalyzing reaction frequency at each time step...")

    total_reactions = []
    all_histograms = []

    for it in iterations:
        data = load_particle_data_2d(DIAG_DIR, it, ION_SPECIES)

        if data is None or len(data['x']) == 0:
            print(f"  Iteration {it}: No ions yet")
            total_reactions.append(0)
            all_histograms.append(None)
            continue

        n_particles = len(data['x'])
        total_weight = np.sum(data['w'])
        total_reactions.append(total_weight)

        print(f"  Iteration {it}: {n_particles:,} macroparticles, "
              f"{total_weight:.3e} total ions")

        # Compute histogram
        hist, edges = compute_2d_histogram(
            data, NBINS_X, NBINS_Z,
            [XMIN, XMAX], [ZMIN, ZMAX]
        )
        all_histograms.append((hist, edges))

    # Plot time evolution
    print("\nGenerating temporal evolution plot...")
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(iterations, total_reactions, 'b-o', linewidth=2, markersize=6)
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Total Ionization Events (cumulative)')
    ax.set_title('Cumulative Ionization vs Time')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('reaction_analysis_2d/ionization_vs_time.png', dpi=150)
    print("  Saved: reaction_analysis_2d/ionization_vs_time.png")

    # Plot spatial distributions at key time points
    print("\nGenerating spatial distribution plots...")

    # Find indices for early, middle, and late times
    n_its = len(iterations)
    key_indices = [
        (n_its // 4, 'early'),
        (n_its // 2, 'middle'),
        (3 * n_its // 4, 'late'),
        (-1, 'final')
    ]

    # Get max value for consistent color scale
    vmax = max(h[0].max() if h is not None else 0 for h in all_histograms)

    for idx, label in key_indices:
        if idx >= len(all_histograms):
            continue

        hist_data = all_histograms[idx]
        if hist_data is None:
            continue

        hist, edges = hist_data
        iteration = iterations[idx]

        plot_2d_map(
            hist, edges,
            f'Ionization Frequency (iteration {iteration}, {label} time)',
            output_file=f'reaction_analysis_2d/ionization_map_{label}.png',
            vmax=vmax
        )

    # Create animation data if multiple time steps
    if len([h for h in all_histograms if h is not None]) > 1:
        print("\nSaving animation data...")
        np.savez('reaction_analysis_2d/animation_data.npz',
                 iterations=np.array(iterations),
                 total_reactions=np.array(total_reactions),
                 histograms=[h[0] if h is not None else None
                            for h in all_histograms],
                 x_edges=all_histograms[-1][1][0] if all_histograms[-1] else None,
                 z_edges=all_histograms[-1][1][1] if all_histograms[-1] else None)
        print("  Saved: reaction_analysis_2d/animation_data.npz")

    # Summary statistics
    print("\n" + "="*70)
    print("Summary Statistics")
    print("="*70)

    final_total = total_reactions[-1]
    print(f"\nTotal ionization events: {final_total:.3e}")

    if all_histograms[-1] is not None:
        final_hist, final_edges = all_histograms[-1]
        peak_idx = np.unravel_index(final_hist.argmax(), final_hist.shape)

        x_edges, z_edges = final_edges
        peak_x = 0.5 * (x_edges[peak_idx[0]] + x_edges[peak_idx[0] + 1])
        peak_z = 0.5 * (z_edges[peak_idx[1]] + z_edges[peak_idx[1] + 1])

        print(f"\nPeak ionization location:")
        print(f"  X = {peak_x*1e6:.2f} μm")
        print(f"  Z = {peak_z*1e6:.2f} μm")
        print(f"  Peak cell count: {final_hist.max():.3e}")

        # Cell area
        dx = (XMAX - XMIN) / NBINS_X
        dz = (ZMAX - ZMIN) / NBINS_Z
        cell_area = dx * dz

        print(f"\nIonization density:")
        print(f"  Peak: {final_hist.max() / cell_area:.3e} m⁻²")
        print(f"  Average: {final_hist.mean() / cell_area:.3e} m⁻²")

    print("\n" + "="*70)
    print("Analysis complete!")
    print("="*70)
    print("\nOutput files:")
    print("  reaction_analysis_2d/ionization_vs_time.png")
    print("  reaction_analysis_2d/ionization_map_*.png")
    print("  reaction_analysis_2d/animation_data.npz")
    print("="*70)


if __name__ == '__main__':
    main()
