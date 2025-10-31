#!/usr/bin/env python3
"""
Analysis Script: Position-Dependent Reaction Frequency

This script processes WarpX output to determine where reactions occur
most frequently in space. It works by:
1. Loading product particle data (ions created by ionization)
2. Binning particles by position to create spatial histograms
3. Calculating reaction rates per cell
4. Visualizing the results

The reaction frequency is inferred from the spatial distribution of
product particles, weighted by their creation weight.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import glob
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

DIAG_DIR = 'diags/particle_diag'
PRODUCT_SPECIES = 'ions'  # Ions are created by ionization reactions

# Spatial binning
NBINS_X = 64
NBINS_Y = 64
NBINS_Z = 64

# Domain bounds (should match simulation)
XMIN, XMAX = 0.0, 1.0e-3
YMIN, YMAX = 0.0, 1.0e-3
ZMIN, ZMAX = 0.0, 1.0e-3

# ##############################################################################
# Functions
# ##############################################################################

def find_latest_iteration(diag_dir):
    """Find the most recent diagnostic output iteration."""
    series = io.Series(f"{diag_dir}/openpmd_%T.h5", io.Access.read_only)
    iterations = list(series.iterations)
    series.close()
    return iterations[-1] if iterations else None


def load_particle_data(diag_dir, iteration, species):
    """
    Load particle position and weight data for a given species and iteration.

    Returns:
        dict with keys 'x', 'y', 'z', 'w' (weight) as numpy arrays
    """
    series = io.Series(f"{diag_dir}/openpmd_%T.h5", io.Access.read_only)

    if iteration not in series.iterations:
        series.close()
        raise ValueError(f"Iteration {iteration} not found")

    it = series.iterations[iteration]

    try:
        particles = it.particles[species]
    except KeyError:
        series.close()
        raise ValueError(f"Species '{species}' not found in iteration {iteration}")

    # Load position data
    x = particles['position']['x'].load_chunk()
    y = particles['position']['y'].load_chunk()
    z = particles['position']['z'].load_chunk()

    # Load weight (number of real particles per macroparticle)
    w = particles['weighting'][io.Mesh_Record_Component.SCALAR].load_chunk()

    series.flush()

    data = {
        'x': np.array(x),
        'y': np.array(y),
        'z': np.array(z),
        'w': np.array(w)
    }

    series.close()
    return data


def compute_reaction_frequency_3d(data, nbins_x, nbins_y, nbins_z,
                                   xrange, yrange, zrange):
    """
    Compute 3D histogram of reaction frequency by position.

    The histogram is weighted by particle weight, so each bin represents
    the total number of real particles (reactions) that occurred there.

    Returns:
        freq_3d: 3D array of reaction counts per bin
        edges: tuple of (x_edges, y_edges, z_edges)
    """
    x, y, z, w = data['x'], data['y'], data['z'], data['w']

    # Create 3D histogram weighted by particle weight
    freq_3d, edges = np.histogramdd(
        np.column_stack([x, y, z]),
        bins=[nbins_x, nbins_y, nbins_z],
        range=[xrange, yrange, zrange],
        weights=w
    )

    return freq_3d, edges


def compute_reaction_frequency_2d(data, nbins_a, nbins_b,
                                   arange, brange, plane='xy'):
    """
    Compute 2D histogram of reaction frequency in a given plane.

    Args:
        plane: 'xy', 'xz', or 'yz'

    Returns:
        freq_2d: 2D array of reaction counts per bin
        edges: tuple of (a_edges, b_edges)
    """
    w = data['w']

    if plane == 'xy':
        a, b = data['x'], data['y']
    elif plane == 'xz':
        a, b = data['x'], data['z']
    elif plane == 'yz':
        a, b = data['y'], data['z']
    else:
        raise ValueError(f"Invalid plane: {plane}")

    freq_2d, a_edges, b_edges = np.histogram2d(
        a, b,
        bins=[nbins_a, nbins_b],
        range=[arange, brange],
        weights=w
    )

    return freq_2d, (a_edges, b_edges)


def plot_2d_reaction_frequency(freq_2d, edges, plane='xy', title=None,
                                 output_file=None):
    """Create 2D visualization of reaction frequency."""
    a_edges, b_edges = edges
    a_centers = 0.5 * (a_edges[:-1] + a_edges[1:])
    b_centers = 0.5 * (b_edges[:-1] + b_edges[1:])

    fig, ax = plt.subplots(figsize=(10, 8))

    # Use log scale if there's significant dynamic range
    vmin = freq_2d[freq_2d > 0].min() if np.any(freq_2d > 0) else 1
    vmax = freq_2d.max()

    if vmax / vmin > 10:
        norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = None

    im = ax.pcolormesh(a_centers * 1e3, b_centers * 1e3, freq_2d.T,
                       shading='auto', cmap='hot', norm=norm)

    cbar = plt.colorbar(im, ax=ax, label='Reaction Count (weighted)')

    # Labels
    axis_labels = {'x': 'X', 'y': 'Y', 'z': 'Z'}
    ax.set_xlabel(f'{axis_labels[plane[0]]} Position (mm)')
    ax.set_ylabel(f'{axis_labels[plane[1]]} Position (mm)')

    if title is None:
        title = f'Reaction Frequency in {plane.upper()} Plane'
    ax.set_title(title)

    ax.set_aspect('equal')
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_file}")

    return fig, ax


def plot_1d_profiles(freq_3d, edges, output_file=None):
    """Create 1D profiles by integrating over two dimensions."""
    x_edges, y_edges, z_edges = edges
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])

    # Integrate over dimensions
    freq_x = freq_3d.sum(axis=(1, 2))
    freq_y = freq_3d.sum(axis=(0, 2))
    freq_z = freq_3d.sum(axis=(0, 1))

    fig, axes = plt.subplots(3, 1, figsize=(10, 10))

    axes[0].plot(x_centers * 1e3, freq_x, 'b-', linewidth=2)
    axes[0].set_xlabel('X Position (mm)')
    axes[0].set_ylabel('Integrated Reaction Count')
    axes[0].set_title('Reaction Frequency Profile (X)')
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(y_centers * 1e3, freq_y, 'r-', linewidth=2)
    axes[1].set_xlabel('Y Position (mm)')
    axes[1].set_ylabel('Integrated Reaction Count')
    axes[1].set_title('Reaction Frequency Profile (Y)')
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(z_centers * 1e3, freq_z, 'g-', linewidth=2)
    axes[2].set_xlabel('Z Position (mm)')
    axes[2].set_ylabel('Integrated Reaction Count')
    axes[2].set_title('Reaction Frequency Profile (Z)')
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_file}")

    return fig, axes


def analyze_temporal_evolution(diag_dir, species, iterations):
    """
    Track how reaction frequency evolves over time.

    Returns total reactions at each iteration.
    """
    total_reactions = []

    series = io.Series(f"{diag_dir}/openpmd_%T.h5", io.Access.read_only)

    for it_num in iterations:
        if it_num not in series.iterations:
            continue

        it = series.iterations[it_num]

        try:
            particles = it.particles[species]
            w = particles['weighting'][io.Mesh_Record_Component.SCALAR].load_chunk()
            series.flush()
            total_reactions.append(np.sum(w))
        except KeyError:
            total_reactions.append(0)

    series.close()
    return np.array(total_reactions)


# ##############################################################################
# Main Analysis
# ##############################################################################

def main():
    print("="*70)
    print("WarpX Reaction Frequency Analysis")
    print("="*70)

    # Check if diagnostic directory exists
    if not os.path.exists(DIAG_DIR):
        print(f"\nError: Diagnostic directory not found: {DIAG_DIR}")
        print("Please run the simulation first:")
        print("  python reaction_frequency_tracking_example.py")
        return

    # Find latest iteration
    print(f"\nSearching for data in: {DIAG_DIR}")
    iteration = find_latest_iteration(DIAG_DIR)

    if iteration is None:
        print("Error: No iterations found in diagnostic output")
        return

    print(f"Latest iteration: {iteration}")

    # Load product particle data
    print(f"\nLoading {PRODUCT_SPECIES} particle data...")
    try:
        data = load_particle_data(DIAG_DIR, iteration, PRODUCT_SPECIES)
    except Exception as e:
        print(f"Error loading data: {e}")
        return

    n_particles = len(data['x'])
    total_weight = np.sum(data['w'])

    print(f"  Macroparticles: {n_particles:,}")
    print(f"  Total weighted count: {total_weight:.3e}")
    print(f"  (This represents the total number of reactions)")

    # Compute 3D reaction frequency
    print("\nComputing 3D reaction frequency histogram...")
    freq_3d, edges_3d = compute_reaction_frequency_3d(
        data, NBINS_X, NBINS_Y, NBINS_Z,
        [XMIN, XMAX], [YMIN, YMAX], [ZMIN, ZMAX]
    )

    print(f"  Grid: {NBINS_X} × {NBINS_Y} × {NBINS_Z}")
    print(f"  Max reactions in single cell: {freq_3d.max():.3e}")
    print(f"  Min reactions in single cell: {freq_3d.min():.3e}")

    # Create output directory
    os.makedirs('reaction_analysis', exist_ok=True)

    # Plot 2D slices
    print("\nGenerating 2D reaction frequency maps...")

    for plane in ['xy', 'xz', 'yz']:
        print(f"  {plane.upper()} plane...")

        if plane == 'xy':
            freq_2d, edges_2d = compute_reaction_frequency_2d(
                data, NBINS_X, NBINS_Y,
                [XMIN, XMAX], [YMIN, YMAX], plane
            )
        elif plane == 'xz':
            freq_2d, edges_2d = compute_reaction_frequency_2d(
                data, NBINS_X, NBINS_Z,
                [XMIN, XMAX], [ZMIN, ZMAX], plane
            )
        else:  # yz
            freq_2d, edges_2d = compute_reaction_frequency_2d(
                data, NBINS_Y, NBINS_Z,
                [YMIN, YMAX], [ZMIN, ZMAX], plane
            )

        plot_2d_reaction_frequency(
            freq_2d, edges_2d, plane,
            output_file=f'reaction_analysis/reaction_freq_{plane}.png'
        )

    # Plot 1D profiles
    print("\nGenerating 1D reaction profiles...")
    plot_1d_profiles(
        freq_3d, edges_3d,
        output_file='reaction_analysis/reaction_profiles_1d.png'
    )

    # Summary statistics
    print("\n" + "="*70)
    print("Summary Statistics")
    print("="*70)

    # Find peak reaction location
    peak_idx = np.unravel_index(freq_3d.argmax(), freq_3d.shape)
    x_edges, y_edges, z_edges = edges_3d

    peak_x = 0.5 * (x_edges[peak_idx[0]] + x_edges[peak_idx[0] + 1])
    peak_y = 0.5 * (y_edges[peak_idx[1]] + y_edges[peak_idx[1] + 1])
    peak_z = 0.5 * (z_edges[peak_idx[2]] + z_edges[peak_idx[2] + 1])

    print(f"\nPeak reaction location:")
    print(f"  X = {peak_x*1e3:.4f} mm")
    print(f"  Y = {peak_y*1e3:.4f} mm")
    print(f"  Z = {peak_z*1e3:.4f} mm")
    print(f"  Reactions in peak cell: {freq_3d.max():.3e}")

    # Cell volumes
    dx = (XMAX - XMIN) / NBINS_X
    dy = (YMAX - YMIN) / NBINS_Y
    dz = (ZMAX - ZMIN) / NBINS_Z
    cell_volume = dx * dy * dz

    print(f"\nReaction density (reactions per m³):")
    print(f"  Peak: {freq_3d.max() / cell_volume:.3e} m⁻³")
    print(f"  Average: {freq_3d.mean() / cell_volume:.3e} m⁻³")

    # Save data to file
    print("\nSaving processed data...")
    np.savez('reaction_analysis/reaction_frequency_data.npz',
             freq_3d=freq_3d,
             x_edges=edges_3d[0],
             y_edges=edges_3d[1],
             z_edges=edges_3d[2],
             total_reactions=total_weight,
             iteration=iteration)
    print("  Saved: reaction_analysis/reaction_frequency_data.npz")

    print("\n" + "="*70)
    print("Analysis complete!")
    print("="*70)
    print("\nOutput files:")
    print("  reaction_analysis/reaction_freq_xy.png")
    print("  reaction_analysis/reaction_freq_xz.png")
    print("  reaction_analysis/reaction_freq_yz.png")
    print("  reaction_analysis/reaction_profiles_1d.png")
    print("  reaction_analysis/reaction_frequency_data.npz")
    print("="*70)


if __name__ == '__main__':
    main()
