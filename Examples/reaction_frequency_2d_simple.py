#!/usr/bin/env python3
"""
Simplified 2D Example: Reaction Frequency Tracking

This is a minimal 2D example for quick testing of position-dependent
reaction frequency tracking. It uses field ionization (ADK model) which
is simpler than MCC collisions and doesn't require cross-section files.

The setup creates a spatially-varying electric field that ionizes atoms
more in some regions than others, allowing us to visualize position-dependent
ionization rates.
"""

import numpy as np
from pywarpx import picmi

# ##############################################################################
# Simulation Parameters (2D for speed)
# ##############################################################################

# Reduced 2D domain
nx = 256
nz = 256

xmin = -10.0e-6  # -10 μm
xmax = 10.0e-6   # +10 μm
zmin = -10.0e-6
zmax = 10.0e-6

# Grid (2D)
grid = picmi.Cartesian2DGrid(
    number_of_cells=[nx, nz],
    lower_bound=[xmin, zmin],
    upper_bound=[xmax, zmax],
    lower_boundary_conditions=['periodic', 'periodic'],
    upper_boundary_conditions=['periodic', 'periodic'],
    warpx_max_grid_size=64
)

# Timestep and duration
dt = 1.0e-16  # 0.1 fs (small for field ionization)
max_steps = 200
diagnostic_interval = 20

# ##############################################################################
# Species Setup
# ##############################################################################

# Neutral helium atoms (will be ionized)
# Gaussian distribution centered in domain with varying density
helium_density = 1.0e24  # High density for visible ionization

helium_dist = picmi.AnalyticDistribution(
    density_expression=f"{helium_density} * exp(-((x)^2 + (z)^2) / (5.0e-6)^2)",
    upper_bound=[xmax, 0, zmax],
    lower_bound=[xmin, 0, zmin],
    directed_velocity=[0., 0., 0.],
    rms_velocity=[0., 0., 0.]  # Cold atoms
)

helium_atoms = picmi.Species(
    particle_type='He',
    name='he_atoms',
    charge_state=0,  # Neutral
    initial_distribution=helium_dist
)

# He+ ions (created by ionization)
helium_ions = picmi.Species(
    particle_type='He',
    name='he_ions',
    charge_state=+1,
    initial_distribution=None  # Created by ionization
)

# Electrons (created by ionization)
electrons = picmi.Species(
    particle_type='electron',
    name='electrons',
    initial_distribution=None  # Created by ionization
)

# ##############################################################################
# External Field for Ionization
# ##############################################################################

# Create a spatially-varying electric field using laser
# This will cause non-uniform ionization across space

# Gaussian laser pulse with focal spot
laser = picmi.GaussianLaser(
    wavelength=800.0e-9,  # 800 nm (Ti:Sapphire)
    waist=5.0e-6,  # 5 μm focal spot
    duration=50.0e-15,  # 50 fs pulse
    focal_position=[0, 0, 0],
    centroid_position=[0, 0, -5.0e-6],  # Starts behind focus
    propagation_direction=[0, 0, 1],
    polarization_direction=[1, 0, 0],
    a0=2.0,  # Normalized vector potential (strong field)
    fill_in=False
)

laser_antenna = picmi.LaserAntenna(
    position=[0., 0., zmin],
    normal_vector=[0, 0, 1]
)

# ##############################################################################
# Field Ionization (ADK Model)
# ##############################################################################

# ADK ionization for He -> He+ + e-
ionization = picmi.FieldIonization(
    model="ADK",
    ionized_species=helium_atoms,
    product_species=[helium_ions, electrons]
)

# ##############################################################################
# Diagnostics
# ##############################################################################

# Track product particles (ions and electrons from ionization)
particle_diag = picmi.ParticleDiagnostic(
    name='particle_diag',
    period=diagnostic_interval,
    species=[helium_atoms, helium_ions, electrons],
    data_list=['ux', 'uy', 'uz', 'w'],
    write_dir='diags_2d',
    warpx_format='openpmd',
    warpx_openpmd_backend='h5'
)

# Field diagnostic to see the laser
field_diag = picmi.FieldDiagnostic(
    name='field_diag',
    grid=grid,
    period=diagnostic_interval,
    data_list=['E', 'B', 'rho_electrons', 'rho_he_ions'],
    write_dir='diags_2d',
    warpx_format='openpmd',
    warpx_openpmd_backend='h5'
)

# ##############################################################################
# Simulation Assembly
# ##############################################################################

solver = picmi.ElectromagneticSolver(
    grid=grid,
    method='Yee',
    cfl=0.999
)

sim = picmi.Simulation(
    solver=solver,
    time_step_size=dt,
    max_steps=max_steps,
    verbose=1
)

# Add species
sim.add_species(
    helium_atoms,
    layout=picmi.GriddedLayout(
        n_macroparticle_per_cell=[4, 4],
        grid=grid
    )
)

sim.add_species(
    helium_ions,
    layout=picmi.GriddedLayout(
        n_macroparticle_per_cell=[4, 4],
        grid=grid
    )
)

sim.add_species(
    electrons,
    layout=picmi.GriddedLayout(
        n_macroparticle_per_cell=[4, 4],
        grid=grid
    )
)

# Add laser
sim.add_laser(laser, laser_antenna)

# Add ionization
sim.add_interaction(ionization)

# Add diagnostics
sim.add_diagnostic(particle_diag)
sim.add_diagnostic(field_diag)

# ##############################################################################
# Run
# ##############################################################################

if __name__ == '__main__':
    print("="*70)
    print("2D Field Ionization with Position-Dependent Reaction Tracking")
    print("="*70)
    print(f"\nSimulation parameters:")
    print(f"  Domain: {(xmax-xmin)*1e6:.1f} × {(zmax-zmin)*1e6:.1f} μm²")
    print(f"  Grid: {nx} × {nz}")
    print(f"  Timestep: {dt*1e15:.2f} fs")
    print(f"  Total steps: {max_steps}")
    print(f"  Total time: {dt*max_steps*1e15:.1f} fs")
    print("\nRunning simulation...")

    sim.step(max_steps)

    print("\n" + "="*70)
    print("Simulation complete!")
    print("="*70)
    print("\nTo analyze 2D reaction frequency, run:")
    print("  python analyze_reaction_frequency_2d.py")
    print("="*70)
