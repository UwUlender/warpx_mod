#!/usr/bin/env python3
"""
Example: Position-Dependent Reaction Frequency Tracking in WarpX

This example demonstrates how to:
1. Set up ionization reactions using PICMI
2. Track reaction frequency as a function of position
3. Visualize spatial distribution of reactions

The simulation uses MCC (Monte Carlo Collisions) ionization to create
electron-ion pairs from a background gas, then analyzes where reactions
occur most frequently.
"""

import numpy as np
from pywarpx import picmi

# Physical constants
constants = picmi.constants

# ##############################################################################
# Simulation Parameters
# ##############################################################################

# Domain setup
nx = 128
ny = 128
nz = 128

xmin = 0.0
xmax = 1.0e-3  # 1 mm
ymin = 0.0
ymax = 1.0e-3
zmin = 0.0
zmax = 1.0e-3

# Grid
grid = picmi.Cartesian3DGrid(
    number_of_cells=[nx, ny, nz],
    lower_bound=[xmin, ymin, zmin],
    upper_bound=[xmax, ymax, zmax],
    lower_boundary_conditions=['periodic', 'periodic', 'periodic'],
    upper_boundary_conditions=['periodic', 'periodic', 'periodic'],
    warpx_max_grid_size=64
)

# Timestep
dt = 1.0e-11  # 10 ps
max_steps = 100
diagnostic_interval = 10

# ##############################################################################
# Particle Species Setup
# ##############################################################################

# Electron distribution (seed electrons for ionization)
uniform_dist = picmi.UniformDistribution(
    density=1.0e16,  # 10^16 m^-3 seed electrons
    upper_bound=[xmax, ymax, zmax],
    lower_bound=[xmin, ymin, zmin],
    directed_velocity=[0., 0., 0.]
)

# Initial electrons
electrons = picmi.Species(
    particle_type='electron',
    name='electrons',
    initial_distribution=uniform_dist,
    warpx_save_particles_at_eb=0
)

# Ions created by ionization (initially empty)
ions = picmi.Species(
    particle_type='He',
    name='ions',
    charge_state=1,
    initial_distribution=None  # Created by ionization
)

# ##############################################################################
# Background Gas and MCC Collisions
# ##############################################################################

# Background helium gas density (spatially varying)
# Higher density in center creates more reactions there
background_density_expression = """
    1.0e20 * (1.0 + 2.0 * exp(-((x-5.0e-4)^2 + (y-5.0e-4)^2 + (z-5.0e-4)^2) / (2.0e-4)^2))
"""
# This creates a Gaussian density peak at the center

# MCC collision configuration
# Cross-section files would normally be provided; here we use simplified data
mcc_collisions = picmi.MCCCollisions(
    name='mcc_ionization',
    species=electrons,
    background_density=background_density_expression,
    background_temperature=300.0,  # 300 K background gas
    background_mass=constants.m_p * 4.0,  # Helium mass
    ndt=1,  # Apply collisions every timestep
    max_background_density=3.0e20,  # Maximum for normalization
    scattering_processes={
        'elastic': {
            # In real use, provide path to cross-section file
            # 'cross_section': 'he_elastic.dat'
        },
        'ionization': {
            # In real use, provide path to cross-section file
            # 'cross_section': 'he_ionization.dat',
            'energy': 24.6,  # Ionization energy of He in eV
            'species': ions  # Product species
        }
    }
)

# ##############################################################################
# Diagnostics Setup
# ##############################################################################

# Particle diagnostics - track all particles
particle_diag = picmi.ParticleDiagnostic(
    name='particle_diag',
    period=diagnostic_interval,
    species=[electrons, ions],
    data_list=['ux', 'uy', 'uz', 'w'],  # Momentum and weight
    write_dir='diags',
    warpx_format='openpmd',
    warpx_openpmd_backend='h5'
)

# Field diagnostics - optional, to see EM fields
field_diag = picmi.FieldDiagnostic(
    name='field_diag',
    grid=grid,
    period=diagnostic_interval,
    data_list=['E', 'B', 'J', 'rho'],
    write_dir='diags',
    warpx_format='openpmd',
    warpx_openpmd_backend='h5'
)

# ##############################################################################
# Simulation Setup
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
    warpx_collisions=[mcc_collisions],
    verbose=1
)

sim.add_species(
    electrons,
    layout=picmi.GriddedLayout(
        n_macroparticle_per_cell=[2, 2, 2],
        grid=grid
    )
)

sim.add_species(
    ions,
    layout=picmi.GriddedLayout(
        n_macroparticle_per_cell=[2, 2, 2],
        grid=grid
    )
)

sim.add_diagnostic(particle_diag)
sim.add_diagnostic(field_diag)

# ##############################################################################
# Run Simulation
# ##############################################################################

if __name__ == '__main__':
    sim.step(max_steps)

    print("\n" + "="*70)
    print("Simulation complete!")
    print("="*70)
    print("\nTo analyze reaction frequency by position, run:")
    print("  python analyze_reaction_frequency.py")
    print("\nThis will:")
    print("  1. Load ion particle data (ions only exist from ionization)")
    print("  2. Bin ions by creation position")
    print("  3. Generate spatial reaction frequency maps")
    print("  4. Create visualizations (2D slices and 3D volume)")
    print("="*70)
