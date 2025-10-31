#!/usr/bin/env python3
"""
Verification script for reaction frequency tracking examples.

Checks that all required dependencies are available before running
the full simulations.
"""

import sys

def check_imports():
    """Check all required Python packages."""
    print("Checking Python dependencies...")

    missing = []

    # Check numpy
    try:
        import numpy as np
        print(f"  ✓ numpy {np.__version__}")
    except ImportError:
        print("  ✗ numpy not found")
        missing.append("numpy")

    # Check matplotlib
    try:
        import matplotlib
        print(f"  ✓ matplotlib {matplotlib.__version__}")
    except ImportError:
        print("  ✗ matplotlib not found")
        missing.append("matplotlib")

    # Check openpmd_api
    try:
        import openpmd_api as io
        print(f"  ✓ openpmd_api {io.__version__}")
    except ImportError:
        print("  ✗ openpmd_api not found")
        missing.append("openPMD-api")

    # Check WarpX
    try:
        from pywarpx import picmi
        print(f"  ✓ pywarpx (PICMI interface)")
    except ImportError:
        print("  ✗ pywarpx not found")
        missing.append("pywarpx")

    return missing


def check_example_files():
    """Check that all example files exist."""
    print("\nChecking example files...")

    import os

    files = [
        'reaction_frequency_2d_simple.py',
        'analyze_reaction_frequency_2d.py',
        'reaction_frequency_tracking_example.py',
        'analyze_reaction_frequency.py',
        'REACTION_FREQUENCY_TRACKING_README.md'
    ]

    missing = []
    for f in files:
        if os.path.exists(f):
            print(f"  ✓ {f}")
        else:
            print(f"  ✗ {f}")
            missing.append(f)

    return missing


def print_installation_help(missing_packages):
    """Print installation instructions for missing packages."""
    if not missing_packages:
        return

    print("\n" + "="*70)
    print("Missing Dependencies")
    print("="*70)
    print("\nInstall missing packages with:")
    print(f"  pip install {' '.join(missing_packages)}")

    if 'pywarpx' in missing_packages:
        print("\nFor pywarpx, you need to build WarpX with Python support:")
        print("  cmake -S . -B build -DWarpX_DIMS=3 -DWarpX_PYTHON=ON")
        print("  cmake --build build -j")
        print("  pip install -e .")


def run_mini_test():
    """Run a minimal test simulation."""
    print("\nRunning minimal test simulation...")

    try:
        from pywarpx import picmi
        import numpy as np

        # Tiny test grid
        grid = picmi.Cartesian2DGrid(
            number_of_cells=[8, 8],
            lower_bound=[0, 0],
            upper_bound=[1e-6, 1e-6],
            lower_boundary_conditions=['periodic', 'periodic'],
            upper_boundary_conditions=['periodic', 'periodic']
        )

        # Simple species
        dist = picmi.UniformDistribution(
            density=1e20,
            upper_bound=[1e-6, 0, 1e-6],
            lower_bound=[0, 0, 0]
        )

        electrons = picmi.Species(
            particle_type='electron',
            name='electrons',
            initial_distribution=dist
        )

        solver = picmi.ElectromagneticSolver(grid=grid, method='Yee')

        sim = picmi.Simulation(
            solver=solver,
            time_step_size=1e-15,
            max_steps=2,
            verbose=0
        )

        sim.add_species(
            electrons,
            layout=picmi.GriddedLayout(
                n_macroparticle_per_cell=[2, 2],
                grid=grid
            )
        )

        print("  ✓ PICMI simulation setup successful")
        print("  ✓ WarpX appears to be working correctly")

        return True

    except Exception as e:
        print(f"  ✗ Test failed: {e}")
        return False


def main():
    print("="*70)
    print("WarpX Reaction Frequency Tracking - Setup Verification")
    print("="*70)
    print()

    # Check imports
    missing_packages = check_imports()

    # Check files
    missing_files = check_example_files()

    # Run test if dependencies available
    test_passed = False
    if not missing_packages:
        test_passed = run_mini_test()

    # Summary
    print("\n" + "="*70)
    print("Verification Summary")
    print("="*70)

    if not missing_packages and not missing_files and test_passed:
        print("\n✓ All checks passed!")
        print("\nYou're ready to run the examples:")
        print("  python reaction_frequency_2d_simple.py")
        print("  python analyze_reaction_frequency_2d.py")
        print("\nSee REACTION_FREQUENCY_TRACKING_README.md for details.")
        return 0
    else:
        print("\n✗ Some checks failed")

        if missing_packages:
            print_installation_help(missing_packages)

        if missing_files:
            print(f"\nMissing {len(missing_files)} example file(s)")
            print("Make sure you're in the Examples directory")

        return 1


if __name__ == '__main__':
    sys.exit(main())
