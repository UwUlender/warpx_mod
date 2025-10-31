# Position-Dependent Reaction Frequency Tracking in WarpX

This directory contains examples demonstrating how to track and visualize the spatial distribution of reactions (collisions, ionization, fusion, etc.) in WarpX using the Python/PICMI interface.

## Overview

WarpX supports various particle reactions:
- **Field Ionization** (ADK model)
- **MCC Collisions** (Monte Carlo Collisions with background gas)
- **Nuclear Fusion** (D-D, D-T, p-B, etc.)
- **Coulomb Collisions**
- **QED Processes** (photon emission, pair generation, Schwinger)
- **Bremsstrahlung**

While WarpX doesn't have a built-in diagnostic for spatially-resolved reaction rates, we can track reaction frequency by position through post-processing of product particle data.

## Methodology

### How It Works

1. **Simulation Setup**: Configure a reaction process (e.g., ionization, fusion) using PICMI
2. **Particle Output**: Output product particles that are created by reactions
3. **Post-Processing**: Analyze product particle positions and weights to determine where reactions occurred
4. **Visualization**: Create spatial histograms showing reaction frequency vs. position

### Key Principle

Product particles only exist because of reactions. By binning these particles by position (weighted by their particle weight), we create a map of where reactions occurred most frequently.

## Example Files

### 1. Simple 2D Example (Recommended Starting Point)

**Simulation**: `reaction_frequency_2d_simple.py`
- 2D field ionization simulation using ADK model
- Laser pulse ionizes helium atoms
- Fast runtime (~minutes on typical hardware)
- No external data files required

**Analysis**: `analyze_reaction_frequency_2d.py`
- Processes 2D simulation output
- Creates spatial ionization frequency maps
- Shows temporal evolution

**Usage**:
```bash
# Run simulation
python reaction_frequency_2d_simple.py

# Analyze results
python analyze_reaction_frequency_2d.py
```

**Output**:
- `reaction_analysis_2d/ionization_map_*.png` - Spatial distribution at different times
- `reaction_analysis_2d/ionization_vs_time.png` - Temporal evolution
- `reaction_analysis_2d/animation_data.npz` - Data for creating animations

### 2. Full 3D MCC Example

**Simulation**: `reaction_frequency_tracking_example.py`
- 3D Monte Carlo Collision (MCC) simulation
- Background gas ionization with spatially-varying density
- More realistic but slower to run

**Analysis**: `analyze_reaction_frequency.py`
- Full 3D spatial analysis
- Multiple 2D slice visualizations (XY, XZ, YZ planes)
- 1D integrated profiles along each axis

**Usage**:
```bash
# Run simulation
python reaction_frequency_tracking_example.py

# Analyze results
python analyze_reaction_frequency.py
```

**Output**:
- `reaction_analysis/reaction_freq_xy.png` - XY plane view
- `reaction_analysis/reaction_freq_xz.png` - XZ plane view
- `reaction_analysis/reaction_freq_yz.png` - YZ plane view
- `reaction_analysis/reaction_profiles_1d.png` - 1D profiles
- `reaction_analysis/reaction_frequency_data.npz` - Raw data

## Requirements

### Python Packages

```bash
pip install numpy matplotlib openpmd-api
```

### WarpX Installation

These examples require WarpX with Python/PICMI support. Build with:

```bash
cmake -S . -B build -DWarpX_DIMS=3 -DWarpX_PYTHON=ON
cmake --build build -j
```

For 2D examples:
```bash
cmake -S . -B build_2d -DWarpX_DIMS=2 -DWarpX_PYTHON=ON
cmake --build build_2d -j
```

## Adapting to Your Needs

### Changing Reaction Type

#### Field Ionization (ADK)
```python
ionization = picmi.FieldIonization(
    model="ADK",
    ionized_species=neutral_atoms,
    product_species=[ions, electrons]
)
sim.add_interaction(ionization)
```

#### MCC Collisions
```python
mcc = picmi.MCCCollisions(
    name='mcc_coll',
    species=electrons,
    background_density=1.0e20,
    background_temperature=300.0,
    scattering_processes={
        'ionization': {
            'cross_section': 'path/to/cross_section.dat',
            'energy': 24.6,  # eV
            'species': ions
        }
    }
)
sim.warpx_collisions = [mcc]
```

#### Coulomb Collisions
```python
coulomb = picmi.CoulombCollisions(
    name='coll_ei',
    species=[electrons, ions],
    CoulombLog=15.0
)
sim.warpx_collisions = [coulomb]
```

#### Nuclear Fusion (via input file parameters)
```python
# Note: Nuclear fusion not yet in PICMI, use input file
sim.warpx_set_key('collisions.collision_names', 'fusion')
sim.warpx_set_key('fusion.species', 'deuterium deuterium')
sim.warpx_set_key('fusion.product_species', 'tritium proton')
sim.warpx_set_key('fusion.type', 'nuclearfusion')
```

### Adjusting Spatial Resolution

Change binning in analysis scripts:
```python
NBINS_X = 128  # Increase for finer resolution
NBINS_Y = 128
NBINS_Z = 128
```

### Tracking Different Product Species

In analysis script:
```python
PRODUCT_SPECIES = 'ions'  # or 'electrons', 'alpha', etc.
```

## Understanding the Output

### 2D Spatial Maps

Color intensity represents reaction frequency:
- **Hot spots** (bright) = high reaction rate
- **Cold spots** (dark) = low reaction rate
- Units: weighted particle count per spatial bin

### 1D Profiles

Shows integrated reaction frequency along each axis:
- Sum all reactions perpendicular to the axis
- Useful for identifying peak locations

### Temporal Evolution

Cumulative reaction count vs. simulation time:
- Should generally increase (reactions accumulate)
- Rate of increase shows reaction rate over time

## Tips and Best Practices

### Performance

1. **Start with 2D**: Much faster than 3D for testing
2. **Reduce resolution**: Use coarser grids during development
3. **Fewer timesteps**: Run short simulations first
4. **Diagnostic frequency**: Output every N steps, not every step

### Accuracy

1. **Sufficient statistics**: Need enough macroparticles
   - Typical: 2-8 macroparticles per cell minimum
   - More for better statistics in low-reaction regions

2. **Spatial resolution**: Balance between:
   - Resolution (more bins = finer detail)
   - Statistics (fewer bins = more particles per bin)

3. **Temporal resolution**:
   - Output frequently enough to capture dynamics
   - Not so frequent that you generate excessive data

### Debugging

If you see no reactions:
1. Check product species exist in output
2. Verify reaction conditions are met (field strength, density, etc.)
3. Ensure particle species are configured correctly
4. Check diagnostic output frequency

## Advanced: Real-Time Reaction Rate

For instantaneous (not cumulative) reaction rate:

```python
# In analysis script, track particle creation time
# and bin only newly created particles between outputs

def compute_instantaneous_rate(data_t1, data_t2):
    """Compute reactions between two time steps"""
    # Particles in t2 but not in t1 are newly created
    new_particles = find_new_particles(data_t1, data_t2)
    return bin_by_position(new_particles)
```

## References

### WarpX Documentation
- [Collision Models](https://warpx.readthedocs.io/en/latest/usage/parameters.html#collision-models)
- [PICMI Interface](https://warpx.readthedocs.io/en/latest/usage/picmi.html)
- [Diagnostics](https://warpx.readthedocs.io/en/latest/usage/parameters.html#diagnostics)

### Example Simulations
- `Examples/Tests/collision/` - Coulomb collision tests
- `Examples/Tests/nuclear_fusion/` - Fusion reaction tests
- `Examples/Tests/field_ionization/` - ADK ionization tests
- `Examples/Physics_applications/capacitive_discharge/` - MCC examples

## Contributing

To extend these examples:

1. **New reaction types**: Add example in new file following existing pattern
2. **Better visualization**: Enhance analysis scripts with new plots
3. **Performance**: Optimize binning and I/O operations
4. **Native diagnostic**: Consider implementing C++ diagnostic in WarpX core

## Support

For issues or questions:
1. Check WarpX documentation: https://warpx.readthedocs.io
2. WarpX GitHub discussions: https://github.com/ECP-WarpX/WarpX/discussions
3. WarpX Slack channel

## License

These examples follow WarpX licensing (BSD-3-Clause-LBNL).
