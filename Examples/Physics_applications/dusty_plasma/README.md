# Dusty Plasma Simulation with WarpX

This example demonstrates the simulation of dusty plasmas using WarpX's new dust particle interaction capabilities.

## Overview

Dusty plasmas (also called complex plasmas) consist of electrons, ions, and micron-sized dust particles. The dust particles acquire charge through collection of plasma particles and interact with the plasma through:

1. **Dynamic charging**: Electrons and ions collected on dust surface (OML theory)
2. **Electrostatic interactions**: Charged dust particles in plasma electric fields
3. **Momentum transfer**: Ion drag force from plasma-dust collisions
4. **Coulomb scattering**: Long-range interactions

## Physics Implementation

### Dust Particle Species

Dust particles are implemented as a special particle species with:
- Standard PIC attributes: position, momentum, weight
- Additional runtime attributes:
  - `dust_charge`: Current charge in units of elementary charge
  - `dust_radius`: Physical size of grain (meters)

### Charging Model

The implementation uses **Orbital Motion Limited (OML) theory**:

**Electron collection**:
- Enhanced cross-section due to attractive potential (for negative dust)
- σ_e = π a² (1 + 2 e |φ_d| / (m_e v²))

**Ion collection**:
- Reduced cross-section due to repulsive potential (for negative dust)
- σ_i = π a² (1 - e φ_d / (m_i v²))

**Equilibrium**: At steady-state, electron and ion currents balance, determining the dust charge.

### Collision Types

Two new collision types are available:

1. **`electrondustcollision`**: Electron-dust interactions
   - Electron collection and charging
   - Momentum transfer to dust
   - Enhanced cross-section for attractive potential

2. **`iondustcollision`**: Ion-dust interactions
   - Ion collection and charging
   - Ion drag force (momentum transfer)
   - Coulomb scattering
   - Reduced collection due to potential barrier

## Running the Example

### Basic Usage

```bash
mpirun -np 4 warpx.3d inputs_dusty_plasma_3d
```

### Input File Structure

The input file `inputs_dusty_plasma_3d` contains:

1. **Species definitions**:
   ```
   particles.species_names = electrons ions dust

   dust.species_type = dust
   dust.addRealAttributes = dust_charge dust_radius
   dust.attribute.dust_charge(x,y,z,ux,uy,uz,t) = 0.0
   dust.attribute.dust_radius(x,y,z,ux,uy,uz,t) = 1.0e-6
   ```

2. **Collision setup**:
   ```
   collisions.collision_names = electron_dust ion_dust

   electron_dust.species = electrons dust
   electron_dust.type = electrondustcollision
   electron_dust.dust_radius = 1.0e-6
   ```

### Important Parameters

**Dust properties**:
- `dust.mass`: Dust particle mass (kg)
- `dust.charge`: Initial charge (usually 0)
- `dust_radius`: Physical radius (m)
- `dust_density`: Material density (kg/m³)

**Collision parameters**:
- `dust_radius`: Grain radius for cross-section calculation
- `collection_efficiency`: 0-1, fraction of collisions that result in collection
- `enable_charging`: Enable/disable dynamic charging
- `enable_momentum_transfer`: Enable/disable ion drag
- `enable_coulomb_scattering`: Enable/disable Coulomb scattering (ions only)

## Analysis

### Using the Analysis Script

```bash
python analysis_dusty_plasma.py ./diags/diag1/
```

The script produces:
1. **Dust charging evolution**: Time evolution of dust charge vs. OML prediction
2. **Velocity distribution**: Dust particle velocity distribution
3. **Summary statistics**: Comparison with theory

### Expected Results

For typical parameters (n_e = 10¹⁶ m⁻³, T_e = 3 eV, a = 1 μm):
- Equilibrium dust charge: ~ -1000 to -5000 e
- Charging time: ~ 1-10 μs
- Surface potential: ~ -2 to -5 V

## Physical Parameters

### Typical Dusty Plasma Regimes

| Regime | n_e (m⁻³) | T_e (eV) | a (μm) | Q_dust (e) |
|--------|-----------|----------|--------|------------|
| RF discharge | 10¹⁵-10¹⁷ | 1-5 | 0.1-10 | -100 to -10⁴ |
| DC discharge | 10¹⁴-10¹⁶ | 1-3 | 1-10 | -500 to -5×10⁴ |
| Fusion edge | 10¹⁸-10¹⁹ | 10-100 | 0.01-1 | -10³ to -10⁵ |

### Time Scales

- Electron plasma period: 1/ω_pe ~ 10⁻¹¹ s
- Ion plasma period: 1/ω_pi ~ 10⁻⁹ s
- Dust charging time: τ_c ~ 10⁻⁶ to 10⁻³ s
- Dust plasma period: 1/ω_pd ~ 10⁻⁴ to 10⁻² s

## Validation

### Test Cases

1. **Single dust charging**: Compare equilibrium charge with OML theory
2. **Dust in electric field**: Verify F = Q E dynamics
3. **Ion drag**: Measure drag force vs. streaming velocity
4. **Dust acoustic waves**: Dispersion relation verification

### Known Limitations

1. Current implementation assumes spherical dust grains
2. OML theory valid for: a << λ_D (Debye length)
3. Collisionless plasma assumption (mean free path >> λ_D)
4. No dust-dust collisions (hard sphere)
5. No secondary electron emission (for high-energy impacts)

## References

1. Goree, J. (1994). "Charging of particles in a plasma". *Plasma Sources Sci. Technol.* 3, 400.
2. Khrapak, S. A. & Morfill, G. E. (2009). "Basic processes in complex (dusty) plasmas". *Contrib. Plasma Phys.* 49, 148.
3. Fortov, V. E., et al. (2005). "Dusty plasmas". *Phys. Rep.* 421, 1.
4. Shukla, P. K. & Mamun, A. A. (2002). *Introduction to Dusty Plasma Physics*. IOP Publishing.

## Future Enhancements

Potential extensions to this implementation:

1. **Capacitance model**: Alternative simplified charging model
2. **Secondary emission**: Include electron emission from dust surface
3. **Thermionic emission**: For high-temperature dust
4. **Photoemission**: For UV-illuminated dust
5. **Dust-dust interactions**: Hard sphere collisions
6. **Non-spherical grains**: Ellipsoids, irregular shapes
7. **Grain size distribution**: Polydisperse dust populations
8. **Dust growth/ablation**: Dynamic grain size evolution

## Contact

For questions or issues related to dusty plasma simulations in WarpX, please:
- Open an issue on the WarpX GitHub repository
- Consult the WarpX documentation at https://warpx.readthedocs.io/

## License

This example is part of WarpX and is distributed under the BSD-3-Clause-LBNL license.
