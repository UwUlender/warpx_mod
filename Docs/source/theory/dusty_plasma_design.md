# Dusty Plasma Simulation Design Document

## 1. Physics Background

### 1.1 Dusty Plasma Characteristics

Dusty plasmas (also called complex plasmas) consist of:
- **Electrons**: Light, highly mobile, negative charge
- **Ions**: Heavy, moderately mobile, positive charge
- **Dust particles**: Very heavy (10^6-10^15 amu), radius ~0.1-10 μm, highly charged (10^2-10^5 elementary charges)

Key physical processes:
1. **Dust charging**: Collection of plasma particles on dust surface
2. **Electrostatic interactions**: Via electric field (long-range)
3. **Collisional drag**: Momentum transfer from plasma to dust
4. **Dust-dust interactions**: Through plasma screening

### 1.2 Dust Charging Theory

#### Orbital Motion Limited (OML) Theory

For a spherical dust grain of radius `a` in a collisionless plasma:

**Electron collection current:**
```
I_e = -e * n_e * π * a² * sqrt(8 k_B T_e / (π m_e)) * exp(e φ_d / (k_B T_e))
```

**Ion collection current:**
```
I_i = +e * n_i * π * a² * sqrt(8 k_B T_i / (π m_i)) * (1 - e φ_d / (k_B T_i))
```

where:
- `φ_d` = dust surface potential (negative for typical plasmas)
- `n_e, n_i` = electron and ion densities
- `T_e, T_i` = electron and ion temperatures
- `m_e, m_i` = electron and ion masses
- `a` = dust grain radius

**Equilibrium charge**: At steady-state, `I_e + I_i = 0`, which determines `φ_d` and the dust charge:
```
Q_d = 4π ε_0 a φ_d
```

#### Capacitance Model (Simplified)

For quick estimates:
```
Q_d ≈ -4π ε_0 a T_e * Z_d
```
where `Z_d` is the dimensionless charge number, typically:
```
Z_d ≈ a T_e / (e λ_D)
```
and `λ_D` is the Debye length.

### 1.3 Forces on Dust Particles

1. **Electric field force**: `F_E = Q_d * E`
2. **Ion drag force**: Momentum transfer from streaming ions
   - Collection force: `F_coll ∝ n_i v_i² σ_coll`
   - Coulomb force: `F_Coul ∝ n_i v_i² σ_Coul`
3. **Neutral drag**: `F_n = -m_d ν_dn (v_d - v_n)` where `ν_dn` is collision frequency
4. **Gravity**: `F_g = m_d g`
5. **Thermophoresis**: Force due to temperature gradients

### 1.4 Collision Cross-Sections

**Electron-dust collision:**
- Cross-section: `σ_ed = π a² (1 + 2 e |φ_d| / (m_e v_e²))`
- Includes geometric + Coulomb focusing terms

**Ion-dust collision:**
- Collection: `σ_coll = π a² (1 + 2 e |φ_d| / (m_i v_i²))`
- Coulomb scattering: `σ_Coul ≈ π b_90²` where `b_90 = e Q_d / (4π ε_0 m_i v_i²)`

## 2. Implementation Strategy

### 2.1 Dust Particle Species

We extend WarpX to support dust particles as a special species with:
- Standard PIC attributes: position, momentum, weight
- Additional attributes:
  - `dust_radius`: Physical size of grain
  - `dust_charge`: Current charge (in units of e)
  - `dust_temperature`: Internal temperature (for thermal emission, optional)

### 2.2 Charging Model Implementation

Create a new class `DustChargingFunc` that:
1. Computes local plasma density and temperature from nearby PIC particles
2. Calculates electron and ion collection currents using OML theory
3. Updates dust charge: `dQ/dt = I_e + I_i`
4. Accounts for secondary electron emission (optional, for high energies)

### 2.3 Collision Kernels

Implement collision functors following WarpX patterns:

**`ElectronDustCollisionFunc`:**
- Calculates enhanced cross-section due to attractive potential
- Removes electron (deposited on dust)
- Updates dust charge
- Transfers momentum to dust particle

**`IonDustCollisionFunc`:**
- Handles both collection and Coulomb scattering
- For collection: removes ion, updates charge, transfers momentum
- For scattering: deflects ion trajectory

### 2.4 Field Interaction

Dust particles participate in:
- **Charge deposition**: Deposit charge to grid like other species
- **Field gather**: Experience E and B fields
- **Modified pusher**: Include additional forces (drag, gravity)

### 2.5 Integration into WarpX

1. **Add dust species type**:
   ```cpp
   enum struct PhysicalSpecies {
       ..., dust_micron, dust_submicron
   };
   ```

2. **Extend PhysicalParticleContainer**:
   - Add runtime attributes for dust-specific properties
   - Modify initialization to handle dust parameters

3. **Collision handler**:
   - Register new collision types: `electrondustcollision`, `iondustcollision`
   - Implement charging update in collision loop

4. **Input parameters**:
   ```
   particles.species_names = electrons ions dust

   dust.species_type = dust
   dust.dust_radius = 1e-6  # meters
   dust.dust_initial_charge = -1000  # elementary charges
   dust.mass = 1e-15  # kg

   collisions.coll_names = electron_dust_coll ion_dust_coll
   electron_dust_coll.species = electrons dust
   electron_dust_coll.type = electrondustcollision
   ion_dust_coll.species = ions dust
   ion_dust_coll.type = iondustcollision
   ```

## 3. Physics Validation Tests

### 3.1 Dust Charging Test
- Single dust grain in uniform plasma
- Verify charge reaches equilibrium value predicted by OML theory
- Compare charging time constant

### 3.2 Dust Dynamics Test
- Dust grain in uniform electric field
- Verify acceleration matches `a = Q E / m`
- Test with time-varying charge

### 3.3 Ion Drag Test
- Streaming plasma past stationary dust
- Measure drag force vs. ion velocity
- Compare with theoretical predictions

### 3.4 Dusty Plasma Waves
- Dust acoustic waves (DAW)
- Dust ion-acoustic waves (DIAW)
- Verify dispersion relations

## 4. Computational Considerations

### 4.1 Time Scales

Wide separation of time scales:
- Electron plasma frequency: `ω_pe ~ 10^10 s^-1`
- Ion plasma frequency: `ω_pi ~ 10^8 s^-1`
- Dust plasma frequency: `ω_pd ~ 10^2-10^4 s^-1`
- Dust charging time: `τ_c ~ 10^-6-10^-3 s`

**Strategy**:
- Subcycle dust charging updates
- Use larger time step for dust motion if appropriate

### 4.2 Spatial Scales

- Debye length: `λ_D ~ 10^-5-10^-3 m`
- Dust grain size: `a ~ 10^-7-10^-5 m`
- Grid resolution: Must resolve `λ_D` for accurate field solution

### 4.3 Performance

- Dust particles typically much fewer than electrons/ions (ratio ~10^-6)
- Collision detection: Use spatial binning (already implemented in WarpX)
- Charging calculation: Can be relatively expensive due to local plasma parameter computation

## 5. References

1. Goree, J. (1994). "Charging of particles in a plasma". *Plasma Sources Sci. Technol.* 3, 400.
2. Khrapak, S. A. & Morfill, G. E. (2009). "Basic processes in complex (dusty) plasmas". *Contrib. Plasma Phys.* 49, 148.
3. Fortov, V. E., et al. (2005). "Dusty plasmas". *Phys. Rep.* 421, 1.
4. Shukla, P. K. & Mamun, A. A. (2002). *Introduction to Dusty Plasma Physics*. IOP Publishing.
