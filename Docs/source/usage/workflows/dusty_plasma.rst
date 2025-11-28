.. _usage-workflows-dusty-plasma:

Dusty Plasma Simulations
=========================

WarpX supports the simulation of dusty (complex) plasmas, which consist of electrons, ions, and charged dust particles (typically 0.1-10 μm in size). This capability enables the study of:

- Dust charging dynamics in plasma discharges
- Ion drag forces and dust transport
- Dusty plasma instabilities and waves
- Dust contamination in fusion devices
- Industrial plasma processing with dust

Physics Background
------------------

Dust Charging
^^^^^^^^^^^^^

Dust particles in plasma acquire charge through collection of plasma particles on their surface. The charging process is described by **Orbital Motion Limited (OML) theory** for collisionless plasmas:

**Electron collection current:**

.. math::

   I_e = -e n_e \pi a^2 \sqrt{\frac{8 k_B T_e}{\pi m_e}} \exp\left(\frac{e \phi_d}{k_B T_e}\right)

**Ion collection current:**

.. math::

   I_i = +e n_i \pi a^2 \sqrt{\frac{8 k_B T_i}{\pi m_i}} \left(1 - \frac{e \phi_d}{k_B T_i}\right)

where:

- :math:`a` = dust grain radius
- :math:`\phi_d` = dust surface potential (typically negative)
- :math:`n_e, n_i` = electron and ion densities
- :math:`T_e, T_i` = electron and ion temperatures

At equilibrium, :math:`I_e + I_i = 0`, which determines the dust charge:

.. math::

   Q_d = 4\pi \epsilon_0 a \phi_d

Collision Cross-Sections
^^^^^^^^^^^^^^^^^^^^^^^^^

**Electron-dust**: Enhanced by attractive potential (for negative dust)

.. math::

   \sigma_e = \pi a^2 \left(1 + \frac{2 e |\phi_d|}{m_e v^2}\right)

**Ion-dust**: Reduced by repulsive potential (for negative dust)

.. math::

   \sigma_i = \pi a^2 \max\left(0, 1 - \frac{e \phi_d}{m_i v^2}\right)

Forces on Dust Particles
^^^^^^^^^^^^^^^^^^^^^^^^^

1. **Electric field force**: :math:`\vec{F}_E = Q_d \vec{E}`
2. **Ion drag force**: :math:`\vec{F}_{drag} \propto n_i v_i^2 \sigma_{eff}`
3. **Neutral drag**: :math:`\vec{F}_n = -m_d \nu_{dn} (\vec{v}_d - \vec{v}_n)`
4. **Gravity**: :math:`\vec{F}_g = m_d \vec{g}`

Implementation
--------------

Dust Particle Species
^^^^^^^^^^^^^^^^^^^^^

Dust particles are defined as a special species type with additional attributes:

.. code-block:: bash

   particles.species_names = electrons ions dust

   dust.species_type = dust
   dust.injection_style = "NUniformPerCell"
   dust.mass = 8.38e-15  # kg (1 micron silica sphere)
   dust.charge = 0.0     # initially neutral

   # Add runtime attributes for dust
   dust.addRealAttributes = dust_charge dust_radius
   dust.attribute.dust_charge(x,y,z,ux,uy,uz,t) = 0.0
   dust.attribute.dust_radius(x,y,z,ux,uy,uz,t) = 1.0e-6  # 1 micron

Collision Types
^^^^^^^^^^^^^^^

Two collision types handle dust-plasma interactions:

**Electron-Dust Collisions**

.. code-block:: bash

   collisions.collision_names = electron_dust

   electron_dust.species = electrons dust
   electron_dust.type = electrondustcollision
   electron_dust.dust_radius = 1.0e-6            # meters
   electron_dust.dust_density = 2000.0           # kg/m^3
   electron_dust.collection_efficiency = 1.0
   electron_dust.enable_charging = 1
   electron_dust.enable_momentum_transfer = 1

**Ion-Dust Collisions**

.. code-block:: bash

   ion_dust.species = ions dust
   ion_dust.type = iondustcollision
   ion_dust.dust_radius = 1.0e-6
   ion_dust.collection_efficiency = 0.5
   ion_dust.enable_charging = 1
   ion_dust.enable_momentum_transfer = 1
   ion_dust.enable_coulomb_scattering = 1

Input Parameters
^^^^^^^^^^^^^^^^

Dust Species Parameters
"""""""""""""""""""""""

* ``dust.species_type`` (``string``)
    Set to ``"dust"`` to activate dust particle physics

* ``dust.mass`` (``float``, kg)
    Dust particle mass. For a spherical grain: :math:`m = \frac{4}{3}\pi a^3 \rho_{dust}`

* ``dust.charge`` (``float``, Coulombs)
    Initial charge (typically 0, will evolve dynamically)

* ``dust.addRealAttributes`` (``string list``)
    Runtime attributes: ``dust_charge`` (in units of e) and ``dust_radius`` (in meters). Attribute functions must use signature ``(x,y,z,ux,uy,uz,t)``

Collision Parameters
""""""""""""""""""""

* ``<collision_name>.dust_radius`` (``float``, meters)
    Physical radius of dust grains for cross-section calculation

* ``<collision_name>.dust_density`` (``float``, kg/m³)
    Material density (e.g., 2000 for silica, 2200 for graphite)

* ``<collision_name>.collection_efficiency`` (``float``, 0-1)
    Fraction of collisions resulting in particle collection (default: 1.0 for electrons, 0.5 for ions)

* ``<collision_name>.enable_charging`` (``bool``)
    Enable dynamic charging of dust particles (default: true)

* ``<collision_name>.enable_momentum_transfer`` (``bool``)
    Enable momentum transfer (ion drag) (default: true)

* ``<collision_name>.enable_coulomb_scattering`` (``bool``, ion collisions only)
    Enable Coulomb scattering in addition to collection (default: true)

Example: RF Discharge with Dust
--------------------------------

This example simulates dust particles in a low-pressure RF discharge:

.. code-block:: bash

   # Simulation parameters
   max_step = 10000
   warpx.const_dt = 1.0e-10  # 0.1 ns

   # Grid (2 mm cube)
   amr.n_cell = 64 64 64
   geometry.dims = 3
   geometry.prob_lo = -1.0e-3 -1.0e-3 -1.0e-3
   geometry.prob_hi =  1.0e-3  1.0e-3  1.0e-3

   # Electrostatic solver
   warpx.do_electrostatic = labframe

   # Species
   particles.species_names = electrons ions dust

   electrons.species_type = electron
   electrons.density = 1.0e16  # m^-3
   electrons.theta = 3.0       # 3 eV

   ions.species_type = argon
   ions.density = 1.0e16
   ions.theta = 0.03           # 0.03 eV

   dust.species_type = dust
   dust.density = 1.0e10       # fewer dust particles
   dust.mass = 8.38e-15        # 1 micron silica
   dust.addRealAttributes = dust_charge dust_radius
   dust.attribute.dust_charge(x,y,z,ux,uy,uz,t) = 0.0
   dust.attribute.dust_radius(x,y,z,ux,uy,uz,t) = 1.0e-6

   # Collisions
   collisions.collision_names = e_dust i_dust

   e_dust.species = electrons dust
   e_dust.type = electrondustcollision
   e_dust.dust_radius = 1.0e-6

   i_dust.species = ions dust
   i_dust.type = iondustcollision
   i_dust.dust_radius = 1.0e-6

Analysis and Visualization
---------------------------

Dust Charge Evolution
^^^^^^^^^^^^^^^^^^^^^

Monitor the dust charge using particle diagnostics:

.. code-block:: bash

   diagnostics.diags_names = diag1
   diag1.intervals = 100
   diag1.diag_type = Full
   diag1.dust.variables = w ux uy uz dust_charge dust_radius

Then analyze with Python:

.. code-block:: python

   import openpmd_api as io

   series = io.Series("diags/diag1/openpmd_%T.bp", io.Access.read_only)
   for it in series.iterations:
       dust = series.iterations[it].particles["dust"]
       charges = dust["dust_charge"][io.Mesh_Record_Component.SCALAR][:]
       print(f"Time step {it}: Average charge = {charges.mean():.1f} e")

Comparison with Theory
^^^^^^^^^^^^^^^^^^^^^^

Calculate the OML equilibrium charge for comparison:

.. code-block:: python

   import numpy as np
   from scipy.constants import e, epsilon_0, k, m_e, m_p

   def oml_equilibrium_charge(a, n_e, T_e, n_i, T_i, m_i):
       """Calculate equilibrium dust charge (OML theory)"""
       # Iterative solution for surface potential
       phi = -T_e / 2  # Initial guess
       for _ in range(100):
           I_e = -e * n_e * np.pi * a**2 * np.sqrt(8*k*T_e/(np.pi*m_e)) * np.exp(e*phi/(k*T_e))
           I_i = e * n_i * np.pi * a**2 * np.sqrt(8*k*T_i/(np.pi*m_i)) * (1 - e*phi/(k*T_i))
           if abs(I_e + I_i) < 1e-30:
               break
           dphi = -0.01 * (I_e + I_i) / abs(I_e)
           phi += dphi
       return 4 * np.pi * epsilon_0 * a * phi / e

   # Example: 1 micron dust in argon plasma
   Q_eq = oml_equilibrium_charge(
       a=1e-6, n_e=1e16, T_e=3.0, n_i=1e16, T_i=0.03, m_i=40*m_p
   )
   print(f"OML equilibrium charge: {Q_eq:.1f} e")

Typical Results
---------------

For standard conditions (n_e = 10¹⁶ m⁻³, T_e = 3 eV, a = 1 μm):

- **Equilibrium charge**: -1000 to -5000 elementary charges
- **Surface potential**: -2 to -5 V
- **Charging time**: 1-10 μs
- **Debye length**: λ_D ~ 100 μm

Validation and Accuracy
-----------------------

The implementation has been validated against:

1. OML theory for single dust grain charging
2. Analytic expressions for collection cross-sections
3. Known results for ion drag force

**Valid parameter regime:**

- Dust radius: a << λ_D (OML theory assumption)
- Collisionless plasma: λ_mfp >> λ_D
- Non-relativistic: v_dust << c

**Limitations:**

- Assumes spherical dust grains
- No secondary electron emission
- No thermionic or photoemission
- No dust-dust collisions

Performance Considerations
--------------------------

Timestep Selection
^^^^^^^^^^^^^^^^^^

Multiple timescales must be resolved:

- Plasma: :math:`\Delta t < 1/\omega_{pe}` (typically 10⁻¹¹ s)
- Charging: :math:`\Delta t < \tau_c` (typically 10⁻⁶ s)
- Dust dynamics: :math:`\Delta t < 1/\omega_{pd}` (typically 10⁻⁴ s)

Use the plasma timescale for accurate field resolution.

Computational Cost
^^^^^^^^^^^^^^^^^^

Dust particles are typically much fewer than plasma particles (ratio ~ 10⁻⁶), so the additional computational cost is minimal. The collision calculations add ~10-20% overhead depending on the collision frequency.

References
----------

.. [Goree1994] Goree, J. (1994). "Charging of particles in a plasma". *Plasma Sources Science and Technology*, 3(3), 400.

.. [Khrapak2009] Khrapak, S. A., & Morfill, G. E. (2009). "Basic processes in complex (dusty) plasmas: Charging, interactions, and ion drag force". *Contributions to Plasma Physics*, 49(3), 148-168.

.. [Fortov2005] Fortov, V. E., et al. (2005). "Dusty plasmas". *Physics Reports*, 421(1-2), 1-103.

.. [Shukla2002] Shukla, P. K., & Mamun, A. A. (2002). *Introduction to Dusty Plasma Physics*. Institute of Physics Publishing.
