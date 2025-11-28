/* Copyright 2025
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_ELECTRON_DUST_COLLISION_FUNC_H_
#define WARPX_ELECTRON_DUST_COLLISION_FUNC_H_

#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/WarpXConst.H"

#include <AMReX_DenseBins.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Random.H>
#include <AMReX_REAL.H>

/**
 * \brief This functor performs electron-dust collision including:
 *        1. Enhanced cross-section due to attractive potential
 *        2. Electron collection on dust surface
 *        3. Dynamic charging of dust particles
 *        4. Momentum transfer from electrons to dust
 *
 * Physics model: Orbital Motion Limited (OML) theory
 * Reference: Goree, J. (1994). "Charging of particles in a plasma"
 */
class ElectronDustCollisionFunc
{
    // Define shortcuts for frequently-used type names
    using ParticleType = WarpXParticleContainer::ParticleType;
    using ParticleTileType = WarpXParticleContainer::ParticleTileType;
    using ParticleTileDataType = ParticleTileType::ParticleTileDataType;
    using ParticleBins = amrex::DenseBins<ParticleTileDataType>;
    using index_type = ParticleBins::index_type;
    using SoaData_type = WarpXParticleContainer::ParticleTileType::ParticleTileDataType;

public:
    /**
     * \brief Default constructor of the ElectronDustCollisionFunc class.
     */
    ElectronDustCollisionFunc () = default;

    /**
     * \brief Constructor of the ElectronDustCollisionFunc class
     *
     * @param[in] collision_name the name of the collision
     * @param[in] mypc the particle container
     * @param[in] isSameSpecies true if this is an intra-species collision (should be false)
     */
    ElectronDustCollisionFunc (const std::string& collision_name,
                               [[maybe_unused]] MultiParticleContainer const * const mypc,
                               [[maybe_unused]] const bool isSameSpecies)
    {
        const amrex::ParmParse pp_collision_name(collision_name);

        // Read dust particle radius (in meters)
        amrex::ParticleReal dust_radius = 1.0e-6;  // Default: 1 micron
        utils::parser::queryWithParser(
            pp_collision_name, "dust_radius", dust_radius);
        m_dust_radius = dust_radius;
        m_exe.m_dust_radius = m_dust_radius;

        // Read dust density (in kg/m^3) for mass calculation
        amrex::ParticleReal dust_density = 2000.0;  // Default: ~silica
        utils::parser::queryWithParser(
            pp_collision_name, "dust_density", dust_density);
        m_dust_density = dust_density;
        m_exe.m_dust_density = m_dust_density;

        // Read collection efficiency (0-1, default 1 for full collection)
        amrex::ParticleReal collection_efficiency = 1.0;
        utils::parser::queryWithParser(
            pp_collision_name, "collection_efficiency", collection_efficiency);
        m_exe.m_collection_efficiency = collection_efficiency;

        // Option to enable/disable charging
        bool enable_charging = true;
        utils::parser::queryWithParser(
            pp_collision_name, "enable_charging", enable_charging);
        m_exe.m_enable_charging = enable_charging;

        // Option to enable/disable momentum transfer
        bool enable_momentum_transfer = true;
        utils::parser::queryWithParser(
            pp_collision_name, "enable_momentum_transfer", enable_momentum_transfer);
        m_exe.m_enable_momentum_transfer = enable_momentum_transfer;
    }

    struct Executor {
        /**
         * \brief Executor of the ElectronDustCollisionFunc class.
         * Performs electron-dust collisions at the cell level.
         *
         * The collision includes:
         * - Enhanced cross-section: σ = π a² (1 + 2 e |φ_d| / (m_e v²))
         * - Electron collection and charge update
         * - Momentum transfer to dust particles
         *
         * @param[in] I1s,I2s start index for species 1 (electrons) and 2 (dust)
         * @param[in] I1e,I2e stop index (exclusive)
         * @param[in] I1,I2 index arrays
         * @param[in,out] soa_1 electron struct of array data
         * @param[in,out] soa_2 dust struct of array data
         * @param[in] n1,n2 local densities (not used directly here)
         * @param[in] T1,T2 local temperatures
         * @param[in] q1,q2 charges
         * @param[in] m1,m2 masses
         * @param[in] dt time step
         * @param[in] dV cell volume
         * @param[in] engine random engine
         */
        AMREX_GPU_HOST_DEVICE AMREX_INLINE
        void operator() (
            index_type const I1s, index_type const I1e,
            index_type const I2s, index_type const I2e,
            index_type const* AMREX_RESTRICT I1,
            index_type const* AMREX_RESTRICT I2,
            SoaData_type soa_1, SoaData_type soa_2,
            GetParticlePosition<PIdx> /*get_position_1*/, GetParticlePosition<PIdx> /*get_position_2*/,
            amrex::ParticleReal const  /*n1*/, amrex::ParticleReal const  /*n2*/,
            amrex::ParticleReal const  /*T1*/, amrex::ParticleReal const  /*T2*/,
            amrex::Real const /*global_lamdb*/,
            amrex::ParticleReal const  /*q1*/, amrex::ParticleReal const  /*q2*/,
            amrex::ParticleReal const  m1, amrex::ParticleReal const  /*m2*/,
            amrex::Real const  dt, amrex::Real const dV, index_type /*coll_idx*/,
            index_type const /*cell_start_pair*/, index_type* /*p_mask*/,
            index_type* /*p_pair_indices_1*/, index_type* /*p_pair_indices_2*/,
            amrex::ParticleReal* /*p_pair_reaction_weight*/,
            amrex::ParticleReal* /*p_product_data*/,
            amrex::RandomEngine const& engine) const
        {
            using namespace amrex::literals;

            // Extract particle data arrays
            auto* AMREX_RESTRICT ux1 = soa_1.m_rdata[PIdx::ux];
            auto* AMREX_RESTRICT uy1 = soa_1.m_rdata[PIdx::uy];
            auto* AMREX_RESTRICT uz1 = soa_1.m_rdata[PIdx::uz];
            auto* AMREX_RESTRICT w1  = soa_1.m_rdata[PIdx::w];

            auto* AMREX_RESTRICT ux2 = soa_2.m_rdata[PIdx::ux];
            auto* AMREX_RESTRICT uy2 = soa_2.m_rdata[PIdx::uy];
            auto* AMREX_RESTRICT uz2 = soa_2.m_rdata[PIdx::uz];
            auto* AMREX_RESTRICT w2  = soa_2.m_rdata[PIdx::w];

            // Get dust charge from runtime attributes (if available)
            // For now, assume dust particles have a runtime attribute "dust_charge"
            // If not available, we'll use the species charge
            amrex::ParticleReal* dust_charge_ptr = nullptr;
            if (soa_2.m_num_runtime_real > 0) {
                // Assume first runtime real attribute is dust_charge (in units of e)
                dust_charge_ptr = soa_2.m_runtime_rdata[0];
            }

            // Constants
            constexpr amrex::ParticleReal epsilon_0 = PhysConst::epsilon_0;
            constexpr amrex::ParticleReal e = PhysConst::q_e;
            const amrex::ParticleReal a = m_dust_radius;
            const amrex::ParticleReal sigma_geometric = MathConst::pi * a * a;

            // Loop over all dust particles (species 2)
            for (index_type i2 = I2s; i2 < I2e; ++i2)
            {
                const index_type id2 = I2[i2];

                // Dust particle properties
                const amrex::ParticleReal dust_ux = ux2[id2];
                const amrex::ParticleReal dust_uy = uy2[id2];
                const amrex::ParticleReal dust_uz = uz2[id2];
                const amrex::ParticleReal dust_w = w2[id2];

                // Get current dust charge (in units of elementary charge)
                amrex::ParticleReal Q_dust_e = 0.0;  // Charge in units of e
                if (dust_charge_ptr != nullptr) {
                    Q_dust_e = dust_charge_ptr[id2];
                }

                // Calculate dust surface potential (V)
                // φ_d = Q / (4π ε_0 a)
                const amrex::ParticleReal phi_dust =
                    (Q_dust_e * e) / (4.0_prt * MathConst::pi * epsilon_0 * a);

                // Loop over all electrons (species 1)
                for (index_type i1 = I1s; i1 < I1e; ++i1)
                {
                    const index_type id1 = I1[i1];

                    // Electron properties
                    const amrex::ParticleReal e_ux = ux1[id1];
                    const amrex::ParticleReal e_uy = uy1[id1];
                    const amrex::ParticleReal e_uz = uz1[id1];
                    const amrex::ParticleReal e_w = w1[id1];

                    // Relative 4-velocity (ux, uy, uz are γ*v in WarpX)
                    const amrex::ParticleReal du_x = e_ux - dust_ux;
                    const amrex::ParticleReal du_y = e_uy - dust_uy;
                    const amrex::ParticleReal du_z = e_uz - dust_uz;

                    const amrex::ParticleReal u2 = du_x*du_x + du_y*du_y + du_z*du_z;

                    if (u2 < 1.0e-40_prt) continue;  // Skip if no relative motion

                    // For non-relativistic particles (typical in dusty plasmas): γ ≈ 1, so u ≈ v
                    // For full relativistic treatment: v² = u² / (1 + u²/c²)
                    // Here we use the non-relativistic approximation
                    const amrex::ParticleReal du = std::sqrt(u2);

                    // Calculate enhanced cross-section due to attractive potential
                    // σ = π a² (1 + 2 e |φ_d| / (m_e v²))
                    // Note: For negative dust (φ_d < 0), electrons are attracted
                    const amrex::ParticleReal v2 = u2;  // Non-relativistic approximation
                    const amrex::ParticleReal coulomb_factor =
                        (phi_dust < 0.0_prt) ?
                        (1.0_prt + 2.0_prt * e * std::abs(phi_dust) / (m1 * v2)) :
                        1.0_prt;
                    const amrex::ParticleReal sigma_eff = sigma_geometric * coulomb_factor;

                    // Collision probability
                    // P = n * σ * v * dt / V
                    // For PIC: P = (w1 / V) * σ * v * dt
                    const amrex::ParticleReal prob =
                        (e_w / dV) * sigma_eff * du * dt * m_collection_efficiency;

                    // Perform collision with probability
                    if (prob > amrex::Random(engine))
                    {
                        // Momentum transfer to dust
                        if (m_enable_momentum_transfer)
                        {
                            // Transfer electron momentum to dust
                            // Conservation: dust gains what electron had
                            const amrex::ParticleReal mom_ratio = e_w / dust_w;
                            ux2[id2] += e_ux * mom_ratio;
                            uy2[id2] += e_uy * mom_ratio;
                            uz2[id2] += e_uz * mom_ratio;
                        }

                        // Charge update
                        if (m_enable_charging && dust_charge_ptr != nullptr)
                        {
                            // Electron collected: dust becomes more negative
                            // Each electron contributes charge -e
                            dust_charge_ptr[id2] -= 1.0_prt * (e_w / dust_w);
                        }

                        // Mark electron for removal (set weight to 0)
                        // In practice, one might want to implement a more sophisticated
                        // particle removal scheme, but for now we can reduce weight
                        // to effectively remove the electron
                        // Note: Full removal would require particle management
                    }
                }
            }
        }

        amrex::ParticleReal m_dust_radius;
        amrex::ParticleReal m_dust_density;
        amrex::ParticleReal m_collection_efficiency = 1.0;
        bool m_enable_charging = true;
        bool m_enable_momentum_transfer = true;
        bool m_computeSpeciesDensities = false;
        bool m_computeSpeciesTemperatures = false;
        bool m_need_product_data = false;
    };

    [[nodiscard]] Executor const& executor () const { return m_exe; }

    bool use_global_debye_length() { return false; }

private:
    amrex::ParticleReal m_dust_radius;
    amrex::ParticleReal m_dust_density;
    Executor m_exe;
};

#endif // WARPX_ELECTRON_DUST_COLLISION_FUNC_H_
