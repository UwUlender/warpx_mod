/* Copyright 2025
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDNEUTRALDRAG_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDNEUTRALDRAG_H_

#include "Particles/MultiParticleContainer.H"
#include "Particles/Collision/CollisionBase.H"

#include <AMReX_REAL.H>
#include <AMReX_Parser.H>

#include <string>

/**
 * \brief Class for Neutral Gas Drag (Epstein Drag) on particles (e.g. Dust).
 *
 * Applies a drag force F = - delta * v
 * where delta = (8/3) * sqrt(2*pi) * a^2 * m_n * n_n * v_th_n
 * a: particle radius
 * m_n: neutral mass
 * n_n: neutral density (spatially varying)
 * v_th_n: neutral thermal velocity (sqrt(kB * T_n / m_n))
 */
class BackgroundNeutralDrag final
    : public CollisionBase
{
public:
    BackgroundNeutralDrag (const std::string& collision_name);

    ~BackgroundNeutralDrag () override = default;

    BackgroundNeutralDrag ( BackgroundNeutralDrag const &)             = delete;
    BackgroundNeutralDrag& operator= ( BackgroundNeutralDrag const & ) = delete;
    BackgroundNeutralDrag ( BackgroundNeutralDrag&& )                  = delete;
    BackgroundNeutralDrag& operator= ( BackgroundNeutralDrag&& )       = delete;

    /** Perform the drag calculation
     *
     * @param cur_time Current time
     * @param dt Time step size
     * @param mypc Container of species involved
     *
     */
    void doCollisions (amrex::Real cur_time, amrex::Real dt, MultiParticleContainer* mypc) override;

    /** Perform the drag calculation within a tile
     *
     * @param pti particle iterator
     * @param dt time step size
     * @param t current simulation time
     * @param species_mass mass of the active species
     */
    void doNeutralDragWithinTile (WarpXParIter& pti, amrex::Real dt, amrex::Real t,
                                  amrex::ParticleReal species_mass);

private:

    amrex::ParticleReal m_neutral_mass;
    amrex::ParticleReal m_dust_radius;

    amrex::Parser m_neutral_density_parser;
    amrex::Parser m_neutral_temperature_parser;

    amrex::ParserExecutor<4> m_neutral_density_func;
    amrex::ParserExecutor<4> m_neutral_temperature_func;

};

#endif // WARPX_PARTICLES_COLLISION_BACKGROUNDNEUTRALDRAG_H_
