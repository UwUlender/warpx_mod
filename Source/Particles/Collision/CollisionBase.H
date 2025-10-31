/* Copyright 2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_COLLISIONBASE_H_
#define WARPX_PARTICLES_COLLISION_COLLISIONBASE_H_

#include "Particles/MultiParticleContainer_fwd.H"

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <string>

class CollisionBase
{
public:

    CollisionBase (const std::string& collision_name);

    virtual void doCollisions (amrex::Real /*cur_time*/, amrex::Real /*dt*/, MultiParticleContainer* /*mypc*/ ){}

    CollisionBase(CollisionBase const &) = delete;
    CollisionBase(CollisionBase &&) = delete;
    CollisionBase & operator=(CollisionBase const &) = delete;
    CollisionBase & operator=(CollisionBase const &&) = delete;

    virtual ~CollisionBase() = default;

    [[nodiscard]] int get_ndt() const {return m_ndt;}

    [[nodiscard]] bool use_global_debye_length() const {return m_use_global_debye_length;}

    /**
     * \brief Get reaction counts MultiFab for a specific level (virtual, overridden by BinaryCollision)
     * @param[in] lev The refinement level
     * @return Pointer to MultiFab, or nullptr if not available
     */
    [[nodiscard]] virtual amrex::MultiFab* getReactionCountsAtLevel([[maybe_unused]] int lev) const { return nullptr; }

    /**
     * \brief Get reaction rates MultiFab for a specific level (virtual, overridden by BinaryCollision)
     * @param[in] lev The refinement level
     * @return Pointer to MultiFab, or nullptr if not available
     */
    [[nodiscard]] virtual amrex::MultiFab* getReactionRatesAtLevel([[maybe_unused]] int lev) const { return nullptr; }

    /**
     * \brief Get reaction names (virtual, overridden by BinaryCollision)
     * @return Vector of reaction names
     */
    [[nodiscard]] virtual const amrex::Vector<std::string>& getReactionNames() const {
        static const amrex::Vector<std::string> empty;
        return empty;
    }

protected:

    amrex::Vector<std::string> m_species_names;
    int m_ndt;

    bool m_use_global_debye_length = false;

};

#endif // WARPX_PARTICLES_COLLISION_COLLISIONBASE_H_
