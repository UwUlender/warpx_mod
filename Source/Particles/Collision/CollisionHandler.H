/* Copyright 2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_COLLISIONHANDLER_H_
#define WARPX_PARTICLES_COLLISION_COLLISIONHANDLER_H_

#include "CollisionBase.H"

#include "Particles/MultiParticleContainer_fwd.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <memory>
#include <string>

/* \brief CollisionHandler is a light weight class that contains the
 * list of collisions to be done.
 */
class CollisionHandler
{
public:
    CollisionHandler (const MultiParticleContainer*  mypc);

    /* Perform all of the collisions */
    void doCollisions (int step, amrex::Real cur_time, amrex::Real dt, MultiParticleContainer* mypc);

    /**
     * \brief Get collision by name
     * @param[in] name Name of the collision
     * @return Pointer to CollisionBase, or nullptr if not found
     */
    [[nodiscard]] CollisionBase* getCollisionByName(const std::string& name) const;

    /**
     * \brief Get all collision names
     * @return Vector of collision names
     */
    [[nodiscard]] const amrex::Vector<std::string>& getCollisionNames() const { return collision_names; }

    /**
     * \brief Get all collisions
     * @return Vector of collision pointers
     */
    [[nodiscard]] const amrex::Vector<std::unique_ptr<CollisionBase>>& getAllCollisions() const { return allcollisions; }

private:

    amrex::Vector<std::string> collision_names;
    amrex::Vector<std::string> collision_types;
    amrex::Vector< std::unique_ptr<CollisionBase> > allcollisions;

    bool m_use_global_debye_length = false;

};

#endif // WARPX_PARTICLES_COLLISION_COLLISIONHANDLER_H_
