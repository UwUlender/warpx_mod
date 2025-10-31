/* Copyright 2025 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ReactionRateFunctor.H"

#include "WarpX.H"

#include <AMReX.H>
#include <AMReX_MultiFab.H>

ReactionRateFunctor::ReactionRateFunctor (int lev,
                                          amrex::IntVect crse_ratio,
                                          const std::string& collision_name,
                                          int reaction_component,
                                          bool convertRZmodes2cartesian,
                                          int ncomp)
    : ComputeDiagFunctor(ncomp, crse_ratio),
      m_lev(lev),
      m_collision_name(collision_name),
      m_reaction_component(reaction_component),
      m_convertRZmodes2cartesian(convertRZmodes2cartesian)
{
}

void
ReactionRateFunctor::operator() (amrex::MultiFab& mf_dst, int dcomp, int /*i_buffer*/) const
{
    auto& warpx = WarpX::GetInstance();
    auto& mypc = warpx.GetPartContainer();

    // Get the collision handler
    if (!mypc.collisionhandler) {
        // No collisions defined, fill with zeros
        mf_dst.setVal(0.0, dcomp, nComp(), 0);
        return;
    }

    // Initialize output to zero
    mf_dst.setVal(0.0, dcomp, nComp(), 0);

    if (!m_collision_name.empty()) {
        // Output rates for a specific collision
        CollisionBase* collision = mypc.collisionhandler->getCollisionByName(m_collision_name);
        if (!collision) {
            // Collision not found, output remains zero
            return;
        }

        amrex::MultiFab* rates_mf = collision->getReactionRatesAtLevel(m_lev);
        if (!rates_mf) {
            // No reaction rate data for this collision at this level
            return;
        }

        // Copy data from the collision's MultiFab to output
        if (m_reaction_component >= 0) {
            // Copy specific component
            if (m_reaction_component < rates_mf->nComp()) {
                amrex::MultiFab::Copy(mf_dst, *rates_mf, m_reaction_component, dcomp, 1, 0);
            }
        } else {
            // Sum all components
            for (int comp = 0; comp < rates_mf->nComp(); ++comp) {
                amrex::MultiFab::Add(mf_dst, *rates_mf, comp, dcomp, 1, 0);
            }
        }
    } else {
        // Sum rates from all collisions
        const auto& all_collisions = mypc.collisionhandler->getAllCollisions();

        for (const auto& collision : all_collisions) {
            amrex::MultiFab* rates_mf = collision->getReactionRatesAtLevel(m_lev);
            if (!rates_mf) {
                continue; // Skip collisions without reaction data
            }

            // Sum all components from this collision
            for (int comp = 0; comp < rates_mf->nComp(); ++comp) {
                amrex::MultiFab::Add(mf_dst, *rates_mf, comp, dcomp, 1, 0);
            }
        }
    }

    // Apply coarsening if needed (for compatibility with other functors)
    // Note: For now, we assume reaction rate data is already on the diagnostic level
    // If interpolation is needed, it can be added here similar to RhoFunctor
}
