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
    CollisionHandler* collision_handler = mypc.GetCollisionHandler();
    if (!collision_handler) {
        // No collisions defined, fill with zeros
        mf_dst.setVal(0.0, dcomp, nComp(), 0);
        return;
    }

    // Create a temporary MultiFab with the simulation grid structure at this level
    // We'll accumulate reaction rates here, then interpolate to the diagnostic grid
    const amrex::BoxArray& ba = warpx.boxArray(m_lev);
    const amrex::DistributionMapping& dm = warpx.DistributionMap(m_lev);
    std::unique_ptr<amrex::MultiFab> rates_tmp = std::make_unique<amrex::MultiFab>(ba, dm, 1, 0);
    rates_tmp->setVal(0.0);

    if (!m_collision_name.empty()) {
        // Output rates for a specific collision
        CollisionBase* collision = collision_handler->getCollisionByName(m_collision_name);
        if (!collision) {
            // Collision not found, interpolate zeros
            InterpolateMFForDiag(mf_dst, *rates_tmp, dcomp, dm, m_convertRZmodes2cartesian);
            return;
        }

        amrex::MultiFab* rates_mf = collision->getReactionRatesAtLevel(m_lev);
        if (!rates_mf) {
            // No reaction rate data for this collision at this level
            InterpolateMFForDiag(mf_dst, *rates_tmp, dcomp, dm, m_convertRZmodes2cartesian);
            return;
        }

        // Copy/sum data from the collision's MultiFab to temporary MultiFab
        if (m_reaction_component >= 0) {
            // Copy specific component
            if (m_reaction_component < rates_mf->nComp()) {
                amrex::MultiFab::Copy(*rates_tmp, *rates_mf, m_reaction_component, 0, 1, 0);
            }
        } else {
            // Sum all components
            for (int comp = 0; comp < rates_mf->nComp(); ++comp) {
                amrex::MultiFab::Add(*rates_tmp, *rates_mf, comp, 0, 1, 0);
            }
        }
    } else {
        // Sum rates from all collisions
        const auto& all_collisions = collision_handler->getAllCollisions();

        for (const auto& collision : all_collisions) {
            amrex::MultiFab* rates_mf = collision->getReactionRatesAtLevel(m_lev);
            if (!rates_mf) {
                continue; // Skip collisions without reaction data
            }

            // Sum all components from this collision
            for (int comp = 0; comp < rates_mf->nComp(); ++comp) {
                amrex::MultiFab::Add(*rates_tmp, *rates_mf, comp, 0, 1, 0);
            }
        }
    }

    // Interpolate from simulation grid to diagnostic grid
    InterpolateMFForDiag(mf_dst, *rates_tmp, dcomp, dm, m_convertRZmodes2cartesian);
}
