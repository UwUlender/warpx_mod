/* Copyright 2025 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_REACTIONRATEFUNCTOR_H_
#define WARPX_REACTIONRATEFUNCTOR_H_

#include "ComputeDiagFunctor.H"

#include <AMReX_BaseFwd.H>

#include <string>

/**
 * \brief Functor to compute reaction rates (reactions per unit volume) into mf_out
 *
 * This functor outputs the reaction rate data stored in binary collision MultiFabs.
 * It can output either total reaction rates (summed across all collision types)
 * or rates for a specific collision.
 */
class ReactionRateFunctor final : public ComputeDiagFunctor
{

public:

    /**
     * \brief Constructor
     *
     * \param[in] lev level of MultiFab
     * \param[in] crse_ratio coarsening ratio for interpolation of field values
     *                       from simulation MultiFabs to the output MultiFab mf_dst
     * \param[in] collision_name Name of collision to dump rates for (optional).
     *                          If empty, sums rates from all collisions.
     * \param[in] reaction_component Component index of the reaction subtype (optional).
     *                              -1 means sum all reaction subtypes for this collision.
     * \param[in] convertRZmodes2cartesian if true, all RZ modes are averaged into one component
     * \param[in] ncomp optional number of component of source MultiFab mf_src
     *                  to be cell-centered in output MultiFab mf_dst
     */
    ReactionRateFunctor (int lev,
                amrex::IntVect crse_ratio,
                const std::string& collision_name = "",
                int reaction_component = -1,
                bool convertRZmodes2cartesian = true,
                int ncomp = 1);

    /**
     * \brief Compute reaction rate directly into mf_dst
     *
     * \param[out] mf_dst output MultiFab where the result is written
     * \param[in]  dcomp  first component of mf_dst in which cell-centered data are stored
     */
    void operator() ( amrex::MultiFab& mf_dst, int dcomp, int /*i_buffer=0*/ ) const override;

private:

    // Level on which source MultiFab is defined
    int m_lev;

    // Name of collision to output rates for (empty = sum all collisions)
    std::string m_collision_name;

    // Reaction component index (-1 = sum all components for this collision)
    int m_reaction_component;

    // Whether to average all modes into one component in RZ geometry
    bool m_convertRZmodes2cartesian;
};

#endif // WARPX_REACTIONRATEFUNCTOR_H_
