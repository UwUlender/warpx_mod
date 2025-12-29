/* Copyright 2025
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "BackgroundNeutralDrag.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/WarpXConst.H"

#include <AMReX_ParmParse.H>

using namespace amrex;

BackgroundNeutralDrag::BackgroundNeutralDrag (const std::string& collision_name)
    : CollisionBase(collision_name)
{
    const ParmParse pp_collision_name(collision_name);

    // Read neutral gas properties
    utils::parser::getWithParser(
        pp_collision_name, "neutral_mass", m_neutral_mass);
    
    // Read dust radius (parameter for now, could be attribute)
    utils::parser::getWithParser(
        pp_collision_name, "dust_radius", m_dust_radius);

    // Read density and temperature parsers
    std::string neutral_density_string;
    utils::parser::Store_parserString(
        pp_collision_name, "neutral_density", neutral_density_string);
    m_neutral_density_parser = utils::parser::makeParser(
        neutral_density_string, {"x", "y", "z", "t"});
    m_neutral_density_func = m_neutral_density_parser.compile<4>();

    std::string neutral_temperature_string;
    utils::parser::Store_parserString(
        pp_collision_name, "neutral_temperature", neutral_temperature_string);
    m_neutral_temperature_parser = utils::parser::makeParser(
        neutral_temperature_string, {"x", "y", "z", "t"});
    m_neutral_temperature_func = m_neutral_temperature_parser.compile<4>();
}

void BackgroundNeutralDrag::doCollisions (amrex::Real cur_time, amrex::Real dt, MultiParticleContainer* mypc)
{
    const int n_species = m_species_names.size();
    for (int i = 0; i < n_species; ++i)
    {
        const std::string& species_name = m_species_names[i];
        WarpXParticleContainer& species = mypc->GetParticleContainerFromName(species_name);
        amrex::ParticleReal mass = species.getMass();

        for (WarpXParIter pti(species, 0); pti.isValid(); ++pti)
        {
            doNeutralDragWithinTile(pti, dt, cur_time, mass);
        }
    }
}

void BackgroundNeutralDrag::doNeutralDragWithinTile (WarpXParIter& pti, amrex::Real dt, amrex::Real t,
                                                     amrex::ParticleReal species_mass)
{
    using namespace amrex::literals;

    const long np = pti.numParticles();
    
    auto& attribs = pti.GetAttribs();
    amrex::ParticleReal* AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
    amrex::ParticleReal* AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
    amrex::ParticleReal* AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();

    const auto GetPosition = GetParticlePosition<PIdx>(pti);

    // Capture members
    amrex::ParticleReal neutral_mass = m_neutral_mass;
    amrex::ParticleReal dust_radius = m_dust_radius;
    auto neutral_density_func = m_neutral_density_func;
    auto neutral_temperature_func = m_neutral_temperature_func;

    // Constants
    constexpr amrex::ParticleReal pi = MathConst::pi;
    constexpr amrex::ParticleReal kb = PhysConst::kb;

    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (long ip)
    {
        amrex::ParticleReal xp, yp, zp;
        GetPosition(ip, xp, yp, zp);

        // Evaluate neutral properties at particle position
        amrex::ParticleReal n_n = neutral_density_func(xp, yp, zp, t);
        amrex::ParticleReal T_n = neutral_temperature_func(xp, yp, zp, t);

        if (n_n <= 0.0 || T_n <= 0.0) return;

        // Calculate thermal velocity of neutrals
        // v_th = sqrt(kB * T / m)
        amrex::ParticleReal v_th_n = std::sqrt(kb * T_n / neutral_mass);

        // Calculate Epstein drag coefficient delta
        // delta = (8/3) * sqrt(2*pi) * a^2 * m_n * n_n * v_th_n
        amrex::ParticleReal delta = (8.0_prt / 3.0_prt) * std::sqrt(2.0_prt * pi) *
                                    dust_radius * dust_radius * neutral_mass * n_n * v_th_n;

        // Apply drag force
        // F = -delta * v
        // du/dt = F/m = -(delta/m) * v
        // Implicit update: u_new = u_old / (1 + (delta/m)*dt)
        // Explicit update: u_new = u_old * (1 - (delta/m)*dt)
        // Using implicit for stability if drag is strong
        
        // Note: ux, uy, uz are gamma*v. For non-relativistic dust, gamma ~ 1.
        // F = -delta * (u/gamma)
        // du/dt = F/m = -(delta/m) * u/gamma
        // We assume gamma=1 for dust (heavy, slow).
        
        amrex::ParticleReal damping_factor = 1.0_prt / (1.0_prt + (delta / species_mass) * dt);
        
        ux[ip] *= damping_factor;
        uy[ip] *= damping_factor;
        uz[ip] *= damping_factor;
    });
}
