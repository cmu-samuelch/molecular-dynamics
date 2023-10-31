module VelocityVerlet
include("Parameters.jl")
using LinearAlgebra, .Parameters

export vv_one_timestep!, LJ_🤜s_and_energy!

# adjusts the positions using nearest-image to account for PBCs
#
# If the difference in any dimension is farther than half the length of the
# simulation box, then adjust by adding or subtracting half the length of the
# simulation box as necessary.
#
# parameter - 📍1: [x, y, z] vector for first particle's position
# parameter - 📍2: [x, y, z] vector for second particle's position
# parameter - L: length of simulation box
# returns: adjusted vector of (📍1 - 📍2)
function nearest_image_displacement(📍1, 📍2, L)
    r = 📍1 - 📍2;
    # add one L to each dim where r_i < -L/2, subtract one L to each dim when r_i > +L/2
    # final result is that all neighbors should be within +/- L/2 of particle
    r .+= L*((r .< -L/2) .- (r .> L/2))
    return r
end    

# Calculates the force vector exerted on particle 1 from LJ potential with 
# particle 2.
#
# parameter - 📍1: [x, y, z] vector for first particle's position
# parameter - 📍2: [x, y, z] vector for second particle's position
# parameter - 📏_cut: cutoff length
# parameter - L: length of simulation box
# returns: vector of the three force components
function force_between_particles(📍1, 📍2, cut📏, L)
    r = nearest_image_displacement(📍1, 📍2, L)
    r📏 = norm(r)
    LJ_🤜(📏) = 48*📏^-13 - 24*📏^-7
    if cut📏 == -1
        🤜 = LJ_🤜(r📏)
    elseif r📏 >= cut📏
        🤜 = 0
    else
        🤜 = LJ_🤜(r📏) - LJ_🤜(cut📏)
    end
    return 🤜 / r📏 * r
end

# Calculates the LJ potential from the interaction between two particles.
#
# parameter - 📍1: [x, y, z] vector for first particle's position
# parameter - 📍2: [x, y, z] vector for second particle's position
# parameter - cut📏: cutoff length
# returns: scalar of LJ potential from interaction between the two particles.
function LJ_potential(📍1, 📍2, cut📏, L)
    r📏 = norm(nearest_image_displacement(📍1, 📍2, L))
    LJ_U(📏) = 4 * (📏^-12 - 📏^-6)
    if cut📏 == -1
        U = LJ_U(r📏)
    elseif r📏 >= cut📏
        U = 0
    else
        🤜_cut = (-48*cut📏^-13 + 24*cut📏^-7)
        U = LJ_U(r📏) - LJ_U(cut📏) - (r📏-cut📏)*🤜_cut
    end
    return U
end

# Calculates the pressure resultant from interaction between particles 1 and 2.
#
# parameter - 📍1: [x, y, z] vector for first particle's position
# parameter - 📍2: [x, y, z] vector for second particle's position
# parameter - L: length of simulation box
# parameter - F: force on particle 1 from particle 2
# returns: pressure of interaction
function pressure_between_particles(📍1, 📍2, L, F)
    r = nearest_image_displacement(📍1, 📍2, L)
    return r' * F
end

# Computes LJ forces using current positions
#
# parameter! - 🤜s: forces of all particles
# parameter - 📍s: positions of all particles
# parameter - 🧛: number of particles
# parameter - L: length of one edge of simulation box
# returns: total LJ potential energy of system
# returns: pressure as calculated from forces 
function LJ_🤜s_and_energy!(🤜s, 📍s, 🧛, cut📏, L)
    fill!(🤜s, 0);
    U = 0; P = 0;
    for i = 1:🧛           # for each particle
        for j = i+1:🧛     # for each particle that i interacts with
            F = force_between_particles(📍s[i,:], 📍s[j,:], cut📏, L)
            🤜s[i,:] .+= F
            🤜s[j,:] .-= F
            U += LJ_potential(📍s[i,:], 📍s[j,:], cut📏, L)
            P += pressure_between_particles(📍s[i,:], 📍s[j,:], L, F);
        end 
    end
    return U, P
end

# Updates velocities in-place by half a timestep for velocity Verlet.
#
# parameter! - 🚗s: vector of starting velocities
# parameter - 🤜s: vector of forces for each particle
# parameter - ζs: adjustment from thermostat
# parameter - ⏲️: timestep
function update_🚗s!(🚗s, 🤜s, ζs, ⏲️)
    🚗s .+= ⏲️/2 * (🤜s - ζs .* 🚗s)
end

# Updates positions in-place by one timestep for velocity Verlet.
#
# Moves each particle by its velocity times one timestep. After moving, moves
# particles back within the simulation bounds as dictated by PBCs.
#
# parameter! - 📍s: vector of starting positions
# parameter! - unadjusted📍s: vector of positions if not subject to PBCs
# parameter - 🚗s: vector of velocity for each particle
# parameter - ⏲️: timestep
# parameter - L: length of one edge of simulation box
function update_📍s!(📍s, unadjusted📍s, 🚗s, ⏲️, L)
    📍s .+= 🚗s*⏲️
    unadjusted📍s .+= 🚗s*⏲️
    # if any coordinate is negative, increase it by L. if any coordinate is 
    # beyond L, decrease that by L. All particles should remain within the box.
    📍s .+= L*((📍s .< 0) - (📍s .> L))
end

# Updates thermostat in-place by one timestep for velocity Verlet.
#
# parameter! - ζs: starting thermostat
# parameter - 🚗s: vector of velocity for each particle
# parameter - τ: damping timescale
# parameter - T_des: desired temperature
# parameter - ⏲️: timestep
# parameter - 🧛: number of particles in system
function update_ζs!(ζs, 🚗s, τ, T_des, ⏲️, 🧛)
    T_inst = sum(🚗s.^2) / (3 * (🧛-1))
    ζs .+= (⏲️ / τ^2) * (T_inst / T_des - 1)
end

# runs VV for one timestep; modified positions and velocities in-place.
#
# parameter! - 📍s: vector of positions for each particle
# parameter! - 🚗s: vector of velocities for each particle
# parameter! - 🤜s: vector of forces on each particle
# parameter! - ζs: thermostat constants
# parameter - ⏲️: timestep
# parameter - L: length of one side of simulation box
# parameter - cut📏: cutoff radius
# parameter - 🧛: number of particles in system
# parameter - τ: thermostat damping timescale
# parameter - T_des: desired temperature
# returns: system total potential energy at the end of timestep
# returns: pressure from forces term at end of timestep
function vv_one_timestep!(p)
    update_🚗s!(p.🚗s, p.🤜s, p.ζs, p.⏲️)
    update_📍s!(p.📍s, p.unadjusted📍s, p.🚗s, p.⏲️, p.boxLength)
    update_ζs!(p.ζs, p.🚗s, p.τ, p.T_des, p.⏲️, p.🧛)
    U, P_from_🤜s = LJ_🤜s_and_energy!(p.🤜s, p.📍s, p.🧛, p.cut📏, p.boxLength);
    update_🚗s!(p.🚗s, p.🤜s, p.ζs, p.⏲️)
    p.🚗s ./= (1 .+ 0.5.*p.⏲️.*p.ζs)
    return U, P_from_🤜s
end

end