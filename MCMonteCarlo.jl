module MCMC
include("Parameters.jl")
using LinearAlgebra, Random, .Parameters

export MCMCtrial!, systemUandP

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

# Computes U_potential for a single particle
#
# parameter - 📍s: positions of all particles
# parameter - i: index of particle to perturb
# parameter - 🧛: number of particles
# parameter - cut📏: cutoff radius
# parameter - L: length of one edge of simulation box
# returns: potential from the selected particle's interactions
function Uparticle(📍s, i, 🧛, cut📏, L)
    U = 0
    for j = 1:🧛
        if j != i
            U += LJ_potential(📍s[i,:], 📍s[j,:], cut📏, L)
        end
    end
    return U
end

# Computes total U_potential and pressure from interactions of system using LJ
#
# parameter - 📍s: positions of all particles
# parameter - 🧛: number of particles
# parameter - cut📏: cutoff radius
# parameter - L: length of one edge of simulation box
# return - U: total U_potential for the system
# return - P: total pressure from interparticle interactions for the system
function systemUandP(📍s, 🧛, cut📏, L)
    U, P = 0, 0
    for i = 1:🧛           # for each particle
        for j = i+1:🧛     # for each particle that i interacts with
            U += LJ_potential(📍s[i,:], 📍s[j,:], cut📏, L)
            F = force_between_particles(📍s[i,:], 📍s[j,:], cut📏, L)
            P += pressure_between_particles(📍s[i,:], 📍s[j,:], L, F);
        end
    end
    return U, P / (3*L^3)
end

# Perturbs a single particle, and accepts the perturbation as appropriate.
#
# parameter! - 📍s: positions of all particles
# parameter - i: index of particle to perturb
# parameter - β: (k_B T)^(-1), determines perturbation probability
# parameter - 🧛: number of particles
# parameter - cut📏: cutoff radius
# parameter - L: length of one edge of simulation box
# parameter - 🫨max: maximum perturbation in any dim
# returns: whether this perturbation was accepted or rejected
function 🫨1!(📍s, i, β, 🧛, cut📏, L, 🫨max)
    U = Uparticle(📍s, i, 🧛, cut📏, L)
    🫨 = (rand(Float64, 3) .- 0.5) .* 🫨max;        # generate a perturbation
    📍s[i,:] += 🫨                                  # apply it
    📍s .+= L*((📍s .< 0) - (📍s .> L))             # PBC
    🫨U = Uparticle(📍s, i, 🧛, cut📏, L)
    # reject and revert the perturbation if U increases and stat check failed
    if 🫨U > U && (rand(Float64) > exp(-β*(🫨U-U)))
        📍s[i,:] .-= 🫨
        📍s .-= L*((📍s .< 0) - (📍s .> L)) # undoes PBC
        return 0
    end
    return 1
end

# Carries out one MCMC trial on the entire system.
#
# parameter! - 📍s: positions of all particles
# parameter - U: total potential energy of current state
# parameter - β: (k_B T)^(-1), determines perturbation probability
# parameter - 🧛: number of particles
# parameter - cut📏: cutoff radius
# parameter - L: length of one edge of simulation box
# parameter - 🫨max: maximum perturbation in any dim
# returns: energy of the new state
# returns: pressure of the new state
function MCMCtrial!(📍s, U, β, 🧛, cut📏, L, 🫨max)
    order = shuffle(1:256)
    accept = 0
    for i = order
        accept += 🫨1!(📍s, i, β, 🧛, cut📏, L, 🫨max)
    end
    V = L^3
    U, Pints = systemUandP(📍s, 🧛, cut📏, L)
    P = 🧛 / (β*V) + Pints
    return U, P, accept / 🧛
end

end