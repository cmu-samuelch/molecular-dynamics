module RecordData
export calculate_kinetic, calculate_🌡️_and_P, calculateMSD

# Calculates instantaneous total kinetic energy in the system.
#
# parameter - 🚗s: vector of velocities for each particle
# returns: sum of kinetic energy for the entire system at current time
function calculate_kinetic(🚗s)
    return sum(🚗s.^2) / 2
end

# Calculates instantaneous temperature and pressure in the system.
#
# parameter - 🚗s: vector of velocities for each particle
# parameter - 🧛: number of particles in system
# parameter - V: volume of the system
# parameter - P_from_🤜s: pressure as calculated from forces 
# returns: instantaneous average temperature and pressure for the system
function calculate_🌡️_and_P(🚗s, 🧛, V, P_from_🤜s)
    🌡️ = sum(🚗s.^2) / (3 * (🧛-1))
    P = 🧛 * 🌡️ / V + P_from_🤜s / (3*V)
    return 🌡️, P
end

# calculates mean squared displacement of positions in system
#
# parameter - Δ📍s: vector of displacements from original positions
# parameter - 🧛: number of particles in system
# returns: instantaneous mean squared displacement for system
function calculateMSD(Δ📍s, 🧛)
    return sum((Δ📍s).^2) / 🧛
end

end