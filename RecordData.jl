module RecordData
export calculate_kinetic, calculate_🌡️_and_P

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

end