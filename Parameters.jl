module Parameters
include("ReadWrite.jl");
using .ReadWrite, Random
export setup, SimulationParameters
# initializes the simulation parameters
#
# parameter - inputpositions: path to file containing initial particle positions
# returns: struct containing parameters to pass to `simulate`.
function setup(inputpositions::String)
    📍📩 = inputpositions
    📍s = read_📩(📍📩, "\t")

    unadjusted = copy(📍s)
    🚗s = zeros(size(📍s))
    🤜s = zeros(size(📍s))
    ζs = zeros(size(📍s))

    timestep = 0.004
    numTimesteps = 10000
    resolution = 10
    cut📏 = 2.5
    L = 6.8

    ε = 1.66e-21
    T_des = 100
    T_des_nondimensional = T_des * 1.38e-23 / ε
    τdamp = 0.05

    # 🚗s = read_📩(🚗📩, ",")
    
    P = SimulationParameters(📍s, unadjusted, 🚗s, 🤜s, ζs, 
        T_des_nondimensional, size(📍s)[1], timestep, numTimesteps,
        cut📏, L, τdamp, resolution)

    initialize🚗s!(P.🚗s, 0, T_des_nondimensional)
    return P
end

# initializes velocities to a certain average
#
# parameter! - 🚗s: particle velocities
# parameter - μ: average velocity
# parameter - 🌡️: desired temperature of system
function initialize🚗s!(🚗s, μ, 🌡️)
    randn!(🚗s)
    🚗s[end,:] = -sum(🚗s[1:end-1, :], dims=1)
    🚗s .*= 🌡️; 🚗s .+= μ;
    return 🚗s
end

struct SimulationParameters         # keep track of parameters
    📍s::Array{Float64}             # positions
    unadjusted📍s::Array{Float64}   # unadjusted positions
    🚗s::Array{Float64}             # velocities
    🤜s::Array{Float64}             # forces
    ζs::Array{Float64}              # thermostat adjustments
    T_des::Float64                  # desired temperature
    🧛::Integer                     # particle count
    ⏲️::Float64                     # timestep
    numTimesteps::Int64
    cut📏::Float64                  # cutoff radius
    boxLength::Float64
    
    τ::Float64                      # damping timescale
    
    frameSaveFrequency::Int64       # how often to write frame to xyz file
end

# initializes the struct containing simulation output data
# function initializeoutput()

struct SimulationData               # where to output the sim's data
    data::Array{Float64}            # array of tracked data
    width::Int                      # how many params to keep track of
    params::Array{Function}         # array of functions for side calculations
    xyzDumpPath::String             # where to output the simulation frames
end

end