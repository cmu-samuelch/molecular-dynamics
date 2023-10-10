using Plots, LinearAlgebra, Random, Statistics, Dates

using .ReadWrite, .VelocityVerlet

# initializes velocities to a certain average
#
# parameter - 📍s: number of particles
# parameter - μ: average velocity
# parameter - 🌡️: desired temperature of system
# returns - 🚗s: vector of velocities
function init_velocities(📍s, μ, 🌡️)
    🚗s = zeros(size(📍s))
    🚗s[1:end-1, :] = randn(size(🚗s)[1]-1, size(🚗s)[2])
    🚗s[end,:] = -sum(🚗s, dims=1)
    🚗s = (🚗s .* 🌡️) .+ μ
    return 🚗s
end

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
# returns: instantaneous average temperature and pressure for the system
function calculate_🌡️_and_P(🚗s, 🧛, V)
    🚗s_squared_mean = sum(🚗s.^2)
    🌡️ = 🚗s_squared_mean / (3 * (🧛-1))
    P = 🧛 * 🌡️ / V
    return 🌡️, P
end

# Simulates particles.
# 
# parameter - 📍s: starting positions
# parameter - 🚗s: starting velocities
# parameter - ⏲️: timestep.
# parameter - duration: timesteps to simulate for.
# parameter - 📭: location where positions get dumped
# parameter - cut📏: cutoff radius
# parameter - resolution: number of timesteps between each time frame is written
#                         to the .xyz output file 
# parameter - L: length of one side of the simulation box
# returns - 📨: table with columns containing timesteps, K, U, and p-components.
function simulate(📍s, 🚗s, ⏲️, cut📏, L, duration, 📭, resolution)
    🧛 = size(📍s)[1]
    📨 = zeros(duration, 8)
    📭_stream = open(📭, "a")

    frame = generate_xyz_frame(📍s, 0)
    write(📭_stream, frame)

    🤜s = zeros(size(📍s));
    _ = LJ_🤜s_and_energy!(🤜s,📍s, 🧛, cut📏, L);
    for i = 1:duration
        # VV forward one timestep
        U = vv_one_timestep!(📍s, 🚗s, 🤜s, ⏲️, L, cut📏, 🧛)
        
        # generate some data to plot later
        t = i*⏲️; K = calculate_kinetic(🚗s)
        🌡️, P = calculate_🌡️_and_P(🚗s, 🧛, L^3)
        📨[i,:] = [t K U sum(🚗s, dims=1) 🌡️ P]
        
        # write current positions to outfile as one frame
        if i % resolution == 0
            frame = generate_xyz_frame(📍s, i)
            write(📭_stream, frame)
        end

        if i % (duration/25) == 0
            println("simulation ", i/duration*100, "% complete; ",
                     (duration-i), "/", duration, " timesteps remaining")
        end
    end

    return 📨
end

# runs everything for the current problem set.
function main()
    print("initializing...")
    # PARAMETERS TO CHANGE
    📩 = "liquid256.txt"
    resolution = 1
    cut📏 = 2.5
    L = 6.8
    🌡️ = 0.8

    📍s = read_📩(📩)
    🚗s = init_velocities(📍s, 0, 🌡️)
    println(calculate_🌡️_and_P(🚗s, 256, 6.8^3))
    
    📭 = "pset-3-2.xyz"
    
    write(📭, "")
    println("done!")
    t0 = now();
    println("[", t0, "] running MD...")
    data = simulate(📍s, 🚗s, 0.004, cut📏, L, 20000, 📭, resolution)
    println(now() - t0, " elapsed during MD simulation")

    write_data(data, "diagnostic.csv")

    p_H = plot(data[:,1], [data[:,2:3] sum(data[:,2:3], dims=2)], labels=["K" "U" "H"], legend=:left, xlabel="time", ylabel="energy")
    p_p = plot(data[:,1], data[:,4:6], labels=["p_x" "p_y" "p_z"], xlabel="time", ylabel="momentum")
    p_T = plot(data[:,1], data[:,7], legend=false, xlabel="time", ylabel="temperature")
    p_P = plot(data[:,1], data[:,8], legend=false, xlabel="time", ylabel="pressure")
    plot(p_H, p_p, p_T, p_P)

end

main()