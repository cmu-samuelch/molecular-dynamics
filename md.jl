using Plots, LinearAlgebra, Random, Statistics, Dates

using .ReadWrite, .VelocityVerlet, .RecordData

# initializes velocities to a certain average
#
# parameter - 📍s: number of particles
# parameter - μ: average velocity
# parameter - 🌡️: desired temperature of system
# returns - 🚗s: vector of velocities
function init_velocities(📍s, μ, 🌡️)
    🚗s = randn(size(📍s))
    🚗s[end,:] = -sum(🚗s[1:end-1, :], dims=1)
    🚗s .*= 🌡️; 🚗s .+= μ;
    return 🚗s
end

# Simulates particles.
# 
# parameter - 📍s: starting positions
# parameter - 🚗s: starting velocities
# parameter - ⏲️: timestep.
# parameter - cut📏: cutoff radius
# parameter - L: length of one side of the simulation box
# parameter - duration: timesteps to simulate for.
# parameter - 📭: location where positions get dumped
# parameter - resolution: number of timesteps between each write to outfile
# returns - 📨: table with all data we want to keep track of.
function simulate(📍s, 🚗s, ⏲️, cut📏, L, duration, 📭, resolution)
    🧛 = size(📍s)[1]
    📨 = zeros(duration, 8)
    write(📭, "")
    📭_stream = open(📭, "a")

    write_xyz_frame(📭_stream, 📍s, 0, resolution)
    🤜s = zeros(size(📍s));

    _ = LJ_🤜s_and_energy!(🤜s,📍s, 🧛, cut📏, L);
    for i = 1:duration
        # VV forward one timestep
        U, P_from_🤜s = vv_one_timestep!(📍s, 🚗s, 🤜s, ⏲️, L, cut📏, 🧛)
        
        # generate some data to plot later
        t = i*⏲️; K = calculate_kinetic(🚗s)
        🌡️, P = calculate_🌡️_and_P(🚗s, 🧛, L^3, P_from_🤜s)
        📨[i,:] = [t K U sum(🚗s, dims=1) 🌡️ P]
        
        # write current positions to outfile as one frame
        write_xyz_frame(📭_stream, 📍s, i, resolution)

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
    resolution = 10
    cut📏 = 2.5
    L = 6.8
    🌡️ = 0.8
    📭 = "256test0.xyz"

    📍s = read_📩(📩)
    🚗s = init_velocities(📍s, 0, 🌡️)
    
    println("done!")
    t0 = now();
    println("[", t0, "] running MD...")
    data = simulate(📍s, 🚗s, 0.004, cut📏, L, 50000, 📭, resolution)
    println(now() - t0, " elapsed during MD simulation")

    write_data(data, "256diag0.csv")

    p_H = plot(data[:,1], [data[:,2:3] sum(data[:,2:3], dims=2)], labels=["K" "U" "H"], legend=:left, xlabel="time", ylabel="energy")
    p_p = plot(data[:,1], data[:,4:6], labels=["p_x" "p_y" "p_z"], xlabel="time", ylabel="momentum")
    p_T = plot(data[:,1], data[:,7], legend=false, xlabel="time", ylabel="temperature")
    p_P = plot(data[:,1], data[:,8], legend=false, xlabel="time", ylabel="pressure")
    plot(p_H, p_p, p_T, p_P)

end

main()