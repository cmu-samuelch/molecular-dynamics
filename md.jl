using Plots, LinearAlgebra, Random, Statistics, Dates

include("ReadWrite.jl"); include("Parameters.jl")
include("RecordData.jl"); include("VelocityVerlet.jl");
using .ReadWrite, .VelocityVerlet, .RecordData, .Parameters

# Simulates particles.
# 
# parameter - p: struct containing simulation parameters
# parameter - 📭: path to file where frames will be stored
# returns - 📨: table with all data we want to keep track of.
function simulate(p, 📭)
    📨 = zeros(Float64, p.numTimesteps, 9)
    write(📭, "")
    📭_stream = open(📭, "a")

    write_xyz_frame(📭_stream, p.📍s, 0, p.frameSaveFrequency)

    📍s0 = copy(p.📍s)

    _ = LJ_🤜s_and_energy!(p.🤜s, p.📍s, p.🧛, p.cut📏, p.boxLength);
    for i = 1:p.numTimesteps
        # VV forward one timestep
        U, P_from_🤜s = vv_one_timestep!(p)
        
        # generate some data to plot later
        t = i*p.⏲️; K = calculate_kinetic(p.🚗s)
        🌡️, P = calculate_🌡️_and_P(p.🚗s, p.🧛, p.boxLength^3, P_from_🤜s)
        MSD = calculateMSD(p.📍s, 📍s0)
        📨[i,:] = [t K U sum(p.🚗s, dims=1) 🌡️ P MSD]
        
        # write current positions to outfile as one frame
        write_xyz_frame(📭_stream, p.📍s, i, p.frameSaveFrequency)

        if i % (p.numTimesteps/10) == 0
            println("simulation ", i/p.numTimesteps*100, "% complete; ",
                (p.numTimesteps-i), "/", p.numTimesteps, " timesteps remaining")
        end
    end

    return 📨
end

# runs everything for the current problem set.
function main()
    print("initializing...")
    # PARAMETERS TO CHANGE
    filename = "4.1"
    📭 = "256p" * filename * ".xyz"
    p = setup("liquid256.txt")
    
    println("done!")

    t0 = now();
    println("[", t0, "] running MD...")
    data = simulate(p, 📭)
    println(now() - t0, " elapsed during MD simulation")

    write_data(data, "256p" * filename * ".csv")

    t = data[:,1]
    p_H = plot(t, [data[:,2:3] sum(data[:,2:3], dims=2)], labels=["K" "U" "H"], legend=:left, xlabel="time", ylabel="energy")
    p_p = plot(t, data[:,4:6], labels=["p_x" "p_y" "p_z"], xlabel="time", ylabel="momentum")
    p_T = plot(t, data[:,7], legend=false, xlabel="time", ylabel="temperature")
    p_P = plot(t, data[:,8], legend=false, xlabel="time", ylabel="pressure")
    p_MSD = plot(t, data[:,9], legend=false, xlabel="time", ylabel="msd")
    plot(p_H, p_p, p_T, p_P)

end

main()