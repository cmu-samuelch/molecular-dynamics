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

    equilibrated = false

    _ = LJ_🤜s_and_energy!(p.🤜s, p.📍s, p.🧛, p.cut📏, p.boxLength);
    for i = 1:p.numTimesteps
        # VV forward one timestep
        if equilibrated
            # nve dynamics
            U, P_from_🤜s = nveVVtimestep!(p)
        else
            # run dynamics with Nosé-Hoover thermostat
            U, P_from_🤜s = nhVVtimestep!(p)
        end
        # generate some data to plot later
        t = i*p.⏲️; K = calculate_kinetic(p.🚗s)
        🌡️, P = calculate_🌡️_and_P(p.🚗s, p.🧛, p.boxLength^3, P_from_🤜s)
        MSD = calculateMSD(p.Δ📍s, p.🧛)
        📨[i,:] = [t K U sum(p.🚗s, dims=1) 🌡️ P MSD]
        
        if i > 100 && !equilibrated     # do eq check if not yet eq
            equilibrated = equilibrationcheck(📨[i-100:i, 7], p.T_des)
            if equilibrated
                println("reached equilibration at ", t, " non-dimensional time (",
                i, " timesteps)")
            end
        end
        # write current positions to outfile as one frame
        write_xyz_frame(📭_stream, p.📍s, i, p.frameSaveFrequency)
        
        if i % (p.numTimesteps/10) == 0
            println("simulation ", i/p.numTimesteps*100, "% complete; ",
                (p.numTimesteps-i), "/", p.numTimesteps, " timesteps remaining")
        end
    end

    return 📨
end

# runs the eq check.
#
# Returns whether ≥ 80% of the previous 100 temperature readings were within 5%
# of the equilibrium temperature.
function equilibrationcheck(data, Tdes)
    return sum(0.95*Tdes .< data .< 1.05*Tdes) >= 80
end

# runs everything for the current problem set.
function main()
    print("initializing...")
    # PARAMETERS TO CHANGE
    filename = "4.1-30k"
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
    plot(p_H, p_p, p_T, p_P, p_MSD)

end

main()