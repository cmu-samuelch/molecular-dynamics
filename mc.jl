using Plots, LinearAlgebra, Random, Statistics, Dates

include("ReadWrite.jl"); include("Parameters.jl"); include("MCMonteCarlo.jl")
using .ReadWrite, .Parameters, .MCMC

# Simulates particles.
# 
# parameter - p: struct containing simulation parameters
# parameter - 📭: path to file where frames will be stored
# returns - 📨: table with all data we want to keep track of.
function simulate(p, 📭)
    📨 = zeros(Float64, p.numTimesteps, 4)
    write(📭, "")
    📭_stream = open(📭, "a")

    write_xyz_frame(📭_stream, p.📍s, 0, p.frameSaveFrequency)

    equilibrated = false

    β = 1/p.T_des
    max_🫨 = 0.15

    U = U_system(p.📍s, p.🧛, p.cut📏, p.boxLength);
    for i = 1:p.numTimesteps
        # MC one trial
        U, P, accept = MCMCtrial!(p.📍s, U, β, p.🧛, p.cut📏, p.boxLength, max_🫨)
        # generate some data to plot later
        t = i*p.⏲️;
        📨[i,:] = [t U P accept]

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
    filename = "6.2"
    📭 = "256p" * filename * ".xyz"
    p = setup("liquid256.txt")
    
    println("done!")

    t0 = now();
    println("[", t0, "] running MD...")
    data = simulate(p, 📭)
    println(now() - t0, " elapsed during MD simulation")

    write_data(data, "256p" * filename * ".csv")

    t = data[:,1]
    p_U = plot(t, data[:,2], xlabel="trials", ylabel="⟨U_potential⟩", legend=false)
    savefig("energy.png")
    p_P = plot(t, data[:,3], xlabel="trials", ylabel="⟨P⟩", legend=false)
    savefig("pressure.png")
    plot(p_U, p_P)

    println(mean(data, dims=1))
end

main()