# GLOBAL ASSUMPTIONS: ε, σ, and m are all 1.
# hopefully, this will hold true and I won't have to go back and insert these
# everywhere.

# TODO
# PS #3:
# DONE - randomly initialize particle velocities with zero total momentum
# DONE - implement continuous force/energy with cutoff of 2.5 (dimless)
# DONE - calculates instantaneous temperature, pressure
# applies periodic boundary conditions and the nearest-image convention
#   create side length as a variable set in the code
#
# PS #4:
# ???
#
# LONG-TERM:
# IN PROGRESS - improve variable names

using Plots, Printf, LinearAlgebra, Random, Statistics, Dates

# Reads the contents of the file into a N-by-3 array of positions.
#
# parameter - 📩: path to the file to read
# returns: N-by-3 array of positions
function read_📩(📩)
    text = read(📩, String)
    lines = split(text, "\n")
    📨 = Array{Float64}(undef, length(lines)-1, 3)
    for i = eachindex(lines)
        if lines[i] != ""
            vals = split(lines[i])
            for j = 1:3
                📨[i, j] = parse(Float64, vals[j])
            end
        end
    end
    return 📨
    end

# Writes positions in current state to xyz format
#
# parameter - 📍s: positions to record
# parameter - i: frame number
# returns: string of all positions
function generate_xyz_frame(📍s, i)
    n = size(📍s)[1]
    text = @sprintf("%i\nFrame %i\n", n, i)
    for i = 1:n
        text *= @sprintf("a %f %f %f\n", 📍s[i,1], 📍s[i,2], 📍s[i,3])
    end
    return text
end

# Writes data to a file.
#
# parameter - 📩: matrix of data to store
# parameter - 📭: location to store data
function write_data(📩, 📭)
    (t, cols) = size(📩)
    📨 = ""
    for i = 1:t
        for j = 1:cols
            📨 *= @sprintf("%f,", 📩[i,j])
        end
        📨 *= "\n"
    end
    write(📭, 📨)
end

# initializes velocities to a certain average
#
# parameter - 📍s: number of particles
# parameter - μ: average velocity
# parameter - 🌡️: desired temperature of system
# returns - 🚗s: vector of velocities
function init_velocities(📍s, μ, 🌡️)
    🚗s = zeros(size(📍s))
    🚗s[1:end-1, :] = randn!(🚗s[1:end-1, :]) .* 🌡️
    🚗s[end,:] = -sum(🚗s, dims=1)
    🚗s .+= μ
    return 🚗s
end

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
    r = r .+ L*((r .< -L/2) .- (r .> L/2))
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

# Computes LJ forces using current positions
#
# parameter - 📍s: positions of all particles
# parameter - 🧛: number of particles
# parameter - L: length of one edge of simulation box
# returns: array of forces on each particle
# returns: total LJ potential energy of system
function LJ_🤜s_and_energy(📍s, 🧛, cut📏, L)
    🤜s = zeros(size(📍s))
    U = 0
    for i = 1:🧛           # for each particle
        for j = i+1:🧛     # for each particle that i interacts with
            F = force_between_particles(📍s[i,:], 📍s[j,:], cut📏, L)
            🤜s[i,:] += F
            🤜s[j,:] -= F
            U += LJ_potential(📍s[i,:], 📍s[j,:], cut📏, L)
        end 
    end
    return 🤜s, U
end

# Updates velocities by half a timestep for velocity Verlet.
#
# parameter - 🚗s: vector of starting velocities
# parameter - 🤜s: vector of forces for each particle
# parameter - ⏲️: timestep
# returns - 🚗s: vector of the new velocities
function update_🚗s(🚗s, 🤜s, ⏲️)
    🚗s += 🤜s * ⏲️/2
    return 🚗s
end

# Updates positions by one timestep for velocity Verlet.
#
# Moves each particle by its velocity times one timestep. After moving, moves
# particles back within the simulation bounds as dictated by PBCs.
#
# parameter - 📍s: vector of starting positions
# parameter - 🚗s: vector of velocity for each particle
# parameter - ⏲️: timestep
# parameter - L: length of one edge of simulation box
# returns - 📍s: vector of the new positions
function update_📍s(📍s, 🚗s, ⏲️, L)
    📍s += 🚗s*⏲️
    # if any coordinate is negative, increase it by L. if any coordinate is 
    # beyond L, decrease that by L. All particles should remain within the box.
    📍s = 📍s .+ L*((📍s .< 0) - (📍s .> L))
    return 📍s
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
    🚗s_squared_mean = mean(sum(🚗s.^2, dims=2))
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

    🤜s, _ = LJ_🤜s_and_energy(📍s, 🧛, cut📏, L);
    for i = 1:duration
        # VV forward one timestep
        🚗s = update_🚗s(🚗s, 🤜s, ⏲️)
        📍s = update_📍s(📍s, 🚗s, ⏲️, L)
        🤜s, U = LJ_🤜s_and_energy(📍s, 🧛, cut📏, L);
        🚗s = update_🚗s(🚗s, 🤜s, ⏲️)
        
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
    🌡️ = 1

    📍s = read_📩(📩)
    🚗s = init_velocities(📍s, [0 24 0], 🌡️)
    
    📭 = "pset-3-2.xyz"
    
    write(📭, "")
    println("done!")
    t0 = now();
    println("[", t0, "]", " running MD...")
    data = simulate(📍s, 🚗s, 0.01, cut📏, L, 10000, 📭, resolution)
    println(now() - t0, " elapsed during MD simulation")

    write_data(data, "diagnostic.csv")

    p_H = plot(data[:,1], [data[:,2:3] sum(data[:,2:3], dims=2)], labels=["K" "U" "H"], xlabel="time", ylabel="energy")
    p_p = plot(data[:,1], data[:,4:6], labels=["p_x" "p_y" "p_z"], xlabel="time", ylabel="momentum")
    p_T = plot(data[:,1], data[:,7], legend=false, xlabel="time", ylabel="temperature")
    p_P = plot(data[:,1], data[:,8], legend=false, xlabel="time", ylabel="pressure")
    plot(p_H, p_p, p_T, p_P)

end

main()