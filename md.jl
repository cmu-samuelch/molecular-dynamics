# GLOBAL ASSUMPTIONS: ε, σ, and m are all 1.
# hopefully, this will hold true and I won't have to go back and insert these
# everywhere.

# TODO
# PS #3:
# DONE - randomly initialize particle velocities with zero total momentum
# implement continuous force/energy with cutoff of 2.5 (dimless)
# calculates instantaneous temperature, pressure
# applies periodic boundary conditions and the nearest-image convention
#   create side length as a variable set in the code
#
# PS #4:
# ???
#
# LONG-TERM:
# IN PROGRESS - improve variable names

using Plots, Printf, LinearAlgebra, Random

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
# parameter - 🧛: number of particles
# parameter - μ: average velocity
# returns - 🚗s: vector of velocities
function init_velocities(🧛, μ)
    🚗s = zeros(🧛, 3)
    🚗s[1:end-1, :] = randn!(zeros(🧛-1, 3))
    🚗s[end,:] = -[sum(🚗s[:,1]) sum(🚗s[:,2]) sum(🚗s[:,3])]
    🚗s .+= μ
    return 🚗s
end

# Calculates the force vector exerted on particle 1 from LJ potential with 
# particle 2.
#
# parameter - r1: [x, y, z] vector for first particle's position
# parameter - r2: [x, y, z] vector for second particle's position
# returns: vector of the three force components
function force_between_particles(r1, r2)
    r = r1 - r2
    r_len = norm(r)
    force = 48*r_len^-13 - 24*r_len^-7
    return force / r_len * r
end

# Calculates the LJ potential from the interaction between two particles.
#
# parameter - r1: [x, y, z] vector for first particle's position
# parameter - r2: [x, y, z] vector for second particle's position
# returns: scalar of LJ potential from interaction between the two particles.
function LJ_potential(r1, r2)
    r_len = norm(r1 - r2)
    return 4 * (r_len^-12 - r_len^-6)
end

# Computes LJ forces using current positions
#
# parameter - 📍s: positions of all particles
# parameter - 🧛: number of particles
# returns: array of forces on each particle
# returns: total LJ potential energy of system
function LJ_🤜s_and_energy(📍s, 🧛)
    🤜s = zeros(size(📍s))
    U = 0
    for i = 1:🧛           # for each particle
        for j = i+1:🧛     # for each particle that i interacts with
            F = force_between_particles(📍s[i,:], 📍s[j,:])
            🤜s[i,:] += F
            🤜s[j,:] -= F
            U += LJ_potential(📍s[i,:], 📍s[j,:])
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
# parameter - 📍s: vector of starting positions
# parameter - 🚗s: vector of velocity for each particle
# parameter - ⏲️: timestep
# returns - 📍s: vector of the new positions
function update_📍s(📍s, 🚗s, ⏲️)
    📍s += 🚗s*⏲️
    return 📍s
end

# Updates positions by one timestep for velocity Verlet.
#
# parameter - 🚗s: vector of velocities for each particle
# returns: sum of kinetic energy for the entire system at current time
function calculate_kinetic(🚗s)
    return sum(🚗s.^2) / 2
end

# Simulates particles.
# 
# parameter - 📍s: starting positions
# parameter - 🚗s: starting velocities
# parameter - 🧛: number of particles
# parameter - ⏲️: timestep.
# parameter - duration: timesteps to simulate for.
# parameter - 📭: location where positions get dumped
# returns - 📨: table with columns containing timesteps, K, U, and p-components.
function simulate(📍s, 🚗s, 🧛, ⏲️, duration, 📭, resolution)
    📨 = zeros(duration, 6)
    📭_stream = open(📭, "a")

    🤜s, _ = LJ_🤜s_and_energy(📍s, 🧛);
    for i = 1:duration
        if i % 1000 == 0
            println(i, "/", duration, " timesteps complete...")
        end
        
        # VV forward one timestep
        🚗s = update_🚗s(🚗s, 🤜s, ⏲️)
        📍s = update_📍s(📍s, 🚗s, ⏲️)
        🤜s, U = LJ_🤜s_and_energy(📍s, 🧛);
        🚗s = update_🚗s(🚗s, 🤜s, ⏲️)
        
        # generate some data to plot later
        t = i*⏲️; K = calculate_kinetic(🚗s)
        📨[i,:] = [t K U sum(🚗s[:,1]) sum(🚗s[:,2]) sum(🚗s[:,3])]
        
        # write current positions to outfile as one frame
        if i % resolution == 0
            frame = generate_xyz_frame(📍s, i)
            # println(📨[i,:])
            write(📭_stream, frame)
        end
    end

    return 📨
end

# runs everything for the current problem set.
function main()
    println("running MD...")
    # PARAMETERS TO CHANGE
    📩 = "liquid256.txt"
    resolution = 1

    📍s = read_📩(📩)
    🧛 = size(📍s)[1]
    🚗s = init_velocities(🧛, [0 1 0])

    📭 = "dump-pset-3.xyz"

    write(📭, "")
    data = simulate(📍s, 🚗s, 🧛, 0.002, 1000, 📭, resolution)

    write_data(data, "diagnostic.csv")

    p_H = plot(data[:,1], data[:,2:3])
    p_p = plot(data[:,1], data[:,4:6])
    plot(p_H, p_p)

end

main()