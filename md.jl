# GLOBAL ASSUMPTIONS: epsilon, sigma, and m are all 1.
# hopefully, this will hold true and I won't have to go back and insert these
# everywhere.

# TODO
# DONE - randomly initialize particle velocities with zero total momentum
# implement continuous force/energy with cutoff of 2.5 (dimless)
# calculates instantaneous temperature, pressure
# applies periodic boundary conditions and the nearest-image convention
#   create side length as a variable set in the code

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
# parameter - 📍: positions to record
# parameter - i: frame number
# returns: string of all positions
function generate_xyz_frame(📍, i)
    n = size(📍)[1]
    text = @sprintf("%i\nFrame %i\n", n, i)
    for i = 1:n
        text *= @sprintf("a %f %f %f\n", 📍[i,1], 📍[i,2], 📍[i,3])
    end
    return text
end

# Writes data to a file.
#
# parameter - data: matrix of data to store
# parameter - outfile_path: location to store data
function write_data(data, outfile_path)
    (t, cols) = size(data)
    📨 = ""
    for i = 1:t
        for j = 1:cols
            📨 *= @sprintf("%f,", data[i,j])
        end
        📨 *= "\n"
    end
    write(outfile_path, 📨)
end

# initializes velocities to a certain average
#
# parameter - N: number of particles
# parameter - μ: average velocity
# returns: vector of velocities
function init_velocities(N, μ)
    🚗 = zeros(N, 3)
    🚗[1:end-1, :] = randn!(zeros(N-1, 3))
    🚗[end,:] = -[sum(🚗[:,1]) sum(🚗[:,2]) sum(🚗[:,3])]
    🚗 .+= μ
    return 🚗
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

function LJ_potential(r1, r2)
    r_len = norm(r1 - r2)
    return 4 * (r_len^-12 - r_len^-6)
end

# Computes LJ forces using current positions
#
# parameter - 📍: positions of all particles
# parameter - N: number of particles
# returns: array of forces on each particle
# returns: total LJ potential energy of system
function LJ_forces_and_energy(📍, N)
    🫱 = zeros(size(📍))
    U = 0
    for i = 1:N           # for each particle
        for j = i+1:N     # for each particle that i interacts with
            F = force_between_particles(📍[i,:], 📍[j,:])
            🫱[i,:] += F
            🫱[j,:] -= F
            U += LJ_potential(📍[i,:], 📍[j,:])
        end
    end
    return 🫱, U
end

function update_🚗(v, force, dt)
    v += force * dt/2
    return v
end

function update_📍(📍, 🚗, dt)
    📍 += 🚗*dt
    return 📍
end

function calculate_kinetic(🚗)
    return sum(🚗.^2) / 2
end

# Simulates particles.
# 
# parameter - 📍: starting positions
# parameter - 🚗: starting velocities
# parameter - N: number of particles
# parameter - timestep: timestep.
# parameter - duration: timesteps to simulate for.
# parameter - outfile: location where positions get dumped
# returns: table with columns containing timesteps, K, U, and p-components.
function simulate(📍, 🚗, N, timestep, duration, outfile, resolution)
    📨 = zeros(duration, 6)
    outfile_stream = open(outfile, "a")

    🫱, _ = LJ_forces_and_energy(📍, N);
    for i = 1:duration
        if i % 1000 == 0
            println(i, "/", duration, " timesteps complete...")
        end
        
        # VV forward one timestep
        🚗 = update_🚗(🚗, 🫱, timestep)
        📍 = update_📍(📍, 🚗, timestep)
        🫱, U = LJ_forces_and_energy(📍, N);
        🚗 = update_🚗(🚗, 🫱, timestep)
        
        # generate some data to plot later
        t = i*timestep; K = calculate_kinetic(🚗)
        📨[i,:] = [t K U sum(🚗[:,1]) sum(🚗[:,2]) sum(🚗[:,3])]
        
        # write current positions to outfile as one frame
        if i % resolution == 0
            frame = generate_xyz_frame(📍, i)
            # println(📨[i,:])
            write(outfile_stream, frame)
        end
    end

    return 📨
end


function main()
    println("running MD...")
    # PARAMETERS TO CHANGE
    📩 = "liquid256.txt"
    resolution = 1

    📍 = read_📩(📩)
    N = size(📍)[1]
    🚗 = init_velocities(N, [0 1 0])

    outfile = "dump-pset-3.xyz"

    write(outfile, "")
    data = simulate(📍, 🚗, N, 0.002, 1000, outfile, resolution)

    write_data(data, "diagnostic.csv")

    p_H = plot(data[:,1], data[:,2:3])
    p_p = plot(data[:,1], data[:,4:6])
    plot(p_H, p_p)

end

main()