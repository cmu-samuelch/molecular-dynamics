# GLOBAL ASSUMPTIONS: epsilon, sigma, and m are all 1.
# hopefully, this will hold true and I won't have to go back and insert these
# everywhere.

# TODO
# randomly initialize particle velocities with zero total momentum
# implement continuous force/energy with cutoff of 2.5 (dimless)
# calculates instantaneous temperature, pressure
# applies periodic boundary conditions and the nearest-image convention with
#   side length as a variable set in the code

using Plots, Printf, LinearAlgebra

# Reads the contents of the file into a N-by-3 array of positions.
#
# parameter - infile_path: path to the file to read
# returns: N-by-3 array of positions
function read_infile(infile_path)
    text = read(infile_path, String)
    lines = split(text, "\n")
    output = Array{Float64}(undef, length(lines)-1, 3)
    for i = eachindex(lines)
        if lines[i] != ""
            vals = split(lines[i])
            for j = 1:3
                output[i, j] = parse(Float64, vals[j])
            end
        end
    end
    return output
    end

# Writes positions in current state to xyz format
#
# parameter - rs: positions to record
# parameter - i: frame number
# returns: string of all positions
function generate_xyz_frame(rs, i)
    n = size(rs)[1]
    text = @sprintf("%i\nFrame %i\n", n, i)
    for i = 1:n
        text *= @sprintf("a %f %f %f\n", rs[i,1], rs[i,2], rs[i,3])
    end
    return text
end

# Writes data to a file.
#
# parameter - data: matrix of data to store
# parameter - outfile_path: location to store data
function write_data(data, outfile_path)
    (t, cols) = size(data)
    output = ""
    for i = 1:t
        for j = 1:cols
            output *= @sprintf("%f,", data[i,j])
        end
        output *= "\n"
    end
    write(outfile_path, output)
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
# parameter - rs: positions of all particles
# parameter - N: number of particles
# returns: array of forces on each particle
# returns: total LJ potential energy of system
function LJ_forces_and_energy(rs, N)
    Fs = zeros(size(rs))
    U = 0
    for i = 1:N           # for each particle
        for j = i+1:N     # for each particle that i interacts with
            F = force_between_particles(rs[i,:], rs[j,:])
            Fs[i,:] += F
            Fs[j,:] -= F
            U += LJ_potential(rs[i,:], rs[j,:])
        end
    end
    return Fs, U
end

function update_v(v, force, dt)
    v += force * dt/2
    return v
end

function update_r(rs, vs, dt)
    rs += vs*dt
    return rs
end

function calculate_kinetic(vs)
    return sum(vs.^2) / 2
end

# Simulates molecules.
# 
# parameter - rs: starting positions
# parameter - vs: starting velocities
# parameter - N: number of particles
# parameter - timestep: timestep.
# parameter - duration: timesteps to simulate for.
# parameter - outfile: location where positions get dumped
# returns: table with columns containing timesteps, K, U, and p-components.
function simulate(rs, vs, N, timestep, duration, outfile)
    output = zeros(duration, 6)
    outfile_stream = open(outfile, "a")

    Fs, _ = LJ_forces_and_energy(rs, N);
    for i = 1:duration
        
        # VV forward one timestep
        vs = update_v(vs, Fs, timestep)
        rs = update_r(rs, vs, timestep)
        Fs, U = LJ_forces_and_energy(rs, N);
        vs = update_v(vs, Fs, timestep)
        
        # generate some data to plot later
        t = i*timestep; K = calculate_kinetic(vs)
        output[i,:] = [t K U sum(vs[:,1]) sum(vs[:,2]) sum(vs[:,3])]
        
        # write current positions to outfile as one frame
        frame = generate_xyz_frame(rs, i)
        println(frame)
        write(outfile_stream, frame)
    end

    return output
end


function main()
    rs = read_infile("pset-2/10.txt")
    vs = zeros(size(rs))
    N = size(rs)[1]

    outfile = "dump.xyz"

    write(outfile, "")
    data = simulate(rs, vs, N, 0.002, 1000, outfile)

    write_data(data, "diagnostic.csv")

end

main()