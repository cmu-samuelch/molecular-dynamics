# GLOBAL ASSUMPTIONS: epsilon, sigma, and m are all 1.
# hopefully, this will hold true and I won't have to go back and insert these
# everywhere.

# TODO: fix the output generator so i can get a gif of my particles moving around

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

# Writes positions to a file.
#
# parameter - outfile_path: path to write to
function write_outfile(outfile_path, rs)
    n = size(rs)[1]
    text = "10\n\n"
    for i = 1:n
        text = text * @sprintf("a %f %f %f\n", rs[i,1], rs[i,2], rs[i,3])
    end
    write(outfile_path, text)
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
            Fs[i,:] = Fs[i,:] + F
            Fs[j,:] = Fs[j,:] - F
            U = U + LJ_potential(rs[i,:], rs[j,:])
        end
    end
    return Fs, U
end

function update_v(v, force, dt)
    v = v + force * dt/2
    return v
end

function update_r(rs, vs, dt)
    rs = rs + vs*dt
    return rs
end

function calculate_kinetic(vs)
    return sum(vs.^2) / 2
end

# Simulates molecules.
# 
# parameter - timestep: timestep.
# parameter - duration: timesteps to simulate for.
# output - (to file): positions at end of sim.
# returns: table with columns containing timesteps, K, U, and p-components.
function simulate(rs, vs, N, timestep, duration)
    output = zeros(duration, 6)

    Fs, _ = LJ_forces_and_energy(rs, N);
    for i = 1:duration
        
        # VV forward one timestep
        vs = update_v(vs, Fs, timestep)
        rs = update_r(rs, vs, timestep)
        Fs, U = LJ_forces_and_energy(rs, N);
        vs = update_v(vs, Fs, timestep)
        
        t = i*timestep; K = calculate_kinetic(vs)
        output[i,:] = [t K U sum(vs[:,1]) sum(vs[:,2]) sum(vs[:,3])]
        
        # println(ps[i,2:4])
        # if (i + 5) % 10 == 0
        #     out_path = @sprintf("dump%i.xyz", i)
        #     write_outfile(out_path, rs)
        # end
    end

    out_path = @sprintf("dump%i_final.xyz", duration)
    write_outfile(out_path, rs)

    return output
end


function main()
    rs = read_infile("pset-2/10.txt")
    vs = zeros(size(rs))
    N = size(rs)[1]
    data = simulate(rs, vs, N, 0.002, 1000)

    # momentum plot
    plot(data[:,1], data[:,4:6], label=["p_x" "p_y" "p_z"], title="momentums for system are conserved (at zero)");
    # Hamiltonian plot
    plot(data[:,1], [data[:,2:3] data[:,2]+data[:,3]], label=["K" "U" "H"], title="Hamiltonian is conserved for this system");
end

main()