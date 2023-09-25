# GLOBAL ASSUMPTIONS: epsilon, sigma, and m are all 1.
# hopefully, this will hold true and I won't have to go back and insert these
# everywhere.

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
    r = r2 - r1
    r_len = norm(r)
    # println(r_len)
    force = 4*(12 * r_len^-13 - 6*r_len^-7)
    # println(force)
    return -force / r_len .* r
end

function LJ_potential(r1, r2)
    r_len = norm(r1 - r2)
    return 4 * (r_len^-12 - r_len^-6)
end

# Computes LJ forces using current positions
#
# parameter - rs: positions of all particles
# parameter - particle_ct: number of particles
# returns: array of forces on each particle
# returns: total LJ potential energy of system
function LJ_forces(rs, particle_ct)
    Fs = zeros(size(rs))
    U = 0
    for i = 1:particle_ct           # for each particle
        for j = i+1:particle_ct     # for each particle that i interacts with
            # collects forces
            Fs[i,:] = Fs[i,:] .+ force_between_particles(rs[i,:], rs[j,:])
            Fs[j,:] = Fs[j,:] .- force_between_particles(rs[i,:], rs[j,:])
            U = U + LJ_potential(rs[i,:], rs[j,:])
        end
    end
    return Fs, U
end

# Simulates molecules.
# 
# parameter - timestep: timestep.
# parameter - duration: timesteps to simulate for.
# output - (to file): positions at end of sim.
# plots - several.
# returns: nothing
function simulate(timestep, duration)
    rs = read_infile("10.txt")
    vs = zeros(size(rs))
    particle_ct = size(rs)[1]
    Ks = zeros(2, duration)
    Us = zeros(2, duration)
    ps = zeros(4, duration)
    Ks[1,:] = 1:duration
    for i = 1:duration
        t = i * timestep;

        # VV forward one timestep
        Fs, _ = LJ_forces(rs, particle_ct);
        v_one_half = vs .+ timestep/2 * Fs;
        new_rs = rs .+ timestep.*v_one_half;
        Fs, U = LJ_forces(new_rs, particle_ct);
        new_vs = v_one_half .+ timestep/2 * Fs;
        rs = new_rs;
        vs = new_vs;

        Ks[1,i] = t
        Us[1,i] = t
        ps[1,i] = t
        for i = 1:particle_ct
            Ks[2,i] = Ks[2,i] + norm(vs[i,:])^2/2
            Us[2,i] = Us[2,i] + U
            ps[2:4,i] = ps[2:4,i] + vs[i,:]
        end

        # if (i + 5) % 10 == 0
        #     out_path = @sprintf("dump%i.xyz", i)
        #     write_outfile(out_path, rs)
        # end
    end
    p = plot(ps[1,:], Ks[2,:] + Us[2,:])
    xlabel!("time")
    ylabel!("H")
    display(p)

    out_path = @sprintf("dump%i_final.xyz", duration)
    write_outfile(out_path, rs)
end

simulate(0.002, 1000)





# NOTES
# need N-by-3 arrays to store positions, velocities, forces

# please do not create N-by-N arrays (don't store all atomic separations/pair forces)

# please comment code :)

# background on mol sim file formats - will be useful later on
# for simple applications: XYZ is very simple to work with
# life saving - openBabel, avogadro for converting trajfile formats


# OUTLINE:
# - read 10.txt to initialize the positions
# - run VV algorithm
# - run for however many timesteps we need to run