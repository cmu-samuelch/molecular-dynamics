# compute

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

