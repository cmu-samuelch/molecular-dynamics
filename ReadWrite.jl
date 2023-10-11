module ReadWrite
using Printf
export read_📩, write_xyz_frame, write_data

# Reads the contents of the file into an array of floats.
#
# parameter - 📩: path to the file to read
# parameter - separator: separator used to separate values in infile
# returns: M-by-N array of numeric values.
function read_📩(📩, separator)
    text = read(📩, String)
    lines = split(text, "\n")
    line1 = split(lines[1], separator)
    📨 = Array{Float64}(undef, length(lines)-1, length(line1))
    for i = eachindex(lines)
        if lines[i] != ""
            vals = split(lines[i], separator)
            for j = eachindex(vals)
                📨[i, j] = parse(Float64, vals[j])
            end
        end
    end
    return 📨
    end

# Writes positions in current state to xyz format, if applicable
#
# parameter - 📭: location to store data
# parameter - 📍s: positions to record
# parameter - i: frame number
# parameter - resolution: number of timesteps between each write to outfile
# returns: string of all positions
function write_xyz_frame(📭, 📍s, i, resolution)
    if i % resolution == 0
        n = size(📍s)[1]
        text = @sprintf("%i\nFrame %i\n", n, i)
        for i = 1:n
            text *= @sprintf("a %f %f %f\n", 📍s[i,1], 📍s[i,2], 📍s[i,3])
        end

        write(📭, text)
    end
end

# Writes data to a file.
#
# parameter - 📩: matrix of data to store
# parameter - 📭: location to store data
function write_data(📩, 📭)
    t = size(📩)[1]
    📨 = ""
    for i = 1:t
        📨 *= @sprintf("%s\n", sprint(show, 📩[i,:])[2:end-1])
    end
    write(📭, 📨)
end

end