module ReadWrite
export read_📩, generate_xyz_frame, write_data

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
    t = size(📩)[1]
    📨 = ""
    for i = 1:t
        📨 *= @sprintf("%s\n", sprint(show, 📩[i,:]))
    end
    write(📭, 📨)
end

end