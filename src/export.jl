export write_stl, write_pts, get_conns

#=
binary stl specification:
    - 80 byte header, generally ignored but should not contain `solid` to be confused with ascii stl [80 UINT8]
    - 4 byte unsigned integer (little endian), number of triangles [UINT32]
    - triangle data, each one described by 12 32-bit floating point numbers (little endian)
        - 3 32-bit floating point numbers, normal vector [3 REAL32]
        - 3 32-bit floating point numbers, vertex 1 [3 REAL32]
        - 3 32-bit floating point numbers, vertex 2 [3 REAL32]
        - 3 32-bit floating point numbers, vertex 3 [3 REAL32]
        - 2 byte unsigned integer spacer (always 0) [UINT16]
=#

"""
    write_stl(filename, wing, nchord; sf=1.0)

Write the vector of points contained in `wing` to an STL file. You must supply the `nchord` used to create `wing` so the facet connections can be determined correctly. Keyword argument `sf` is the scale factor and will scale each point in the output by the given amount.
"""
function write_stl(filename, wing, nchord; sf=1.0)
    file = splitext(filename)[1] * ".stl"

    open(file, "w") do io
        # header
        write(io, Vector{UInt8}("UNITS=mm"))
        write(io, zeros(UInt8, 72))

        # connection indices for triangles
        nspan = length(wing) ÷ nchord
        surface_conns = get_conns(nchord, nspan)
        root_face_conns = mapreduce(vcat, 1:nchord÷2-1) do i
            return [(i, i + 1, nchord - (i - 1)), (nchord - (i - 1), nchord - i, i + 1)]
        end
        tip_face_conns = [x .+ nchord * (nspan - 1) for x in root_face_conns]
        connections = [surface_conns; root_face_conns; tip_face_conns]


        # number of triangles
        N = length(connections) |> UInt32
        write(io, N)

        # triangle data
        for i in 1:N
            # normal
            v1, v2, v3 = wing[[connections[i]...]] .* sf # vertices of triangle
            x1 = v2 - v1
            x2 = v3 - v1
            n = [x1[2] * x2[3] - x1[3] * x2[2], x1[3] * x2[1] - x1[1] * x2[3], x1[1] * x2[2] - x1[2] * x2[1]]
            n ./= sqrt(sum(abs2, n))

            # write to file
            write(io, Float32.(n)...)
            write(io, Float32.(v1)...)
            write(io, Float32.(v2)...)
            write(io, Float32.(v3)...)

            # attribute byte count
            write(io, zeros(UInt16, 1))
        end
    end

end


function write_stl(filename, w::Wing; nchord, nspan, sf=1.0)
    wing_pts = wing(w, 0.0, 1.0; nchord, nspan)
    return write_stl(filename, wing_pts, 2nchord; sf)
end

# get indices of the points that make up a triangle
function vertex_conn(i, nchord)
    if rem(i, nchord) != 0
        [(i, i + nchord, i + 1 + nchord), (i + 1 + nchord, i + 1, i)]
    else
        [(i, i + nchord, i + 1), (i + 1, i + 1 - nchord, i)]

    end
end

"""
    get_conns(nchord, nspan)
Return a vector of indices indicating the index into the vector of points of the wing the vertices which make up a triangle/facet. For example, for `nchord = 10` the first two elements will be `(1,2,11), (11,12,2)` which indicate that index `1`, `2`, and `12` index in `wing` will be the vertices forming the first triangular face for a mesh and likewise `11`,`12`,`2` will be the second.

This function can be used to interface with `Meshes.jl`. To convert a `wing` defined in this package to a `Meshes.jl` mesh:

```jldoctest
# convert points
julia> pts = Meshes.Point3.(wing)

# facet connections, nchord and nspan used to create wing
julia> connections = connect.(get_conn(nchord, nspan))

# create mesh
julia> m = SimpleMesh(pts, connections)

# visualise mesh
julia> viz(m)
```

"""
get_conns(nchord, nspan) = reduce(vcat, [vertex_conn(pt + ((i - 1) * nchord), nchord) for i in 1:nspan-1 for pt in 1:nchord])

"""
    write_pts(filename, wing; delim=" ")

Write each point in `wing` to a .txt file with delimiter `delim`.    
"""
function write_pts(filename, wing; delim=" ")
    file = splitext(filename)[1] * ".txt"

    open(file, "w") do io
        for pt in wing
            println(io, join(pt, delim))
        end
    end
end


function write_pts(filename, w::Wing; nchord, nspan, delim=" ")
    wing_pts = wing(w, 0.0, 1.0; nchord, nspan)
    return write_pts(filename, wing_pts; delim)
end