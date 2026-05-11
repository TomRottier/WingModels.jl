"""
    cosine_spacing(η)

Cosine spacing for chordwise point η. Creats a higher clustering of points near the leading and trailing edge so that it can have better resolution.

`η` is the chordwise coordinate [0,1] where 0 is the leading edge and 1 is the trailing edge.

Returns a new chordwise coordinate that is cosine spaced.
"""
cosine_spacing(η) = 0.5 - 0.5cos(π * η)

"""
    chordwise_coordinates(nchord)

Generate a vector of `nchord` chordwise coordinates from 0 to 1 using cosine spacing.

To create the full aerofoil use:
    [chordwise_coordinates(nchord ÷ 2); reverse(chordwise_coordinates(nchord ÷ 2))]

"""
function chordwise_coordinates(nchord)
    xs = range(0, 1, length=nchord)

    return cosine_spacing.(xs)
end
