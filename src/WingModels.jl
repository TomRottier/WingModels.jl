module WingModels

using LinearAlgebra: normalize, ×
using QuadGK: quadgk

export AbstractPlanform, AbstractAerofoil, Wing
export quarter_chord, chord, aerofoil_height
export aerofoil, leading_edge, trailing_edge, planform, wing

const UNDEFINED_INTERFACE_MSG = "reaching here means you haven't defined the appropiate functions for the interface"


"""
    AbstractPlanform

Supertype for all planforms. A planform is defined here as having a quarter chord line and a spanwise chord distribution.
    
To use your own planform with the functions in this package you must define your own subtype of `AbstractPlanform` and define the following methods:

    struct YourPlanform <: AbstractPlanform
        data # data needed to define your planform
    end

    quarter_chord(y, pl::YourPlanform<:AbstractPlanform) = ...    
which gives the x (chordwise) coordinate of quarter chord for given spanwise location `y`.

    chord(y, pl::YourPlanform<:AbstractPlanform) = ... 
which gives the chord length for the given spanwise location `y`.
"""
abstract type AbstractPlanform end

chord(y, p::AbstractPlanform) = error(UNDEFINED_INTERFACE_MSG)
quarter_chord(y, p::AbstractPlanform) = error(UNDEFINED_INTERFACE_MSG)



"""
    AbstractAerofoil

Supertype for all aerofoil distributions. An aerofoil distribution is here defined as a function which which a vector of 3D points describing an aerofoil for a given spanwise location.
    
To use your own aerofoil distribution with the functions in this package you must define your own subtype of `AbstractAerofoil` and the following method:

    struct YourAerofoil <: AbstractAerofoil
        data # data needed to define your aerofoil distrbutions
    end

    aerofoil_height(x, y, af::YourAerofoil<:AbstractAerofoil; upper=true) = ...
which returns the height of the aerofoil at chordwise location `x` and spanwise location `y`. 'upper' determines if this point is on the upper or lower surface of the aerofoil.
"""
abstract type AbstractAerofoil end

aerofoil_height(x, y, p::AbstractAerofoil; upper) = error(UNDEFINED_INTERFACE_MSG)




"""
    aerofoil(y, p::AbstractAerofoil; n=100)

Returns `Vector` of the 2D [x,z] points describing the aerofoil at spanwise location `y`. The ordering of the points goes from the leading edge to trailing edge along the upper surface, the from the trailing edge to the leading edge along the lower surface. The keyword argument `n` gives the number of points along each upper and lower surface such that the total number of points returned is `2n`
"""
aerofoil(y, p::AbstractAerofoil; n=50) = [
    map(x -> [x, aerofoil_height(x, y, p; upper=true)], range(0, 1, length=n));
    map(x -> [x, aerofoil_height(x, y, p; upper=false)], range(1, 0; length=n))
]

"""
    leading_edge(y, pl::AbstractPlanform)

Returns the x (chordwise) coordinate for a given planform `pl` at the given spanwise location `y`.
"""
leading_edge(y, pl::AbstractPlanform) = quarter_chord(y, pl) - 0.25chord(y, pl)

"""
    trailing_edge(y, pl::AbstractPlanform)

Calculate the x (chordwise) coordinate for a given planform `pl` at the given spanwise location `y`.
"""
trailing_edge(y, pl::AbstractPlanform) = quarter_chord(y, pl) + 0.75chord(y, pl)

"""
    planform(pl::AbstractPlanform; n=100)

Returns a `Vector` of 2D [y,x] planform in the y-x plane (z = 0). The ordering of the points goes from the wing root at the leading edge to wing tip leading edge, wing tip trailing edge, and wing root trailing edge. The keyword argument `n` gives the number of points along each leading and trailing edge such that the total number of points returned is `2n`
"""
planform(pl::AbstractPlanform; n=50) = [
    map(y -> [y, leading_edge(y, pl)], range(0, 1, length=n));
    map(y -> [y, trailing_edge(y, pl)], range(1, 0, length=n))
]




"""
    Wing

Container to hold the <:AbstractPlanform and <:AbstractAerofoil for a wing.
"""
struct Wing{A<:AbstractAerofoil,P<:AbstractPlanform}
    aerofoil::A
    planform::P
end

quarter_chord(y, w::Wing) = quarter_chord(y, w.planform)
chord(y, w::Wing) = chord(y, w.planform)
aerofoil_height(x, y, w::Wing; upper) = aerofoil_height(x, y, w.aerofoil; upper=upper)
aerofoil(y, w::Wing; n=50) = aerofoil(y, w.aerofoil; n)
leading_edge(y, w::Wing) = leading_edge(y, w.planform)
trailing_edge(y, w::Wing) = trailing_edge(y, w.planform)
planform(w::Wing; n=50) = planform(w.planform; n)



"""
    wing(pl::AbstractPlanform, af::AbstractAerofoil; nchord=100, nspan=50)

Generates a wing from a given planform and aerofoil distribution.

Returns a `Vector` of 3D [x,y,z] points describing the surface of the wing. The ordering of the points goes from the wing root leading edge to the wing root trailing edge along the upper surface, then from the wing root trailing edge to the wing root leading edge along the lower surface. From here it moves outwards along the spanwise axis to the wing tip following the same ordering pattern. The keyword arguments `nchord` and `nspan` specify the number of points per aerofoil surface and per planform edge, respectively (i.e. there are `2nchord` points per aerofoil and `2nspan` aerofoils along the span giving a total to `4 * nchord * nspan` points per wing).

The other positional arguments `y₀ = 0.0` and `y₁ = 1.0` specify the starting and ending spanwise locations for the wing. This allows portions of the full wing to be easily created.

See also: [`AbstractPlanform`](@ref), [`AbstractAerofoil`](@ref) 
"""
function wing(pl::AbstractPlanform, af::AbstractAerofoil, y₀=0.0, y₁=1.0; nchord=100, nspan=50)
    # spanwise locations
    ys = range(y₀, y₁; length=nspan)

    # create wing
    wing = mapreduce(vcat, ys) do y
        # create normalised aerofoil
        normalised_aerofoil = aerofoil(y, af; n=nchord)

        # fit aerofoils to planform - must scale first, then translate
        return translate_aerofoil(y, scale_aerofoil(y, normalised_aerofoil, pl), pl)
    end

    return wing
end

wing(w::Wing, y₀=0.0, y₁=1.0; nchord=100, nspan=50) = wing(w.planform, w.aerofoil, y₀, y₁; nchord=nchord, nspan=nspan)


"""
    scale_aerofoil(y, aerofoil, pl::AbstractPlanform) 

Scale the `aerofoil` by the local chord length for the spanwise location `y`
"""
scale_aerofoil(y, aerofoil, pl::AbstractPlanform) = map(pt -> pt * chord(y, pl), aerofoil)

"""
    translate_aerofoil(y, aerofoil, pl::AbstractPlanform)

Translate the `aerofoil` at the spanwise location `y` such that its quarter chord is at quarter_chord(y, pl)` and its spanwise location is at `y`
"""
translate_aerofoil(y, aerofoil, pl::AbstractPlanform) = map(pt -> [pt[1] + leading_edge(y, pl), y, pt[2]], aerofoil)



include("utils.jl")
include("properties.jl")

# specific wing parameterisations
include("naca.jl")
include("empirical_wing.jl")
include("geometric.jl")

# liu wings
include("liu_wings.jl")
using .LiuWings
export LiuAerofoil, LiuPlanform, seagull, merganser, teal

# export
include("export.jl")

end # module WingModels
