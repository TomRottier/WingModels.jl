module WingModels

using LinearAlgebra: normalize, ×
using QuadGK: quadgk

export AbstractPlanform, AbstractAerofoil, Wing
export quarter_chord, chord, aerofoil, aerofoil_pt, leading_edge, trailing_edge, wing_planform, wing

"""
    AbstractPlanform

Supertype for all planforms. A planform is defined here as having a quarter chord line and a spanwise chord distribution.
    
To use your own planform with the functions in this package you must define your own subtype of `AbstractPlanform` and define the following methods:

    struct YourPlanform <: AbstractPlanform
        data # data needed to define your planform
    end

    quarter_chord(ξ, pl::YourPlanform<:AbstractPlanform) = ...    
which gives the x (chordwise) coordinate of quarter chord for given spanwise location `ξ`.

    chord(ξ, pl::YourPlanform<:AbstractPlanform) = ... 
which gives the chord length for the given spanwise location `ξ`.
"""
abstract type AbstractPlanform end

"""
    AbstractAerofoil

Supertype for all aerofoil distributions. An aerofoil distribution is here defined as a function which which a vector of 3D points describing an aerofoil for a given spanwise location.
    
To use your own aerofoil distribution with the functions in this package you must define your own subtype of `AbstractAerofoil` and the following method:

    struct YourAerofoil <: AbstractAerofoil
        data # data needed to define your aerofoil distrbutions
    end

    aerofoil(ξ, af::YourAerofoil<:AbstractAerofoil; nchord) = ...

which returns a Vector{Vector{<:Real}} describing a series of points outlining the aerofoil at the given spanwise location `ξ`. The points must go from leading edge to trailing edge along upper surface then back to leading edge along lower surface. The number of points is defined by the keyword argument `nchord`.

You may optionally define the following method:

    aerofoil_pt(η, ξ, af::YourAerofoil<:AbstractAerofoil; upper=true) = ...

which returns a 2D vector for a point (x and z coordinate) on the aerofoil at chordwise location `η` and spanwise location `ξ`. 'upper' determines is this point is on the upper or lower surface of the aerofoil.
"""
abstract type AbstractAerofoil end

# generic functions
"""
    leading_edge(ξ, pl::AbstractPlanform)

Calculate the x (chordwise) coordinate for a given planform `pl` at the given spanwise location `ξ`.
"""
leading_edge(ξ, pl::AbstractPlanform) = quarter_chord(ξ, pl) - 0.25chord(ξ, pl)

"""
    trailing_edge(ξ, pl::AbstractPlanform)

Calculate the x (chordwise) coordinate for a given planform `pl` at the given spanwise location `ξ`.
"""
trailing_edge(ξ, pl::AbstractPlanform) = quarter_chord(ξ, pl) + 0.75chord(ξ, pl)

"""
    wing_planform(pl::AbstractPlanform; n=100)

Generate a vector of points giving the x (chordwise) coordinate of each of `n` points on a given planform `pl`. Points start from wing root leading edge and proceed to wing tip leading edge, wing tip trailing edge and back to wing root trailing edge.
"""
wing_planform(pl::AbstractPlanform; n=100) = [
    map(ξ -> leading_edge(ξ, pl), range(0, 1, length=n ÷ 2));
    map(ξ -> trailing_edge(ξ, pl), range(1, 0, length=n ÷ 2))
]

"""
    wing(pl::AbstractPlanform, af::AbstractAerofoil; nchord=100, nspan=50)
Generate a wing from a given planform and aerofoil distribution. 

A wing is a vector of 3D points on the surface of the wing. The points start from the leading edge of the root aerofoil and go over the upper surface to the trailing edge then along lower surface back to the leading edge before going to the next aerofoil along the span. There are `nchord` points per aerofoil and `nspan` aerofoils along the span, giving a total of `nchord * nspan` points on the wing. Specify `ξ₀` and `ξ₁` as the start and end spanwise locations, respectively.

To use this function with your own defined planform and aerofoil distribution you must define the relevant functions.

See also: [`AbstractPlanform`](@ref), [`AbstractAerofoil`](@ref) 
"""
function wing(pl::AbstractPlanform, af::AbstractAerofoil, ξ₀=0.0, ξ₁=1.0; nchord=100, nspan=50)
    # spanwise locations
    ξs = range(ξ₀, ξ₁, length=nspan)

    # create wing
    wing = mapreduce(vcat, ξs) do ξ
        # create normalised aerofoil
        normalised_aerofoil = aerofoil(ξ, af; nchord=nchord)

        # fit aerofoils to planform - must scale first, then translate
        return translate_aerofoil(ξ, scale_aerofoil(ξ, normalised_aerofoil, pl), pl)
    end

    return wing
end

# scale aerofoil by chord length in x and z components, leaidng edge must be at (0.0, 0.0)
scale_aerofoil(ξ, aerofoil, pl::AbstractPlanform) = map(pt -> pt .* (chord(ξ, pl), 1, chord(ξ, pl)), aerofoil)

# translate aerofoil to quarter chord line - aerofoils all have leading edge at 0, translate to 0.25 (applies translation before scale)
translate_aerofoil(ξ, aerofoil, pl::AbstractPlanform) = aerofoil .+ Ref([quarter_chord(ξ, pl) - 0.25chord(ξ, pl), 0.0, 0.0])


# wing type
"""
    Wing

Container to hold the <:AbstractPlanform and <:AbstractAerofoil for a wing.
"""
struct Wing{A<:AbstractAerofoil,P<:AbstractPlanform}
    aerofoil::A
    planform::P
end
wing(w::Wing, ξ₀=0.0, ξ₁=1.0; nchord=100, nspan=50) = wing(w.planform, w.aerofoil, ξ₀, ξ₁; nchord=nchord, nspan=nspan)
chord(ξ, w::Wing) = chord(ξ, w.planform)

aerofoil(ξ, w::Wing; nchord) = aerofoil(ξ, w.aerofoil; nchord=nchord)
aerofoil_pt(η, ξ, w::Wing; upper) = aerofoil_pt(η, ξ, w.aerofoil; upper=upper)
quarter_chord(ξ, w::Wing) = quarter_chord(ξ, w.planform)
leading_edge(ξ, w::Wing) = leading_edge(ξ, w.planform)
trailing_edge(ξ, w::Wing) = trailing_edge(ξ, w.planform)
wing_planform(w::Wing; n=100) = wing_planform(w.planform; n=n)

# fallback to error if interface not defined properly
const NON_DEFINED_INTERFACE_MSG = "reaching here means you haven't defined the appropiate functions for the interface" # I think 
chord(ξ, p::AbstractPlanform) = leading_edge(ξ, p) - trailing_edge(ξ, p) #error(NON_DEFINED_INTERFACE_MSG)
quarter_chord(ξ, p::AbstractPlanform) = leading_edge(ξ, p) + 0.25chord(ξ, p) #error(NON_DEFINED_INTERFACE_MSG)
aerofoil(ξ, p::AbstractAerofoil; nchord=100) = error(NON_DEFINED_INTERFACE_MSG)
aerofoil_pt(η, ξ, p::AbstractAerofoil; upper, nchord=100) = error(NON_DEFINED_INTERFACE_MSG)

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
