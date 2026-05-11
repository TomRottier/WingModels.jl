export TrapezoidalPlanform, EllipticalPlanform, RectangularPlanform, TriangularPlanform
export RectangularAerofoil

"""
    TrapezoidalPlanform <: AbstractPlanform

A simplfied planform which defines a two jointed quarter chord line and a linear chord distribution. Equivalently, a rectangular upper wing, and triangular lower wing.

    TrapezoidalPlanform(r, ϕ, c₀, x)

Quarter chord defined by two jointed arm model with two parameters:
- `r`: relative length of upper arm
- `ϕ`: angle between upper and lower arm

Chord distribution defined by 2 parameters:
- `c₀`: root chord length normalised by wing length
- `x`: spanwise location where chord distribution begins to decrease linearly

The slope of the linear decrease is determined by the value of `x` such that the chord length is zero when ξ = 1.
"""
struct TrapezoidalPlanform <: AbstractPlanform
    r::Float64
    ϕ::Float64
    c₀::Float64
    x::Float64
end

quarter_chord(ξ, p::TrapezoidalPlanform) = ξ < p.r ? 0.0 : (ξ - p.r) * tan(p.ϕ)
chord(ξ, p::TrapezoidalPlanform) = ξ ≤ p.x ? p.c₀ : p.c₀ - p.c₀ * (ξ - p.x) / (1 - p.x)


"""
    EllipticalPlanform <: AbstractPlanform

A planform defined by its leading and trailing edges. Each edge is described by an elliptical curve given by the following equation:
    (x^2/a^2) + (y^2/b^2) = 1
where a is the extent in the x (spanwise) direction and b is the extent in the y (chordwise) direction. The wings here are normalised by their length so a = 1 and b can be thought of as c₀/2, where c₀ is the ratio of root chord to wing length as used elsewhere. This is only true however when the leading and trailing edges are given by the same elliptical curve.

    EllipticalPlanform(c₁, c₂)

An `EllipticalPlanform` is defined by the following parameters:
- `c₁`: distance from origin to leading edge at wing root normalised by wing length
- `c₂`: distance from origin to trailing edge at wing root normalised by wing length

The root chord length is thus given by `c₀ = c₁ + c₂`
    

"""
struct EllipticalPlanform <: AbstractPlanform
    c₁::Float64
    c₂::Float64
end

__elliptical_leading_edge(ξ, a, b) = sqrt((1 - ξ^2 / a^2) * b^2)
__elliptical_trailing_edge(ξ, c, d) = -sqrt((1 - ξ^2 / c^2) * d^2)
chord(ξ, pl::EllipticalPlanform) = __elliptical_leading_edge(ξ, 1.0, pl.c₁) - __elliptical_trailing_edge(ξ, 1.0, pl.c₂)
quarter_chord(ξ, pl::EllipticalPlanform) = -__elliptical_leading_edge(ξ, 1.0, pl.c₁) + 0.25 * chord(ξ, pl)


"""
    RectangularPlanform <: AbstractPlanform

A rectangular planform with constant chord length along the span and a quarter chord line at 0.25 chord length. Equivalent to a `TrapezoidalPlanform` with `r = 0`, `ϕ = 0`, `x = 0`, `c₀ = c`.

    RectangularPlanform(c)

- `c`: chord length normalised by wing length
"""
struct RectangularPlanform <: AbstractPlanform
    c::Float64
end
chord(ξ, p::RectangularPlanform) = p.c
quarter_chord(ξ, p::RectangularPlanform) = 0.0

"""
    TriangularPlanform <: AbstractPlanform

A triangular planform (delta wing).

    TriangularPlanform(c₀, Λ)

- `c₀`: root chord length normalised by wing length
- `Λ`: leading edge sweep angle

"""
struct TriangularPlanform <: AbstractPlanform
    c₀::Float64
    Λ::Float64
end
__triangular_leading_edge(ξ, pl) = 0.25pl.c₀ - ξ * tan(pl.Λ)
chord(ξ, pl::TriangularPlanform) = pl.c₀ * (1 - ξ)
quarter_chord(ξ, pl::TriangularPlanform) = __triangular_leading_edge(ξ, pl) - 0.25chord(ξ, pl)


"""
    RectangularAerofoil <: AbstractAerofoil

A rectangular aerofoil e.g. for a flat plate.
    RectangularAerofoil(t)
- `t`: thickness of the aerofoil (constant along its chord)
"""
struct RectangularAerofoil <: AbstractAerofoil
    t::Float64
end
aerofoil(ξ, p::RectangularAerofoil; nchord=100) = [fill(p.t, nchord ÷ 2); fill(-p.t, nchord ÷ 2)]
aerofoil_pt(η, ξ, p::RectangularAerofoil; upper) = upper ? p.t : -p.t
