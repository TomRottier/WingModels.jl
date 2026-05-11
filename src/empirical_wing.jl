# determine planform from empirical data
export EmpiricalPlanform


# expects some kind of interpolation object that can be called with ξ and return value
"""
    EmpiricalPlanform{T1,T2} <: AbstractPlanform

Create a planform from the outline of a planform.

    EmpiricalPlanform(le, te)

expects `le` and `te` to be functions that when called a spanwise position returns the x (chordwise) coordinate of the leading and trailing edges repsectively.
"""
struct EmpiricalPlanform{T1,T2} <: AbstractPlanform
    le::T1
    te::T2
end

chord(ξ, p::EmpiricalPlanform) = abs(p.le(ξ) - p.te(ξ))
quarter_chord(ξ, p::EmpiricalPlanform) = p.le(ξ) + 0.25chord(ξ, p)