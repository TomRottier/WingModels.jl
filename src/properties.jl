export area, mean_chord, mean_aerodynamic_chord, aspect_ratio, second_moment_of_area

"""
    mean_chord(pl::AbstractPlanform)
Calculate the (standard) mean chord of a wing as follows:
        
        int_0^b chord(ξ, pl) dξ
(the integral of the chord length over the wing length divided by the wing length (1)).

The answer is in normalised units so can be thought of as the mean chord in multiples of the wing length.
"""
mean_chord(pl::AbstractPlanform) = quadgk(ξ -> chord(ξ, pl), 0, 1)[1]



"""
    mean_aerodynamic_chord(pl::AbstractPlanform)
Calculate the mean aerodynamic chord of a wing as follows:

        int_0^b chord(ξ, pl)² dξ / int_0^b chord(ξ, pl) dξ
(the integral of the chord length squared divided by the integral of the chord length over the wing length).

The answer is in normalised units so can be thought of as the mean aerodynamic chord in multiples of the wing length.

The mean aerodynamic chord is defined such that the total moment on the wing is equivalent to that on a rectangular wing of chord length equal to the mean aerodynamic chord and of equal span (assuming the moment coefficient is constant along the span).
"""
mean_aerodynamic_chord(pl::AbstractPlanform) = quadgk(ξ -> chord(ξ, pl)^2, 0, 1)[1] / area(pl)



"""
    area(pl::AbstractPlanform)
Calculate the projected area of a wing as follows:
    
    int_0^1 chord(ξ, pl) dξ

(the integral from 0 to 1 of the chord length at each normalised spanwise location `ξ`).

The answer is in normalised units so can be thought of as the area in multiples of the wing length squared. For a wing with normalised span, the area is the same as the mean chord.

"""
area(pl::AbstractPlanform) = mean_chord(pl) # area is same as mean chord for wing with normalised span




"""
    aspect_ratio(pl::AbstractPlanform)
Calculate the aspect ratio of a wing as defined as:

    aspect_ratio = span² / area

As these wings are normalised by their span, the aspect ratio is just 1 / area which is equal to 1 / c̄    

"""
aspect_ratio(pl::AbstractPlanform) = 1 / area(pl)



"""
    second_moment_of_area(pl::AbstractPlanform)
Calculate the second moment of area of a wing about a chordwise axis at the root.

    second_moment_of_area = int_A ξ² dξ dη = int_0^1 int_0^c(ξ) ξ² dη dξ = int_0^1 ξ^2 * chord(ξ, pl) dξ

"""
second_moment_of_area(pl::AbstractPlanform) = quadgk(ξ -> ξ^2 * chord(ξ, pl), 0.0, 1.0)[1]


mean_chord(w::Wing) = mean_chord(w.planform)
mean_aerodynamic_chord(w::Wing) = mean_aerodynamic_chord(w.planform)
area(w::Wing) = area(w.planform)
aspect_ratio(w::Wing) = aspect_ratio(w.planform)
second_moment_of_area(w::Wing) = second_moment_of_area(w.planform)