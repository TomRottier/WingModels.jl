# NACA00XX aerofoils
export NACA4, NACA00


# # naca four digit aerofoil constants
# a₀ = 0.2969
# a₁ = -0.126
# a₂ = -0.3516
# a₃ = 0.2843
# a₄ = -0.1036

naca4_aerofoil_thickness(x) = 0.2969 * sqrt(x) - 0.126 * x - 0.3516 * x^2 + 0.2843 * x^3 - 0.1036 * x^4
naca4_aerofoil_camber_front(x, m, p) = m / p^2 * (2p * x - x^2)
naca4_aerofoil_camber_back(x, m, p) = m / (1 - p)^2 * (1 - 2p + 2p * x - x^2)
naca4_aerofoil_camber_gradient_front(x, m, p) = 2m / p^2 * (p - x)
naca4_aerofoil_camber_gradient_back(x, m, p) = 2m / (1 - p)^2 * (p - x)


"""
    NACA4 <: AbstractAerofoil

    NACA4(m, p, xx)

NACA 4 digit series aerofoil with `m` maximum camber (in 10ths of chord), located at `p` along chord, and with maximum thickness `xx`.

A convenience constructor can also be used to create symmetric NACA4 series aerofoils via NACA00(t).

!NOTE!: because NACA aerofoils are defined by a camber line and a thickness perpendicular to that camber line there is no explicit equation for the aerofoil height at chordwise location `x`. As such `aerofoil_height` returns the aerofoil height at a chordwise location given by `x - t * sin(θ)` where `t` and `θ` are the thickness and tangent slope at chordwise location `x`.
"""
struct NACA4 <: AbstractAerofoil
    m::Float64 # maximum camber
    p::Float64 # position of maximum camber
    xx::Float64 # thickness
end

NACA00(t) = NACA4(0, 0, t)

function naca4_aerofoil_height(x, y, m, p, xx; upper)
    # thickness
    t = xx / 100
    yt = 5t * naca4_aerofoil_thickness(x)

    # camber
    m /= 100
    p /= 10
    yc = 0.0 ≤ x ≤ p ? naca4_aerofoil_camber_front(x, m, p) : naca4_aerofoil_camber_back(x, m, p)
    dyc_dx = 0.0 ≤ x ≤ p ? naca4_aerofoil_camber_gradient_front(x, m, p) : naca4_aerofoil_camber_gradient_back(x, m, p)

    # height above camber line
    θ = atan(dyc_dx)
    if upper
        x = x - yt * sin(θ)
        y = yc + yt * cos(θ)
    else
        x = x + yt * sin(θ)
        y = yc - yt * cos(θ)
    end

    return [x, y]
end
aerofoil_height(x, y, p::NACA4; upper) = naca4_aerofoil_height(x, y, p.m, p.p, p.xx; upper)[2]

function naca4_aerofoil(y, m, p, xx; n=50)
    xs = chordwise_coordinates(n)
    upper = map(xs) do x
        tmp = naca4_aerofoil_height(x, y, m, p, xx; upper=true)
        return [tmp[1], tmp[2]]
    end
    lower = map(reverse(xs)) do x
        tmp = naca4_aerofoil_height(x, y, m, p, xx; upper=false)
        return [tmp[1], tmp[2]]
    end

    return [upper; lower]
end

aerofoil(y, p::NACA4; n=50) = naca4_aerofoil(y, p.m, p.p, p.xx; n)
# aerofoil_pt(x, y, p::NACA4; upper) = naca4_aerofoil_pt(x, y, p.m, p.p, p.xx; upper)


# # naca 00XX (symmetric) aerofoil
# """
#     NACA00XX <: AbstractAerofoil

#     NACA00XX(xx)

# NACA 4 digit series symmetrical aerofoil with `xx` thickness as a percentage of chord. Max thickness occurs at 0.3c.

# """
# struct NACA00XX <: AbstractAerofoil
#     xx::Float64 # thickness as a percentage of chord
# end

# # point on surface of aerofoil
# naca00_aerofoil_pt(x, y, xx; upper) = 5 * xx / 100 * naca4_aerofoil_thickness(x) * (2 * upper - 1)

# function naca00_aerofoil(y, xx; N=100)
#     xs = chordwise_coordinates(N ÷ 2)
#     upper = map(x -> [x, y, naca00_aerofoil_pt(x, y, xx; upper=true)], xs)
#     lower = map(x -> [x, y, naca00_aerofoil_pt(x, y, xx; upper=false)], reverse(xs))

#     return [upper; lower]
# end
# aerofoil(y, p::NACA00XX; nchord=100) = naca00_aerofoil(y, p.t; N=nchord)
# aerofoil_pt(x, y, p::NACA00XX; upper) = naca00_aerofoil_pt(x, y, p.t; upper)
