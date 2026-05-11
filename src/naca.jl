# NACA00XX aerofoils
export NACA00XX, NACA4


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

# naca 00XX (symmetric) aerofoil
"""
    NACA00XX <: AbstractAerofoil
    
    NACA00XX(xx)

NACA 4 digit series symmetrical aerofoil with `xx` thickness as a percentage of chord. Max thickness occurs at 0.3c.

"""
struct NACA00XX <: AbstractAerofoil
    xx::Float64 # thickness as a percentage of chord
end

# point on surface of aerofoil
naca00_aerofoil_pt(η, ξ, xx; upper) = 5 * xx / 100 * naca4_aerofoil_thickness(η) * (2 * upper - 1)

function naca00_aerofoil(ξ, xx; N=100)
    ηs = chordwise_coordinates(N ÷ 2)
    upper = map(η -> [η, ξ, naca00_aerofoil_pt(η, ξ, xx; upper=true)], ηs)
    lower = map(η -> [η, ξ, naca00_aerofoil_pt(η, ξ, xx; upper=false)], reverse(ηs))

    return [upper; lower]
end
aerofoil(ξ, p::NACA00XX; nchord=100) = naca00_aerofoil(ξ, p.t; N=nchord)
aerofoil_pt(η, ξ, p::NACA00XX; upper) = naca00_aerofoil_pt(η, ξ, p.t; upper)

"""
    NACA4 <: AbstractAerofoil

    NACA4(m, p, xx)

NACA 4 digit series aerofoil with `m` maximum camber (in 10ths of chord), located at `p` along chord, and with maximum thickness `xx`.
"""
struct NACA4 <: AbstractAerofoil
    m::Float64 # maximum camber
    p::Float64 # position of maximum camber
    xx::Float64 # thickness
end

function naca4_aerofoil_pt(η, ξ, m, p, xx; upper)
    # thickness
    t = xx / 100
    yt = 5t * naca4_aerofoil_thickness(η)

    # camber
    m /= 100
    p /= 10
    yc = 0.0 ≤ η ≤ p ? naca4_aerofoil_camber_front(η, m, p) : p ≤ η ≤ 1.0 ? naca4_aerofoil_camber_back(η, m, p) : error("η out of range")
    dyc_dx = 0.0 ≤ η ≤ p ? naca4_aerofoil_camber_gradient_front(η, m, p) : p ≤ η ≤ 1.0 ? naca4_aerofoil_camber_gradient_back(η, m, p) : error("η out of range")

    # height above camber line
    θ = atan(dyc_dx)
    if upper
        x = η - yt * sin(θ)
        y = yc + yt * cos(θ)
    else
        x = η + yt * sin(θ)
        y = yc - yt * cos(θ)
    end

    return [x, y]
end
function naca4_aerofoil(ξ, m, p, xx; N=100)
    ηs = chordwise_coordinates(N ÷ 2)
    upper = map(ηs) do η
        tmp = naca4_aerofoil_pt(η, ξ, m, p, xx; upper=true)
        return [tmp[1], ξ, tmp[2]]
    end
    lower = map(reverse(ηs)) do η
        tmp = naca4_aerofoil_pt(η, ξ, m, p, xx; upper=false)
        return [tmp[1], ξ, tmp[2]]
    end

    return [upper; lower]
end

aerofoil(ξ, p::NACA4; nchord=100) = naca4_aerofoil(ξ, p.m, p.p, p.xx; N=nchord)
aerofoil_pt(η, ξ, p::NACA4; upper) = naca4_aerofoil_pt(η, ξ, p.m, p.p, p.xx; upper)