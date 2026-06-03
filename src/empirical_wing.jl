# determine planform from empirical data
export EmpiricalAerofoil, EmpiricalPlanform


"""
    smoothing_number(λ, data, d)

Calculate the number of data points to pass through as determined by the smoothing parameter `0 ≤ λ ≤ 1`. The minimum number, when `λ = 0` is given by `d + 1` where `d` is the degree of the interpolating polynomial. The maximum number, when `λ = 1` is given by `size(data, 1)` i.e. the number of data points.
"""
smoothing_number(λ, data, d) = round(Int, λ * (size(data, 1) - (d + 1)) + (d + 1))

"""
    scale_data(x::AbstractMatrix{<:Number})

`data` must be in form:
    - x[:, 1] is the chord/span position
    - x[:, 2] is the upper/leading surface
    - x[:, 3] is the lower/trailing surface

Rescales data such that chord/spanwise position goes from 0 to 1 and surface is normalised by the length.
"""
function scale_data(data::AbstractMatrix{<:Number})
    x0 = data[begin, 1]
    x1 = data[end, 1]
    l = x1 - x0
    out2 = data[:, 2] / l
    out3 = data[:, 3] / l
    out1 = (data[:, 1] .- x0) / l

    return hcat(out1, out2, out3)
end


"""
    EmpiricalAerofoil

Creates an aerofoil by interpolating through a given set of points using a cubic BSpline
"""
struct EmpiricalAerofoil{T1,T2,T3} <: AbstractAerofoil
    data::T1
    itp1::T2
    itp2::T3

    @doc """
        EmpiricalAerofoil(data, λ=1.0)

    `data` must be passed as an `n x 3` matrix of data points where:

    - `x[:, 1]` is the chordwise coordinate
    - `x[:, 2]` is the height on the upper surface at the corresponding chordwise coordinate
    - `x[:, 3]` is the height on the lower surface at the corresponding chordwise coordinate

    `λ` must be between 0 and 1, it controls the number of control points to interpolate through. `λ = 1` interpolates through all data points passed in `data`, whereas `λ = 0` interpolates through the minimum number of control points which is set by the degree of the piecewise polynomials (3).
    """
    function EmpiricalAerofoil(_data::T, λ=1.0) where T
        size(_data, 2) !== 3 && error("data must be an n x 3 matrix, currently size(data) = $(size(_data))")
        data = scale_data(_data)
        !(0.0 ≤ λ ≤ 1.0) && error("0 ≤ λ ≤ 1.0, currently: λ = $λ")
        d = 3 # degree of interpolating polynomial
        n = smoothing_number(λ, data, d)
        us_itp = BSplineApprox(data[:, 2], data[:, 1], d, n, :ArcLen, :Average)
        ls_itp = BSplineApprox(data[:, 3], data[:, 1], d, n, :ArcLen, :Average)

        return new{T,typeof(us_itp),typeof(ls_itp)}(data, us_itp, ls_itp)
    end
end
aerofoil_height(x, y, p::EmpiricalAerofoil; upper) = upper ? p.itp1(x) : p.itp2(x)


"""
    EmpiricalPlanform

Creates a planform by interpolating through a given set of points using a cubic BSpline
"""
struct EmpiricalPlanform{T1,T2,T3} <: AbstractPlanform
    data::T1
    itp1::T2
    itp2::T3

    @doc """
        EmpiricalPlanform(data, λ=1.0)

    `data` must be passed as an `n x 3` matrix of data points where:

    - `x[:, 1]` is the spanwise coordinate
    - `x[:, 2]` is the position of the leading edge at the corresponding spanwise coordinate
    - `x[:, 3]` is the position of the trailing edge at the corresponding spanwise coordinate

    `λ` must be between 0 and 1, it controls the number of control points to interpolate through. `λ = 1` interpolates through all data points passed in `data`, whereas `λ = 0` interpolates through the minimum number of control points which is set by the degree of the piecewise polynomials (3).
    """
    function EmpiricalPlanform(_data::T, λ=1.0) where T
        size(_data, 2) !== 3 && error("data must be an n x 3 matrix, currently size(data) = $(size(_data))")
        data = scale_data(_data)
        !(0.0 ≤ λ ≤ 1.0) && error("0 ≤ λ ≤ 1.0, currently: λ = $λ")
        d = 3 # degree of interpolating polynomial
        n = smoothing_number(λ, data, d)
        le_itp = BSplineApprox(data[:, 2], data[:, 1], d, n, :ArcLen, :Average)
        te_itp = BSplineApprox(data[:, 3], data[:, 1], d, n, :ArcLen, :Average)

        return new{T,typeof(le_itp),typeof(te_itp)}(data, le_itp, te_itp)
    end
end
_leading_edge(y, p::EmpiricalPlanform) = p.itp1(y)
_trailing_edge(y, p::EmpiricalPlanform) = p.itp2(y)

chord(y, p::EmpiricalPlanform) = p.itp2(y) - p.itp1(y)
quarter_chord(y, p::EmpiricalPlanform) = p.itp1(y) + 0.25chord(y, p)
