mutable struct FinitePoints{dim, T} <: Geometry{dim, T}
    coordinates::Vector{Vec{dim, T}}
    q::Quaternion{T}
    xc::Vec{dim, T}
    I::SymmetricSecondOrderTensor{dim, T}
end

FinitePoints(coords::AbstractVector) = Geometry(FinitePoints, coords, _centroid(coords), _moment_of_inertia(coords))

centroid(x::FinitePoints) = x.xc
function _centroid(coords::AbstractVector{Vec{dim, T}}) where {dim, T}
    xc = zero(Vec{dim, T})
    @inbounds @simd for x in coords
        xc += x
    end
    xc / length(coords)
end

moment_of_inertia(x::FinitePoints) = x.I
function _moment_of_inertia(coords::AbstractVector{Vec{dim, T}}) where {dim, T}
    Ic = zero(SymmetricSecondOrderTensor{dim, T})
    @inbounds @simd for x in coords
        Ic += (x⋅x)*I - symmetric(x⊗x, :U)
    end
    Ic / length(coords)
end

function _translate!(geo::FinitePoints, u::Vec)
    geo.xc += u
    geo
end

function _rotate!(geo::FinitePoints, θ::Vec)
    geo.I = rotate(geo.I, rotmat(geo.q))
    geo
end
