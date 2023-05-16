mutable struct FinitePoints{dim, T} <: Geometry{dim, T}
    coordinates::Vector{Vec{dim, T}}
    q::Quaternion{T}
end

FinitePoints(coords::AbstractVector) = Geometry(FinitePoints, coords)

function centroid(geo::FinitePoints{dim, T}) where {dim, T}
    xc = zero(Vec{dim, T})
    @inbounds @simd for x in coordinates(geo)
        xc += x
    end
    xc / num_coordinates(geo)
end

function moment_of_inertia(geo::FinitePoints{dim, T}) where {dim, T}
    Ic = zero(SymmetricSecondOrderTensor{dim, T})
    @inbounds @simd for x in coordinates(geo)
        Ic += (x⋅x)*I - symmetric(x⊗x, :U)
    end
    Ic / num_coordinates(geo)
end
