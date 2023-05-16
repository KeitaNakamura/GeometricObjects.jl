mutable struct FinitePoints{dim, T} <: Geometry{dim, T}
    coordinates::Vector{Vec{dim, T}}
    q::Quaternion{T}
end

FinitePoints(coords::AbstractVector) = Geometry(FinitePoints, coords)

centroid(geo::FinitePoints) = mean(coordinates(geo))

function moment_of_inertia(geo::FinitePoints)
    coords = coordinates(geo)
    sum(x -> (x⋅x)*I-symmetric(x⊗x,:U), coords) / length(coords)
end
