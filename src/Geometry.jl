abstract type Geometry{dim, T} end

function Geometry(::Type{S}, coordinates::AbstractVector{<: Vec{dim, T}}, args...) where {S <: Geometry, dim, T}
    q = one(Quaternion{T})
    S(coordinates, q, args...)
end

realtype(x::Geometry{<: Any, T}) where {T} = T
num_coordinates(x::Geometry) = length(coordinates(x))
coordinates(x::Geometry) = x.coordinates
coordinates(x::Geometry, i::Int) = (@_propagate_inbounds_meta; coordinates(x)[i])
quaternion(x::Geometry) = x.q
moment_of_inertia(x::Geometry) = throw(ArgumentError("$(typeof(x)) is not supported yet."))

"""
    translate!(geo::Geometry, u::Vec)

Translate `geo` by the displacement `u`.
"""
function translate!(geometry::Geometry, u::Vec)
    coords = coordinates(geometry)
    @inbounds @simd for i in eachindex(coords)
        coords[i] += u
    end
    _translate!(geometry, u)
    geometry
end
# for additional treatment in subtypes
_translate!(geometry::Geometry, u::Vec) = nothing

"""
    rotate!(geo::Geometry{2}, θ::Real)
    rotate!(geo::Geometry{3}, θ::Vec)

Rotate `geo` by the angle `θ`.
In 3D, `normalize(θ)` and `norm(θ)` should represent the rotation axis and the angle (radian), respectively.
"""
function rotate!(geometry::Geometry{dim}, θ::Vec{3}, xc::Vec{dim} = centroid(geometry)) where {dim}
    # https://www.ashwinnarayan.com/post/how-to-integrate-quaternions/
    coords = coordinates(geometry)
    q = exp(Quaternion(θ/2))
    @inbounds @simd for i in eachindex(coords)
        coords[i] = xc + Tensorial.resizedim(rotate(coords[i] - xc, q), Val(dim))
    end
    geometry.q = q * geometry.q
    _rotate!(geometry, θ)
    geometry
end
rotate!(geometry::Geometry{2}, θ::Real, x::Vec{2} = centroid(geometry)) = rotate!(geometry, Vec(0,0,θ), x)
# for additional treatment in subtypes
_rotate!(geometry::Geometry, θ::Vec) = nothing

isinside(geo::Geometry) = Base.Fix2(isinside, geo)

Broadcast.broadcastable(geo::Geometry) = (geo,)

function Base.show(io::IO, mime::MIME"text/plain", x::Geometry)
    print(io, typeof(x), ":\n")
    print(io, "  Centroid: ", centroid(x), "\n")
    print(io, "  Quaternion: ", quaternion(x))
end
