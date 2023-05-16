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
function attitude(x::Geometry{dim, T}) where {dim, T}
    v = rotate(Vec{3,T}(1,0,0), quaternion(x))
    @Tensor v[1:dim]
end
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
    geometry
end

"""
    rotate!(geo::Geometry{2}, θ::Real)
    rotate!(geo::Geometry{3}, θ::Vec)

Rotate `geo` by the angle `θ`.
In 3D, `normalize(θ)` and `norm(θ)` should represent the rotation axis and the angle (radian), respectively.
"""
function rotate! end
rotate!(geometry::Geometry{3}, θ::Vec{3}) = _rotate!(geometry, θ)
rotate!(geometry::Geometry{2}, θ::Real) = _rotate!(geometry, Vec(0,0,θ))
function _rotate!(geometry::Geometry{dim}, θ::Vec) where {dim}
    # https://www.ashwinnarayan.com/post/how-to-integrate-quaternions/
    coords = coordinates(geometry)
    q = exp(Quaternion(θ/2))
    xc = centroid(geometry)
    @inbounds @simd for i in eachindex(coords)
        coords[i] = xc + Tensorial.resizedim(rotate(coords[i] - xc, q), Val(dim))
    end
    geometry.q = q * geometry.q
    geometry
end

isinside(geo::Geometry) = Base.Fix2(isinside, geo)

Broadcast.broadcastable(geo::Geometry) = (geo,)

function Base.show(io::IO, mime::MIME"text/plain", x::Geometry)
    print(io, typeof(x), ":\n")
    print(io, "  Coordinates: [", join(coordinates(x), ", "), "]\n")
    print(io, "  Attitude: ", attitude(x))
end
