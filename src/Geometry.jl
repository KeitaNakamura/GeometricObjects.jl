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
    q = quaternion(x)
    v = rotate(Vec{3,T}(1,0,0), q)
    @Tensor v[1:dim]
end
moment_of_inertia(x::Geometry) = throw(ArgumentError("$(typeof(x)) is not supported yet."))

@generated function copy_geometry(geometry::Geometry, coordinates::SVector, q::Quaternion)
    exps = [name in (:coordinates, :q) ? name : :(geometry.$name) for name in fieldnames(geometry)]
    quote
        $geometry($(exps...))
    end
end

function centered(x::Geometry)
    copy_geometry(x, coordinates(x) .- centroid(x), quaternion(x))
end

function enlarge(geometry::Geometry, R::Real)
    copy_geometry(
        geometry,
        centroid(geometry) .+ R * coordinates(centered(geometry)),
        quaternion(geometry),
    )
end

"""
    translate(g::Geometry, u::Vec)

Translate `g` by the displacement `u`.
"""
function translate(geometry::Geometry, u::Vec)
    copy_geometry(geometry, coordinates(geometry) .+ u, quaternion(geometry))
end

"""
    rotate(g::Geometry{2}, θ::Real)
    rotate(g::Geometry{3}, θ::Vec)

Rotate `g` by the angle `θ`.
In 3D, `normalize(θ)` and `norm(θ)` should represent the rotation axis and the angle (radian), respectively.
"""
function rotate end
rotate(geometry::Geometry{3}, θ::Vec{3}) = _rotate(geometry, θ)
rotate(geometry::Geometry{2}, θ::Real) = _rotate(geometry, Vec(0,0,θ))
function _rotate(geometry::Geometry{dim}, θ::Vec) where {dim}
    # https://www.ashwinnarayan.com/post/how-to-integrate-quaternions/
    q = exp(Quaternion(θ/2))
    xc = centroid(geometry)
    copy_geometry(
        geometry,
        (@. xc + $Tensorial.resizedim(rotate($coordinates(geometry) - xc, q), Val(dim))),
        q * quaternion(geometry),
    )
end

for IterType in (:Tuple, :AbstractArray)
    @eval function Base.findall(pred::Base.Fix2{typeof(in), <: Geometry}, iter::$IterType)
        findall(x -> in(x, pred.x), iter)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", x::Geometry)
    print(io, typeof(x), ":\n")
    print(io, "  Coordinates: [", join(coordinates(x), ", "), "]\n")
    print(io, "  Attitude: ", attitude(x))
end
