abstract type Shape{dim, T} end

function Shape(::Type{S}, coordinates::AbstractVector{<: Vec{dim, T}}, args...) where {S <: Shape, dim, T}
    q = one(Quaternion{T})
    S(coordinates, q, args...)
end

coordinates(x::Shape) = x.coordinates
quaternion(x::Shape) = x.q
attitude(x::Shape{2, T}) where {T} = rotate(Vec{2,T}(1,0), quaternion(x))
attitude(x::Shape{3, T}) where {T} = rotate(Vec{3,T}(1,0,0), quaternion(x))

Base.size(x::Shape) = size(coordinates(x))
Base.length(x::Shape) = length(coordinates(x))

Base.eltype(x::Shape{dim, T}) where {dim, T} = Vec{dim, T}
Base.eachindex(x::Shape) = Base.OneTo(length(x))
Base.firstindex(x::Shape) = 1
Base.lastindex(x::Shape) = length(x)

Base.checkbounds(x::Shape, i...) = checkbounds(coordinates(x), i...)
@inline function Base.getindex(x::Shape, i::Int)
    @boundscheck checkbounds(x, i)
    @inbounds coordinates(x)[i]
end

# @inline Base.iterate(x::Shape, i = 1) = (i % UInt) - 1 < length(x) ? (@inbounds x[i], i + 1) : nothing

moment_of_inertia(x::Shape) = throw(ArgumentError("$(typeof(x)) is not supported yet."))

@generated function copy_shape(shape::Shape, coordinates::SVector, q::Quaternion)
    exps = [name in (:coordinates, :q) ? name : :(shape.$name) for name in fieldnames(shape)]
    quote
        $shape($(exps...))
    end
end

function centered(x::Shape)
    copy_shape(x, coordinates(x) .- centroid(x), quaternion(x))
end

function enlarge(shape::Shape, R::Real)
    copy_shape(
        shape,
        centroid(shape) .+ R * coordinates(centered(shape)),
        quaternion(shape),
    )
end

function translate(shape::Shape, u::Vec)
    copy_shape(shape, coordinates(shape) .+ u, quaternion(shape))
end

rotate(shape::Shape{3}, θ::Vec{3}) = _rotate(shape, θ)
rotate(shape::Shape{2}, θ::Real) = _rotate(shape, Vec(0,0,θ))
function _rotate(shape::Shape{dim}, θ::Vec) where {dim}
    # https://www.ashwinnarayan.com/post/how-to-integrate-quaternions/
    q = exp(Quaternion(θ/2))
    xc = centroid(shape)
    copy_shape(
        shape,
        (@. xc + $Tensorial.resizedim(rotate($coordinates(shape) - xc, q), Val(dim))),
        q * quaternion(shape),
    )
end

for IterType in (:Tuple, :AbstractArray)
    @eval function Base.findall(pred::Base.Fix2{typeof(in), <: Shape}, iter::$IterType)
        findall(x -> in(x, pred.x), iter)
    end
end

# isapprox
Base.isapprox(x::Shape, y::AbstractVector; kwargs...) = isapprox(coordinates(x), y; kwargs...)
Base.isapprox(x::AbstractVector, y::Shape; kwargs...) = isapprox(x, coordinates(y); kwargs...)

function Base.show(io::IO, mime::MIME"text/plain", x::Shape)
    print(io, typeof(x), ":\n")
    print(io, "  Coordinates: [", join(coordinates(x), ", "), "]\n")
    print(io, "  Attitude: ", attitude(x))
end
