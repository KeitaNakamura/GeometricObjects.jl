abstract type GeometricObject{dim, T} <: AbstractVector{Vec{dim, T}} end

function GeometricObject(::Type{Obj}, coordinates::Vector{Vec{dim, T}}, args...) where {Obj <: GeometricObject, dim, T}
    m = one(T)
    I = zero(Mat{3, 3, T})
    v = zero(Vec{dim, T})
    θ = zero(Vec{3, T})
    ω = zero(Vec{3, T})
    obj = Obj(coordinates, args..., m, I, v, θ, ω)
    obj.I = moment_of_inertia(obj)
    obj
end

moment_of_inertia(x::GeometricObject) = throw(ArgumentError("$(typeof(x)) is not supported yet."))

coordinates(x::GeometricObject) = x.coordinates

function velocityat(obj::GeometricObject{dim}, x::Vec{dim}) where {dim}
    r = x - center(obj)
    v = obj.v + obj.ω × r
    Vec{dim}(i -> @inbounds v[i])
end

Base.size(x::GeometricObject) = size(coordinates(x))

@inline function Base.getindex(x::GeometricObject, i::Int)
    @boundscheck checkbounds(x, i)
    @inbounds coordinates(x)[i]
end

@inline function Base.setindex!(x::GeometricObject, v, i::Int)
    @boundscheck checkbounds(x, i)
    @inbounds coordinates(x)[i] = v
    x
end
