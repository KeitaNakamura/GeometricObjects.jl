abstract type Shape{dim, T} end

function Shape(::Type{S}, coordinates::AbstractVector{<: Vec{dim, T}}, args...) where {S <: Shape, dim, T}
    q = quaternion(T, 0, Vec(0,0,1))
    S(coordinates, args..., q)
end

coordinates(x::Shape) = x.coordinates
centered(x::Shape) = coordinates(x) .- centroid(x) # call coordinates to keep type of coordinates in broadcast, otherwise always return `Vector`

Base.size(x::Shape) = (length(x),)
Base.length(x::Shape) = length(coordinates(x))

Base.checkbounds(x::Shape, i::Int) = checkbounds(coordinates(x), i)
Base.checkbounds(::Type{Bool}, x::Shape, i::Int) = checkbounds(Bool, coordinates(x), i)

@inline function Base.getindex(x::Shape, i::Int)
    @boundscheck checkbounds(x, i)
    @inbounds coordinates(x)[i]
end

@inline function Base.setindex!(x::Shape, v, i::Int)
    @boundscheck checkbounds(x, i)
    @inbounds coordinates(x)[i] = v
    x
end

Base.:(==)(x::Shape, y::Shape) = coordinates(x) == coordinates(y)

moment_of_inertia(x::Shape) = throw(ArgumentError("$(typeof(x)) is not supported yet."))

@generated function enlarge(shape::Shape, R::Real)
    exps = [name == :coordinates ? :coordinates : :(shape.$name) for name in fieldnames(shape)]
    quote
        coordinates = centroid(shape) .+ R * centered(shape)
        $shape($(exps...))
    end
end

@generated function Base.reverse(shape::Shape)
    exps = [name == :coordinates ? :coordinates : :(shape.$name) for name in fieldnames(shape)]
    quote
        coordinates = reverse(shape.coordinates)
        $shape($(exps...))
    end
end

function translate!(shape::Shape, u::Vec)
    xᵢ = coordinates(shape)
    @. xᵢ = xᵢ + u
    shape
end

function rotate!(shape::Shape, θ::Vec)
    # https://www.ashwinnarayan.com/post/how-to-integrate-quaternions/
    q = exp(Quaternion(θ/2))
    xc = centroid(shape)
    xᵢ = coordinates(shape)
    @. xᵢ = xc + rotate(xᵢ - xc, q)
    shape.q = q * shape.q
    shape
end
