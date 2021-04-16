mutable struct Polygon{dim, T, V <: AbstractVector{Vec{dim, T}}} <: Shape{dim, T}
    coordinates::V
    q::Quaternion{T}
end

_projection(v::Vec{2}, ::Nothing, ::Nothing, ::Nothing) = v
_projection(v::Vec{2}, x::Real, ::Nothing, ::Nothing) = Vec(x, v[1], v[2]) # (y,z)
_projection(v::Vec{2}, ::Nothing, y::Real, ::Nothing) = Vec(v[2], y, v[1]) # (z,x)
_projection(v::Vec{2}, ::Nothing, ::Nothing, z::Real) = Vec(v[1], v[2], z) # (x,y)

function Polygon(coordinates::Vector{<: Vec{2}}; x = nothing, y = nothing, z = nothing)
    coords = _projection.(coordinates, x, y, z)
    Shape(Polygon, coords)
end

function Polygon(coordinates::Vec{2}...; x = nothing, y = nothing, z = nothing)
    coords = _projection.(coordinates, x, y, z)
    Shape(Polygon, MVector(coords))
end

function Rectangle(bottomleft::Vec{2, T}, topright::Vec{2, T}; x = nothing, y = nothing, z = nothing) where {T}
    x0 = bottomleft[1]
    y0 = bottomleft[2]
    x1 = topright[1]
    y1 = topright[2]
    Polygon(@MVector[Vec(x0, y0), Vec(x1, y0), Vec(x1, y1), Vec(x0, y1)]; x, y, z)
end

# handle end+1 index
@inline function Base.getindex(poly::Polygon, i::Int)
    if i == length(poly) + 1
        @inbounds coordinates(poly)[1]
    else
        @boundscheck checkbounds(poly, i)
        @inbounds coordinates(poly)[i]
    end
end

# https://en.wikipedia.org/wiki/Centroid
function centroid(poly::Polygon{2, T}) where {T}
    A = zero(T)
    x_c = zero(T)
    y_c = zero(T)
    for i in 1:length(poly)
        @inbounds begin
            Xᵢ = poly[i]
            Xᵢ₊₁ = poly[i+1]
            xᵢ, yᵢ = Xᵢ[1], Xᵢ[2]
            xᵢ₊₁, yᵢ₊₁ = Xᵢ₊₁[1], Xᵢ₊₁[2]
        end
        a = (xᵢ * yᵢ₊₁ - xᵢ₊₁ * yᵢ)
        x_c += (xᵢ + xᵢ₊₁) * a
        y_c += (yᵢ + yᵢ₊₁) * a
        A += a
    end
    A /= 2
    Vec(x_c/6A, y_c/6A)
end

function area(poly::Polygon{2, T}) where {T}
    A = zero(T)
    for i in 1:length(poly)
        @inbounds begin
            Xᵢ = poly[i]
            Xᵢ₊₁ = poly[i+1]
            xᵢ, yᵢ = Xᵢ[1], Xᵢ[2]
            xᵢ₊₁, yᵢ₊₁ = Xᵢ₊₁[1], Xᵢ₊₁[2]
        end
        a = (xᵢ * yᵢ₊₁ - xᵢ₊₁ * yᵢ)
        A += a
    end
    A /= 2
    A
end

# https://en.wikipedia.org/wiki/List_of_moments_of_inertia
function moment_of_inertia(poly::Polygon{2, T}) where {T}
    num = zero(T)
    den = zero(T)
    for i in 1:length(poly)
        @inbounds begin
            xᵢ = poly[i]
            xᵢ₊₁ = poly[i+1]
        end
        a = norm(xᵢ₊₁ × xᵢ)
        num += a * ((xᵢ ⋅ xᵢ) + (xᵢ ⋅ xᵢ₊₁) + (xᵢ₊₁ ⋅ xᵢ₊₁))
        den += a
    end
    symmetric(@Mat([0 0 0
                    0 0 0
                    0 0 num/den]), :U)
end

@inline function getline(poly::Polygon, i::Int)
    @boundscheck checkbounds(poly, i)
    if i == length(poly)
        @inbounds SLine(poly[i], poly[1])
    else
        @inbounds SLine(poly[i], poly[i+1])
    end
end

function Base.eachline(poly::Polygon)
    (@inbounds(getline(poly, i)) for i in eachindex(poly))
end

"""
    in(x::Vec, ::Polygon; include_bounds = true)

Check if `x` is `in` a polygon.
"""
function Base.in(x::Vec{2}, poly::Polygon{2}; include_bounds::Bool = true)
    I = 0
    for line in eachline(poly)
        x in line && return include_bounds
        I += ray_casting_to_right(line, x)
    end
    isodd(I)
end

function distance(poly::Polygon{2, T}, x::Vec{2, T}, r::Real) where {T}
    dist = nothing
    norm_dist = T(Inf)
    isincontact = false
    for (i, line) in enumerate(eachline(poly))
        d = distance_from_outside(line, x, r)
        if d !== nothing
            norm_d = norm(d)
            if norm_d < norm_dist
                dist = d
                norm_dist = norm_d
            end
        end
    end
    dist
end
