struct Polygon{dim, T, L} <: Geometry{dim, T}
    coordinates::SVector{L, Vec{dim, T}}
    q::Quaternion{T}
end

_projection(v::Vec{2}, ::Nothing, ::Nothing, ::Nothing) = v
_projection(v::Vec{2}, x::Real, ::Nothing, ::Nothing) = Vec(x, v[1], v[2]) # (y,z)
_projection(v::Vec{2}, ::Nothing, y::Real, ::Nothing) = Vec(v[2], y, v[1]) # (z,x)
_projection(v::Vec{2}, ::Nothing, ::Nothing, z::Real) = Vec(v[1], v[2], z) # (x,y)

function Polygon(coordinates::Vec{2}...; x = nothing, y = nothing, z = nothing)
    coords = _projection.(coordinates, x, y, z)
    Geometry(Polygon, SVector(coords))
end

function Rectangle(bottomleft::Vec{2, T}, topright::Vec{2, T}; x = nothing, y = nothing, z = nothing) where {T}
    x0 = bottomleft[1]
    y0 = bottomleft[2]
    x1 = topright[1]
    y1 = topright[2]
    Polygon(Vec(x0, y0), Vec(x1, y0), Vec(x1, y1), Vec(x0, y1); x, y, z)
end

@inline repeated(x::AbstractVector, i::Int) = @inbounds x[(i-1) % length(x) + 1]

# https://en.wikipedia.org/wiki/Centroid
function centroid(poly::Polygon{2, T}) where {T}
    A = zero(T)
    x_c = zero(T)
    y_c = zero(T)
    C = coordinates(poly)
    @simd for i in 1:num_coordinates(poly)
        @inbounds begin
            Xᵢ = repeated(C, i)
            Xᵢ₊₁ = repeated(C, i+1)
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
    C = coordinates(poly)
    @simd for i in 1:num_coordinates(poly)
        @inbounds begin
            Xᵢ = repeated(C, i)
            Xᵢ₊₁ = repeated(C, i+1)
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
    xc = centroid(poly)
    num = zero(T)
    den = zero(T)
    C = coordinates(poly)
    @simd for i in 1:num_coordinates(poly)
        @inbounds begin
            xᵢ = repeated(C, i) - xc
            xᵢ₊₁ = repeated(C, i+1) - xc
        end
        a = norm(xᵢ₊₁ × xᵢ)
        num += a * ((xᵢ ⋅ xᵢ) + (xᵢ ⋅ xᵢ₊₁) + (xᵢ₊₁ ⋅ xᵢ₊₁))
        den += a
    end
    symmetric(@Mat([0 0 0
                    0 0 0
                    0 0 num/6den]), :U)
end

@inline function getline(poly::Polygon, i::Int)
    C = coordinates(poly)
    i == length(C) && return @inbounds Line(C[i], C[1])
    @boundscheck checkbounds(C, i)
    @inbounds Line(C[i], C[i+1])
end

function Base.eachline(poly::Polygon)
    (@inbounds(getline(poly, i)) for i in 1:num_coordinates(poly))
end

"""
    in(x::Vec, ::Polygon; include_bounds = true)

Check if `x` is `in` a polygon.
"""
@inline function Base.in(x::Vec{2}, poly::Polygon{2}; include_bounds::Bool = true)
    I = 0
    @inbounds @simd for i in 1:num_coordinates(poly)
        line = getline(poly, i)
        I += (x in line)
    end
    J = 0
    @inbounds @simd for i in 1:num_coordinates(poly)
        line = getline(poly, i)
        J += ray_casting_to_right(line, x)
    end
    ifelse(I === 1 || I === 2, include_bounds, isodd(J))
end

function distance(poly::Polygon{2, T}, x::Vec{2, T}, r::T) where {T}
    dist = zero(Vec{2, T})
    count = 0
    @inbounds for i in 1:num_coordinates(poly)
        line = getline(poly, i)
        d = distance(line, x, r)
        if d !== nothing
            dist += d
            count += 1
        end
    end
    count != 0 && return dist / count

    for xᵢ in coordinates(poly)
        xᵢ in Circle(x, r) && return xᵢ - x
    end
    nothing
end

function distance(poly::Polygon{2, T}, x::Vec{2, T}, r::T, line_values::AbstractVector{U}) where {T, U}
    @assert num_coordinates(poly) == length(line_values)
    dist = zero(Vec{2, T})
    value = zero(U)
    count = 0
    @inbounds for i in 1:num_coordinates(poly)
        line = getline(poly, i)
        d = distance(line, x, r)
        if d !== nothing
            dist += d
            value += line_values[i]
            count += 1
        end
    end
    count != 0 && return (dist/count, value/count)

    @inbounds for i in 1:num_coordinates(poly)
        xᵢ = coordinates(poly, i)
        if xᵢ in Circle(x, r)
            if i == 1
                return (xᵢ - x), (line_values[end] + line_values[i]) / 2
            else
                return (xᵢ - x), (line_values[i-1] + line_values[i]) / 2
            end
        end
    end
    nothing
end

"""
    intersect(::Polygon, ::Line; [extended = false])

Find the closest intersection point from line to polygon.
Return `nothing` if not found.
"""
function Base.intersect(poly::Polygon{dim, T}, line::Line; extended::Bool = false) where {dim, T}
    output = zero(Vec{dim, T})
    dist = T(Inf)
    @inbounds for i in 1:num_coordinates(poly)
        p = intersect(getline(poly, i), line; extended = (false, extended))
        p === nothing && continue
        v = coordinates(line, 1) - p
        vv = v ⋅ v
        if vv < dist
            output = p
            dist = vv
        end
    end
    isinf(dist) && return nothing
    output
end
