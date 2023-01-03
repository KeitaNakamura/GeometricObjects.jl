struct Polygon{dim, T, L} <: Geometry{dim, T}
    coordinates::SVector{L, Vec{dim, T}}
    q::Quaternion{T}
end

function Polygon(coordinates::Vec{2}...)
    Geometry(Polygon, SVector(coordinates))
end

function Rectangle(bottomleft::Vec{2}, topright::Vec{2})
    x0 = bottomleft[1]
    y0 = bottomleft[2]
    x1 = topright[1]
    y1 = topright[2]
    Polygon(Vec(x0, y0), Vec(x1, y0), Vec(x1, y1), Vec(x0, y1))
end

# https://en.wikipedia.org/wiki/Centroid
function centroid(poly::Polygon{2, T}) where {T}
    A = zero(T)
    x_c = zero(T)
    y_c = zero(T)
    C = SVector(coordinates(poly)..., coordinates(poly, 1))
    @simd for i in 1:num_coordinates(poly)
        @inbounds begin
            Xᵢ = C[i]
            Xᵢ₊₁ = C[i+1]
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
    C = SVector(coordinates(poly)..., coordinates(poly, 1))
    @simd for i in 1:num_coordinates(poly)
        @inbounds begin
            Xᵢ = C[i]
            Xᵢ₊₁ = C[i+1]
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
    C = SVector(coordinates(poly)..., coordinates(poly, 1))
    @simd for i in 1:num_coordinates(poly)
        @inbounds begin
            xᵢ = C[i] - xc
            xᵢ₊₁ = C[i+1] - xc
        end
        a = norm(xᵢ₊₁ × xᵢ)
        num += a * ((xᵢ ⋅ xᵢ) + (xᵢ ⋅ xᵢ₊₁) + (xᵢ₊₁ ⋅ xᵢ₊₁))
        den += a
    end
    num/6den
end

@inline function getline(poly::Polygon, i::Int)
    @_propagate_inbounds_meta
    C = SVector(coordinates(poly)..., coordinates(poly, 1))
    Line(C[i], C[i+1])
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

############
# distance #
############

function distance(poly::Polygon{2}, x::Vec{2}, r::Real)
    d = distance_from_sides((a...)->nothing, poly, x, r)
    d !== nothing && return d
    index = distance_from_vertices(poly, x, r)
    index !== nothing && return coordinates(poly, index) - x
    nothing
end

function distance(poly::Polygon{2}, x::Vec{2}, r::Real, ::Symbol)
    inds = Vector{Int}(undef, num_coordinates(poly))
    n = Ref(0)
    value = distance_from_sides(poly, x, r) do count, i
        n[] = count
        inds[count] = i
    end
    value !== nothing && return (value, resize!(inds, n[]))
    index = distance_from_vertices(poly, x, r)
    index === nothing && return nothing
    resize!(inds, 2)
    @inbounds begin
        inds[1] = ifelse(index==1, num_coordinates(poly), index-1)
        inds[2] = index
    end
    coordinates(poly, index)-x, inds
end

# helpers
function distance_from_sides(f, poly::Polygon{2, T}, x::Vec{2, T}, r::T) where {T}
    val = zero(Vec{2, T})
    count = 0
    @inbounds @simd for i in 1:num_coordinates(poly)
        line = getline(poly, i)
        d = distance(line, x, r)
        if d !== nothing
            val += d
            f(count+=1, i)
        end
    end
    count==0 ? nothing : val/count
end
function distance_from_vertices(poly::Polygon{2, T}, x::Vec{2, T}, r::T) where {T}
    index = 0
    ddmin = T(Inf)
    @inbounds @simd for i in 1:num_coordinates(poly)
        xᵢ = coordinates(poly, i)
        if xᵢ in Circle(x, r)
            d = xᵢ - x
            dd = d ⋅ d
            if dd < ddmin
                ddmin = dd
                index = i
            end
        end
    end
    index==0 ? nothing : index
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
