mutable struct Polyline{dim, T, V <: AbstractVector{Vec{dim, T}}} <: Shape{dim, T}
    coordinates::V
    q::Quaternion{T}
end

Polyline(coordinates::AbstractVector{<: Vec}) = Shape(Polyline, coordinates)

@inline function getline(poly::Polyline, i::Int)
    @boundscheck @assert 1 ≤ i ≤ length(poly)-1
    @inbounds SLine(poly[i], poly[i+1])
end

function distance(poly::Polyline{dim, T}, x::Vec{dim, T}, r::T; include_tips::Bool = true) where {dim, T}
    dist = zero(Vec{dim, T})
    count = 0
    @inbounds for i in 1:length(poly)-1
        line = getline(poly, i)
        d = distance(line, x, r)
        if d !== nothing
            dist += d
            count += 1
        end
    end
    count != 0 && return dist / count

    @inbounds for i in 1:length(poly)
        if !include_tips
            (i == 1 || i == length(poly)) && continue
        end
        xᵢ = poly[i]
        xᵢ in SCircle(x, r) && return xᵢ - x
    end
    nothing
end

# almost copied from `Polygon`, should be unified?
"""
    intersect(::Polyline, ::AbstractLine; [extended = false])

Find the closest intersection point from line to polyline.
Return `nothing` if not found.
"""
function Base.intersect(poly::Polyline{dim, T}, line::AbstractLine; extended::Bool = false) where {dim, T}
    output = zero(Vec{dim, T})
    dist = T(Inf)
    @inbounds for i in 1:length(poly)-1
        p = intersect(getline(poly, i), line; extended = (false, extended))
        p === nothing && continue
        v = line[1] - p
        vv = v ⋅ v
        if vv < dist
            output = p
            dist = vv
        end
    end
    isinf(dist) && return nothing
    output
end
