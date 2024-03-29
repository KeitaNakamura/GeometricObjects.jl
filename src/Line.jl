"""
    Line(a::Vec, b::Vec)
    Line(a::Vec => b::Vec)
"""
mutable struct Line{dim, T} <: Geometry{dim, T}
    coordinates::MVector{2, Vec{dim, T}}
    q::Quaternion{T}
end

Line(a::Vec, b::Vec) = Geometry(Line, MVector(a,b))
Line(pair::Pair) = Line(pair.first, pair.second)

centroid(line::Line) = mean(coordinates(line))
norm(line::Line) = (c=coordinates(line); norm(c[1] - c[2]))

function moment_of_inertia(line::Line{2})
    a, b = coordinates(line)
    v = b - a
    l² = v ⋅ v
    l²/12
end

function moment_of_inertia(line::Line{3, T}) where {T}
    a, b = coordinates(line)
    v = b - a
    l² = v ⋅ v
    I = symmetric(@Mat([l²/12 0     0
                        0     l²/12 0
                        0     0     0]), :U)
    zaxis = Vec{3, T}(0,0,1)
    R = rotmat(zaxis => v/sqrt(l²))
    rotate(I, R)
end

@inline function _distance(line::Line, x::Vec)
    c1, c2 = coordinates(line)
    @inbounds begin
        v = c2 - c1
        a_to_x = x - c1
    end
    scale = (a_to_x ⋅ v) / (v ⋅ v)
    distance = scale*v - a_to_x
    distance, scale
end

"""
    distance(::Line, x::Vec)
    distance(::Line, x::Vec, threshold::Real)

Return the distance vector from `x` to perpendicular foot.
When `threshold` is given, check the contact between line and point `x`,
and return `nothing` if contact is not detected.

```jldoctest
julia> line = Line(Vec(0.0, 0.0) => Vec(1.0, 1.0))
Line{2, Float64}:
  Centroid: [0.5, 0.5]
  Quaternion: 1.0 + 0.0𝙞 + 0.0𝙟 + 0.0𝙠

julia> distance(line, Vec(1.0, 0.0))
2-element Vec{2, Float64}:
 -0.5
  0.5

julia> distance(line, Vec(1.0, 0.0), 1.0)
2-element Vec{2, Float64}:
 -0.5
  0.5

julia> distance(line, Vec(1.0, 0.0), 0.5) === nothing
true
```
"""
function distance(line::Line, x::Vec)
    _distance(line, x)[1]
end

@inline function distance(line::Line, x::Vec, r::Real)
    r² = r^2
    d, scale = _distance(line, x)
    if 0 ≤ scale ≤ 1 # perpendicular foot is on line
        (d ⋅ d) ≤ r² && return d
    end
    nothing
end

"""
    GeometricObjects.perpendicularfoot(::Line, x::Vec)

Return the position of perpendicular foot.
"""
function perpendicularfoot(line::Line, x::Vec)
    x + distance(line, x)
end

function normalunit(line::Line{2})
    c1, c2 = coordinates(line)
    @inbounds begin
        v = c2 - c1
        n = Vec(v[2], -v[1])
    end
    normalize(n)
end

"""
    ison(x::Vec, line::Line)

Check if a point `ison` line.

# Examples
```jldoctest
julia> line = Line(Vec(0.0, 0.0), Vec(2.0, 2.0))
Line{2, Float64}:
  Centroid: [1.0, 1.0]
  Quaternion: 1.0 + 0.0𝙞 + 0.0𝙟 + 0.0𝙠

julia> ison(Vec(1.0, 1.0), line)
true

julia> ison(Vec(1.0, 0.0), line)
false
```
"""
@inline function ison(x::Vec, line::Line)
    c1, c2 = coordinates(line)
    @inbounds begin
        ax = c1 - x
        bx = c2 - x
    end
    ϵ = muladd(dot(ax,ax), dot(bx,bx), _pow2(dot(ax,bx)))
    abs2(ϵ) < eps(typeof(ϵ))^2
end
_pow2(x) = x * abs(x)

# helper function for `ison(x, polygon)`
function ray_casting_to_right(line::Line{2}, X::Vec{2})
    @inbounds begin
        x, y = X[1], X[2]
        a, b = coordinates(line)
        a_x, a_y = a[1], a[2]
        b_x, b_y = b[1], b[2]
    end
    if (a_y ≤ y < b_y) || (b_y ≤ y < a_y) # upward case || downward case
        x_line = muladd((y - a_y) / (b_y - a_y), b_x - a_x, a_x)
        x < x_line && return true
    end
    false
end

"""
    intersect(::Line, ::Line; [extended = (false, false)])

Find intersection point from two lines.
Return `nothing` if not found.
"""
function Base.intersect(line1::Line{2}, line2::Line{2}; extended::Tuple{Bool, Bool} = (false, false))
    x1, y1 = coordinates(line1, 1)
    x2, y2 = coordinates(line1, 2)
    x3, y3 = coordinates(line2, 1)
    x4, y4 = coordinates(line2, 2)
    D = (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4)
    abs(D) < sqrt(eps(typeof(D))) && return nothing
    x1y2_y1x2 = x1*y2 - y1*x2
    x3y4_y3x4 = x3*y4 - y3*x4
    p1 = (x1y2_y1x2*(x3-x4) - (x1-x2)*x3y4_y3x4) / D
    p2 = (x1y2_y1x2*(y3-y4) - (y1-y2)*x3y4_y3x4) / D
    p = Vec(p1, p2)
    ifelse((extended[1] || ison(p, line1)) && (extended[2] || ison(p, line2)), p, nothing)
end

function Base.intersect(line1::Line{3}, line2::Line{3}; extended::Tuple{Bool, Bool} = (false, false))
    c1 = coordinates(line1)
    c2 = coordinates(line2)
    n1 = normalize(c1[2] - c1[1])
    n2 = normalize(c2[2] - c2[1])
    v = c2[1] - c1[1]
    n1_n2 = n1 ⋅ n2
    n1_v  = n1 ⋅ v
    n2_v  = n2 ⋅ v
    d1 = (n1_v - (n1_n2)*(n2_v)) / (1 - n1_n2 * n1_n2)
    d2 = ((n1_n2)*(n1_v) - n2_v) / (1 - n1_n2 * n1_n2)
    p1 = c1[1] + d1 * n1
    p2 = c2[1] + d2 * n2
    p1 ≈ p2 || return nothing
    ifelse((extended[1] || ison(p1, line1)) && (extended[2] || ison(p1, line2)), p1, nothing)
end
