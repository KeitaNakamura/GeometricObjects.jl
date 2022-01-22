"""
    Line(a::Vec, b::Vec)
    Line(a::Vec => b::Vec)
"""
struct Line{dim, T} <: Shape{dim, T}
    coordinates::SVector{2, Vec{dim, T}}
    q::Quaternion{T}
end

Line(a::Vec, b::Vec) = Shape(Line, @SVector[a, b])
Line(pair::Pair) = Line(pair.first, pair.second)

centroid(line::Line) = sum(line) / 2

function moment_of_inertia(line::Line{2})
    a, b = line
    v = b - a
    l² = v ⋅ v
    symmetric(@Mat([0 0 0
                    0 0 0
                    0 0 l²/12]), :U)
end

function moment_of_inertia(line::Line{3, T}) where {T}
    a, b = line
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
    @inbounds begin
        v = line[2] - line[1]
        a_to_x = x - line[1]
    end
    scale = (a_to_x ⋅ v) / (v ⋅ v)
    distance = scale*v - a_to_x
    distance, scale
end

"""
    distance(::Line, x::Vec)
    distance(::Line, x::Vec, threshold::Real)

Compute the distance vector from `x` to perpendicular foot.
When `threshold` is given, check the contact between line and point `x`,
and return `nothing` if contact is not detected.
Note that if the perpendicular foot does not lie on the line,
contact detection is performed using distance between `x` and vertices of line.

```jldoctest
julia> line = Line(@Vec[0.0, 0.0] => @Vec[1.0, 1.0])
Line{2, Float64}:
  Coordinates: [[0.0, 0.0], [1.0, 1.0]]
  Attitude: [1.0, 0.0]

julia> distance(line, @Vec[1.0, 0.0])
2-element Vec{2, Float64}:
 -0.5
  0.5

julia> distance(line, @Vec[1.0, 0.0], 1.0)
2-element Vec{2, Float64}:
 -0.5
  0.5

julia> distance(line, @Vec[1.0, 0.0], 0.5) === nothing
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

Compute the position of perpendicular foot.
"""
function perpendicularfoot(line::Line, x::Vec)
    x + distance(line, x)
end

function normalunit(line::Line{2})
    @inbounds begin
        v = line[2] - line[1]
        n = Vec(v[2], -v[1])
    end
    normalize(n)
end

"""
    in(x::Vec, line::Line)

Check if `x` is `in` line.

# Examples
```jldoctest
julia> line = Line(@Vec[0.0, 0.0], @Vec[2.0, 2.0])
Line{2, Float64}:
  Coordinates: [[0.0, 0.0], [2.0, 2.0]]
  Attitude: [1.0, 0.0]

julia> @Vec[1.0, 1.0] in line
true

julia> @Vec[1.0, 0.0] in line
false
```
"""
@inline function Base.in(x::Vec, line::Line)
    @inbounds begin
        ax = line[1] - x
        bx = line[2] - x
    end
    ϵ = muladd(dot(ax,ax), dot(bx,bx), _pow2(dot(ax,bx)))
    abs2(ϵ) < eps(typeof(ϵ))
end
_pow2(x) = x * abs(x)

# helper function for `in(x, polygon)`
function ray_casting_to_right(line::Line{2}, X::Vec{2})
    @inbounds begin
        x, y = X[1], X[2]
        a, b = line
        a_x, a_y = a[1], a[2]
        b_x, b_y = b[1], b[2]
    end
    if (a_y ≤ y < b_y) || # upward case
       (b_y ≤ y < a_y)    # downward case
        x_line = a_x + (y - a_y) / (b_y - a_y) * (b_x - a_x)
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
    x1, y1 = line1[1]
    x2, y2 = line1[2]
    x3, y3 = line2[1]
    x4, y4 = line2[2]
    D = (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4)
    abs(D) < sqrt(eps(typeof(D))) && return nothing
    x1y2_y1x2 = x1*y2 - y1*x2
    x3y4_y3x4 = x3*y4 - y3*x4
    p1 = (x1y2_y1x2*(x3-x4) - (x1-x2)*x3y4_y3x4) / D
    p2 = (x1y2_y1x2*(y3-y4) - (y1-y2)*x3y4_y3x4) / D
    p = Vec(p1, p2)
    ifelse((extended[1] || p in line1) && (extended[2] || p in line2), p, nothing)
end

function Base.intersect(line1::Line{3}, line2::Line{3}; extended::Tuple{Bool, Bool} = (false, false))
    n1 = normalize(line1[2] - line1[1])
    n2 = normalize(line2[2] - line2[1])
    v = line2[1] - line1[1]
    n1_n2 = n1 ⋅ n2
    n1_v  = n1 ⋅ v
    n2_v  = n2 ⋅ v
    d1 = (n1_v - (n1_n2)*(n2_v)) / (1 - n1_n2 * n1_n2)
    d2 = ((n1_n2)*(n1_v) - n2_v) / (1 - n1_n2 * n1_n2)
    p1 = line1[1] + d1 * n1
    p2 = line2[1] + d2 * n2
    p1 ≈ p2 || return nothing
    ifelse((extended[1] || p1 in line1) && (extended[2] || p1 in line2), p1, nothing)
end
