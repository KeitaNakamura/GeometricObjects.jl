abstract type AbstractLine{dim, T} <: Shape{dim, T} end

centroid(line::AbstractLine) = sum(line) / 2

function moment_of_inertia(line::AbstractLine{2})
    a, b = line
    v = b - a
    l² = v ⋅ v
    symmetric(@Mat([0 0 0
                    0 0 0
                    0 0 l²/12]), :U)
end

function moment_of_inertia(line::AbstractLine{3, T}) where {T}
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

function _distance(line::AbstractLine, x::Vec)
    @inbounds begin
        v = line[2] - line[1]
        a_to_x = x - line[1]
    end
    scale = (a_to_x ⋅ v) / (v ⋅ v)
    distance = scale*v - a_to_x
    distance, scale
end

"""
    distance(::AbstractLine, x::Vec)
    distance(::AbstractLine, x::Vec, threshold::Real)

Compute the distance vector from `x` to perpendicular foot.
When `threshold` is given, check the contact between line and point `x`,
and return `nothing` if contact is not detected.
Note that if the perpendicular foot does not lie on the line,
contact detection is performed using distance between `x` and vertices of line.

```jldoctest
julia> line = Line(@Vec[0.0, 0.0] => @Vec[1.0, 1.0])
2-element Line{2,Float64}:
 [0.0, 0.0]
 [1.0, 1.0]

julia> distance(line, @Vec[1.0, 0.0])
2-element Tensor{Tuple{2},Float64,1,2}:
 -0.5
  0.5

julia> distance(line, @Vec[1.0, 0.0], 1.0)
2-element Tensor{Tuple{2},Float64,1,2}:
 -0.5
  0.5

julia> distance(line, @Vec[1.0, 0.0], 0.5) === nothing
true
```
"""
function distance(line::AbstractLine, x::Vec)
    _distance(line, x)[1]
end

function distance(line::AbstractLine, x::Vec, r::Real)
    r² = r^2
    d, scale = _distance(line, x)
    if 0 ≤ scale ≤ 1 # perpendicular foot is on line
        (d ⋅ d) ≤ r² && return d
    end
    nothing
end

"""
    GeometricObjects.perpendicularfoot(::AbstractLine, x::Vec)

Compute the position of perpendicular foot.
"""
function perpendicularfoot(line::AbstractLine, x::Vec)
    x + distance(line, x)
end

function normalunit(line::AbstractLine{2})
    @inbounds begin
        v = line[2] - line[1]
        n = Vec(v[2], -v[1])
    end
    normalize(n)
end

"""
    in(x::Vec, line::AbstractLine)

Check if `x` is `in` line.

# Examples
```jldoctest
julia> line = Line(@Vec[0.0, 0.0], @Vec[2.0, 2.0])
2-element Line{2,Float64}:
 [0.0, 0.0]
 [2.0, 2.0]

julia> @Vec[1.0, 1.0] in line
true

julia> @Vec[1.0, 0.0] in line
false
```
"""
function Base.in(x::Vec, line::AbstractLine)
    d, scale = _distance(line, x)
    0 ≤ scale ≤ 1 && (d ⋅ d) < eps(eltype(d))
end
# much faster for 2D
function Base.in(X::Vec{2}, line::AbstractLine{2})
    @inbounds begin
        x, y = X[1], X[2]
        a, b = line
        a_x, a_y = a[1], a[2]
        b_x, b_y = b[1], b[2]
    end
    if (a_x ≤ x ≤ b_x) || (b_x ≤ x ≤ a_x)
        if (a_y ≤ y ≤ b_y) || (b_y ≤ y ≤ a_y)
            (x - a_x) * (b_y - a_y) == (y - a_y) * (b_x - a_x) && return true
        end
    end
    false
end

# helper function for `in(x, polygon)`
function ray_casting_to_right(line::AbstractLine{2}, X::Vec{2})
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

# static version for fast computation
struct SLine{dim, T} <: AbstractLine{dim, T}
    a::Vec{dim, T}
    b::Vec{dim, T}
end

coordinates(line::SLine) = SVector(line.a, line.b)

function Base.getindex(line::SLine, i::Int)
    @boundscheck checkbounds(line, i)
    i == 1 && return line.a
    i == 2 && return line.b
    error() # unreachable
end

"""
    Line(a::Vec, b::Vec)
    Line(a::Vec => b::Vec)
"""
mutable struct Line{dim, T} <: AbstractLine{dim, T}
    coordinates::MVector{2, Vec{dim, T}}
    q::Quaternion{T}
end

Line(a::Vec, b::Vec) = Shape(Line, @MVector[a, b])
Line(pair::Pair) = Line(pair.first, pair.second)
