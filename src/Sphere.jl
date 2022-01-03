struct Circle{dim, T} <: Shape{dim, T}
    coordinates::SVector{1, Vec{dim, T}}
    q::Quaternion{T}
    r::T
    reverse::Bool
end

Circle(centroid::Vec{dim}, r::Real) where {dim} = Shape(Circle, @SVector[centroid], r, false)

enlarge(circle::Circle, R::Real) = Circle(coordinates(circle), quaternion(circle), R*radius(circle), circle.reverse)
Base.reverse(circle::Circle) = Circle(coordinates(circle), quaternion(circle), radius(circle), !circle.reverse)

centroid(x::Circle) = @inbounds x[1]
radius(x::Circle) = x.r

function area(x::Circle{2})
    r = radius(x)
    π*r^2
end

# http://hyperphysics.phy-astr.gsu.edu/hbase/tdisc.html
function moment_of_inertia(x::Circle{2})
    r = radius(x)
    symmetric(@Mat([0 0 0
                    0 0 0
                    0 0 r^2/2]), :U)
end

# call methods for Sphere
Base.in(x::Vec{2}, circle::Circle; include_bounds::Bool = true) = in(x, Sphere(circle); include_bounds)
distance(circle::Circle, x::Vec{2}) = distance(Sphere(circle), x)
distance(circle::Circle, x::Vec{2}, r::Real) = distance(Sphere(circle), x, r)


struct Sphere{dim, T} <: Shape{dim, T}
    coordinates::SVector{1, Vec{dim, T}}
    q::Quaternion{T}
    r::T
    reverse::Bool
end

Sphere(centroid::Vec, r::Real) = Shape(Sphere, @SVector[centroid], r, false)
# cannot construct Sphere from Circle
Sphere(circle::Circle) = Sphere(coordinates(circle), circle.q, radius(circle), circle.reverse)

centroid(x::Sphere) = @inbounds x[1]
radius(x::Sphere) = x.r

enlarge(sphere::Sphere, R::Real) = Sphere(coordinates(sphere), sphere.q, R*radius(sphere), sphere.reverse)
Base.reverse(sphere::Sphere) = Sphere(coordinates(sphere), sphere.q, radius(sphere), !sphere.reverse)

function moment_of_inertia(x::Sphere)
    r = radius(x)
    I = 2r^2 / 5
    symmetric(@Mat([I 0 0
                    0 I 0
                    0 0 I]), :U)
end

"""
    in(x::Vec, ::Sphere; include_bounds = true)

Check if `x` is `in` a sphere.

```jldoctest
julia> sphere = Sphere(Vec(1.0,1.0), 1.0)
Sphere{2, Float64}:
  Coordinates: [[1.0, 1.0]]
  Attitude: [1.0, 0.0]

julia> Vec(0.5, 0.5) in sphere
true

julia> Vec(0.0, 0.0) in sphere
false
```
"""
function Base.in(x::Vec{dim}, sphere::Sphere{dim}; include_bounds::Bool = true) where {dim}
    c = centroid(sphere)
    d = x - c
    d² = d ⋅ d
    r² = radius(sphere)^2
    d² == r² ? include_bounds : d² < r²
end

"""
    distance(::Sphere, x::Vec)
    distance(::Sphere, x::Vec, threshold::Real)

Compute the distance vector from `x` to perpendicular foot on surface of sphere.
When `threshold` is given, check the contact between line and point `x`,
and return `nothing` if contact is not detected.

```jldoctest
julia> sphere = Sphere(Vec(1.0,1.0), 1.0)
Sphere{2, Float64}:
  Coordinates: [[1.0, 1.0]]
  Attitude: [1.0, 0.0]

julia> d = distance(sphere, Vec(0.0,0.0))
2-element Vec{2, Float64}:
 0.29289321881345254
 0.29289321881345254

julia> norm(d) ≈ √2 - 1
true
```
"""
function distance(sphere::Sphere{dim}, x::Vec{dim}) where {dim}
    v = centroid(sphere) - x
    norm_v = norm(v)
    d = norm_v - radius(sphere)
    d * (v / norm_v)
end

function distance(sphere::Sphere{dim}, x::Vec{dim}, r::Real) where {dim}
    if sphere.reverse
        v = x - centroid(sphere)
        norm_v = norm(v)
        d = radius(sphere) - norm_v
        d - r ≤ 0 && return d * (v / norm_v)
        nothing
    else
        v = centroid(sphere) - x
        norm_v = norm(v)
        d = norm_v - radius(sphere)
        d - r ≤ 0 && return d * (v / norm_v)
        nothing
    end
end
