struct Circle{dim, T} <: Geometry{dim, T}
    coordinates::SVector{1, Vec{dim, T}}
    q::Quaternion{T}
    r::T
end

Circle(centroid::Vec{dim}, r::Real) where {dim} = Geometry(Circle, @SVector[centroid], r)

enlarge(circle::Circle, R::Real) = Circle(coordinates(circle), quaternion(circle), R*radius(circle))

centroid(x::Circle) = @inbounds only(coordinates(x))
radius(x::Circle) = x.r

function area(x::Circle{2})
    r = radius(x)
    π*r^2
end

# http://hyperphysics.phy-astr.gsu.edu/hbase/tdisc.html
function moment_of_inertia(x::Circle{2})
    r = radius(x)
    r^2/2
end

# call methods for Sphere
isinside(x::Vec{2}, circle::Circle; include_bounds::Bool = true) = isinside(x, Sphere(circle); include_bounds)
distance(circle::Circle, x::Vec{2}) = distance(Sphere(circle), x)
distance(circle::Circle, x::Vec{2}, r::Real) = distance(Sphere(circle), x, r)


struct Sphere{dim, T} <: Geometry{dim, T}
    coordinates::SVector{1, Vec{dim, T}}
    q::Quaternion{T}
    r::T
end

Sphere(centroid::Vec, r::Real) = Geometry(Sphere, @SVector[centroid], r)
Sphere(circle::Circle) = Sphere(coordinates(circle), circle.q, radius(circle))

centroid(x::Sphere) = @inbounds only(coordinates(x))
radius(x::Sphere) = x.r
enlarge(sphere::Sphere, R::Real) = Sphere(coordinates(sphere), sphere.q, R*radius(sphere))

function volume(x::Sphere{3})
    r = radius(x)
    (4π*r^2)/3
end

function moment_of_inertia(x::Sphere)
    r = radius(x)
    I = 2r^2 / 5
    symmetric(@Mat([I 0 0
                    0 I 0
                    0 0 I]), :U)
end

"""
    isinside(x::Vec, ::Sphere; include_bounds = true)

Check if a point `isinside` a sphere.

```jldoctest
julia> sphere = Sphere(Vec(1.0,1.0,1.0), 1.0)
Sphere{3, Float64}:
  Coordinates: [[1.0, 1.0, 1.0]]
  Attitude: [1.0, 0.0, 0.0]

julia> isinside(Vec(0.5, 0.5, 0.5), sphere)
true

julia> isinside(Vec(0.0, 0.0, 0.0), sphere)
false
```
"""
function isinside(x::Vec{dim}, sphere::Sphere{dim}; include_bounds::Bool = true) where {dim}
    c = centroid(sphere)
    d = x - c
    d² = d ⋅ d
    r² = radius(sphere)^2
    d² == r² ? include_bounds : d² < r²
end

"""
    distance(::Sphere, x::Vec)
    distance(::Sphere, x::Vec, threshold::Real)

Return the distance vector from `x` to perpendicular foot on surface of sphere.
When `threshold` is given, check the contact between line and point `x`,
and return `nothing` if contact is not detected.

```jldoctest
julia> sphere = Sphere(Vec(1.0,1.0), 1.0)
Sphere{2, Float64}:
  Coordinates: [[1.0, 1.0]]
  Attitude: [1.0, 0.0]

julia> d = distance(sphere, Vec(0.0,0.0,0.0))
3-element Vec{3, Float64}:
 0.4226497308103742
 0.4226497308103742
 0.4226497308103742

julia> norm(d) ≈ √3 - 1
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
    v = centroid(sphere) - x
    norm_v = norm(v)
    d = norm_v - radius(sphere)
    abs(d) ≤ r && return d * (v / norm_v)
    nothing
end
