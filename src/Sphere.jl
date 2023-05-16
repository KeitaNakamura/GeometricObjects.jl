mutable struct Circle{dim, T} <: Geometry{dim, T}
    coordinates::MVector{1, Vec{dim, T}}
    q::Quaternion{T}
    r::T
end

Circle(centroid::Vec{dim}, r::Real) where {dim} = Geometry(Circle, @MVector[centroid], r)

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
distance(circle::Circle, x::Vec{2}, r::Real; inverse::Bool=false) = distance(Sphere(circle), x, r; inverse)


mutable struct Sphere{dim, T} <: Geometry{dim, T}
    coordinates::MVector{1, Vec{dim, T}}
    q::Quaternion{T}
    r::T
end

Sphere(centroid::Vec, r::Real) = Geometry(Sphere, @MVector[centroid], r)
Sphere(circle::Circle) = Sphere(coordinates(circle), circle.q, radius(circle))

centroid(x::Sphere) = @inbounds only(coordinates(x))
radius(x::Sphere) = x.r

function volume(x::Sphere{3})
    r = radius(x)
    (4π*r^3)/3
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
  Centroid: [1.0, 1.0, 1.0]
  Quaternion: 1.0 + 0.0𝙞 + 0.0𝙟 + 0.0𝙠

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
julia> sphere = Sphere(Vec(1.0,1.0,1.0), 1.0)
Sphere{3, Float64}:
  Centroid: [1.0, 1.0, 1.0]
  Quaternion: 1.0 + 0.0𝙞 + 0.0𝙟 + 0.0𝙠

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

function distance(sphere::Sphere{dim}, x::Vec{dim}, r::Real; inverse::Bool=false) where {dim}
    if inverse
        v = x - centroid(sphere)
        norm_v = norm(v)
        d = radius(sphere) - norm_v
    else
        v = centroid(sphere) - x
        norm_v = norm(v)
        d = norm_v - radius(sphere)
    end
    d ≤ r && return d * (v / norm_v)
    nothing
end
