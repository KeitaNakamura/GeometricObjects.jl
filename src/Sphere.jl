abstract type AbstractCircle{dim, T} <: Shape{dim, T} end

function area(x::AbstractCircle{2})
    r = radius(x)
    π*r^2
end

# http://hyperphysics.phy-astr.gsu.edu/hbase/tdisc.html
function moment_of_inertia(x::AbstractCircle{2})
    r = radius(x)
    symmetric(@Mat([0 0 0
                    0 0 0
                    0 0 r^2/2]), :U)
end

# call methods for SSphere
Base.in(x::Vec{2}, circle::AbstractCircle; include_bounds::Bool = true) = in(x, SSphere(circle); include_bounds)
distance(circle::AbstractCircle, x::Vec{2}) = distance(SSphere(circle), x)
distance(circle::AbstractCircle, x::Vec{2}, r::Real) = distance(SSphere(circle), x, r)

# static version of circle
struct SCircle{dim, T} <: AbstractCircle{dim, T}
    xc::Vec{dim, T}
    r::T
    reverse::Bool
end

SCircle(centroid::Vec, r::Real) = SCircle(centroid, r, false)

coordinates(x::SCircle) = @SVector[centroid(x)]
centroid(x::SCircle) = x.xc
radius(x::SCircle) = x.r

enlarge(circle::SCircle, R::Real) = SCircle(centroid(circle), R*radius(circle), circle.reverse)
Base.reverse(circle::SCircle) = SCircle(centroid(circle), radius(circle), !circle.reverse)

mutable struct Circle{dim, T} <: AbstractCircle{dim, T}
    coordinates::MVector{1, Vec{dim, T}}
    r::T
    reverse::Bool
    q::Quaternion{T}
end

Circle(centroid::Vec{dim}, r::Real) where {dim} = Shape(Circle, @MVector[centroid], r, false)

centroid(x::Circle) = @inbounds x[1]
radius(x::Circle) = x.r

enlarge(circle::Circle, R::Real) = Circle(copy(coordinates(circle)), R*radius(circle), circle.reverse, circle.q)
Base.reverse(circle::Circle) = Circle(copy(coordinates(circle)), radius(circle), !circle.reverse, circle.q)


abstract type AbstractSphere{dim, T} <: Shape{dim, T} end

function moment_of_inertia(x::AbstractSphere)
    r = radius(x)
    I = 2r^2 / 5
    symmetric(@Mat([I 0 0
                    0 I 0
                    0 0 I]), :U)
end

"""
    in(x::Vec, ::AbstractSphere; include_bounds = true)

Check if `x` is `in` a sphere.

```jldoctest
julia> sphere = Sphere(Vec(1.0,1.0), 1.0)
1-element Sphere{2,Float64}:
 [1.0, 1.0]

julia> Vec(0.5, 0.5) in sphere
true

julia> Vec(0.0, 0.0) in sphere
false
```
"""
function Base.in(x::Vec{dim}, sphere::AbstractSphere{dim}; include_bounds::Bool = true) where {dim}
    c = centroid(sphere)
    d = (x - c)
    d² = d ⋅ d
    r² = radius(sphere)^2
    d² == r² && return include_bounds
    d² < r²
end

"""
    distance(::AbstractSphere, x::Vec)
    distance(::AbstractSphere, x::Vec, threshold::Real)

Compute the distance vector from `x` to perpendicular foot on surface of sphere.
When `threshold` is given, check the contact between line and point `x`,
and return `nothing` if contact is not detected.

```jldoctest
julia> sphere = Sphere(Vec(1.0,1.0), 1.0)
1-element Sphere{2,Float64}:
 [1.0, 1.0]

julia> d = distance(sphere, Vec(0.0,0.0))
2-element Tensor{Tuple{2},Float64,1,2}:
 0.29289321881345254
 0.29289321881345254

julia> norm(d) ≈ √2 - 1
true
```
"""
function distance(sphere::AbstractSphere{dim}, x::Vec{dim}) where {dim}
    v = centroid(sphere) - x
    norm_v = norm(v)
    d = norm_v - radius(sphere)
    d * (v / norm_v)
end

function distance(sphere::AbstractSphere{dim}, x::Vec{dim}, r::Real) where {dim}
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

struct SSphere{dim, T} <: AbstractSphere{dim, T}
    xc::Vec{dim, T}
    r::T
    reverse::Bool
end

SSphere(centroid::Vec, r::Real) = SSphere(centroid, r, false)
SSphere(circle::AbstractCircle) = SSphere(centroid(circle), radius(circle), circle.reverse)

coordinates(x::SSphere) = @SVector[centroid(x)]

centroid(x::SSphere) = x.xc
radius(x::SSphere) = x.r

enlarge(sphere::SSphere, R::Real) = SSphere(centroid(sphere), R*radius(sphere), sphere.reverse)
Base.reverse(sphere::SSphere) = SSphere(centroid(sphere), radius(sphere), !sphere.reverse)

mutable struct Sphere{dim, T} <: AbstractSphere{dim, T}
    coordinates::MVector{1, Vec{dim, T}}
    r::T
    reverse::Bool
    q::Quaternion{T}
end

Sphere(centroid::Vec, r::Real) = Shape(Sphere, @MVector[centroid], r, false)
# cannot construct Sphere from SCircle
Sphere(circle::Circle) = Sphere(coordinates(circle), radius(circle), circle.reverse, circle.q)

centroid(x::Sphere) = @inbounds x[1]
radius(x::Sphere) = x.r

enlarge(sphere::Sphere, R::Real) = Sphere(copy(coordinates(sphere)), R*radius(sphere), sphere.reverse, sphere.q)
Base.reverse(sphere::Sphere) = Sphere(copy(coordinates(sphere)), radius(sphere), !sphere.reverse, sphere.q)
