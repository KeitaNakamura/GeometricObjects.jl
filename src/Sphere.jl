mutable struct Sphere{dim, T} <: GeometricObject{dim, T}
    coordinates::Vector{Vec{dim, T}}
    r::T
    m::T
    v::Vec{dim, T}
    ω::Vec{3, T}
    q::Quaternion{T}
end

Sphere(centroid::Vec, r::Real) = GeometricObject(Sphere, [centroid], r)

centroid(x::Sphere) = @inbounds x[1]
radius(x::Sphere) = x.r

function moment_of_inertia(x::Sphere)
    r = radius(x)
    I = 2r^2 / 5
    x.m * symmetric(@Mat([I 0 0
                          0 I 0
                          0 0 I]), :U)
end

"""
    within(x::Vec, ::Sphere; include_bounds = true)

Check if `x` is within a sphere.

```jldoctest
julia> sphere = Sphere(Vec(1.0,1.0), 1.0)
1-element Sphere{2,Float64}:
 [1.0, 1.0]

julia> within(Vec(0.5, 0.5), sphere)
true

julia> within(Vec(0.0, 0.0), sphere)
false
```
"""
function within(x::Vec{dim}, sphere::Sphere{dim}; include_bounds::Bool = true) where {dim}
    c = centroid(sphere)
    d = (x - c)
    d² = d ⋅ d
    r² = radius(sphere)^2
    d² == r² && return include_bounds
    d² < r²
end

"""
    distance(::Sphere, x::Vec)
    distance(::Sphere, x::Vec, threshold::Real)

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
    d - r ≤ 0 && return (d - r) * (v / norm_v)
    nothing
end


mutable struct Disk{dim, T} <: GeometricObject{dim, T}
    coordinates::Vector{Vec{dim, T}}
    r::T
    m::T
    v::Vec{dim, T}
    ω::Vec{3, T}
    q::Quaternion{T}
end

Disk(centroid::Vec{dim}, r::Real) where {dim} = GeometricObject(Disk, [centroid], r)

centroid(x::Disk) = @inbounds x[1]
radius(x::Disk) = x.r

# http://hyperphysics.phy-astr.gsu.edu/hbase/tdisc.html
function moment_of_inertia(x::Disk{2})
    r = radius(x)
    x.m * symmetric(@Mat([0 0 0
                          0 0 0
                          0 0 r^2/2]), :U)
end

to_sphere(disk::Disk) = Sphere(coordinates(disk), radius(disk))
within(disk::Disk, x::Vec{2}; include_bounds::Bool = true) = within(to_sphere(disk), x; include_bounds)
distance(disk::Disk, x::Vec{2}) = distance(to_sphere(disk), x)
distance(disk::Disk, x::Vec{2}, r::Real) = distance(to_sphere(disk), x, r)
