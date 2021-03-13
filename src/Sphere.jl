mutable struct Sphere{dim, T} <: GeometricObject{dim, T}
    coordinates::Vector{Vec{dim, T}}
    r::T
    m::T
    v::Vec{dim, T}
    ω::Vec{3, T}
    q::Quaternion{T}
end

Sphere(center::Vec, r::Real) = GeometricObject(Sphere, [center], r)

center(x::Sphere) = @inbounds x[1]
radius(x::Sphere) = x.r

function moment_of_inertia(x::Sphere)
    r = radius(x)
    I = 2r^2 / 5
    x.m * symmetric(@Mat([I 0 0
                          0 I 0
                          0 0 I]), :U)
end

distance(line::Line, x::Sphere) = distance(line, center(x), radius(x))
distance(poly::Polygon, x::Sphere) = distance(poly, center(x), radius(x))


mutable struct Disk{dim, T} <: GeometricObject{dim, T}
    coordinates::Vector{Vec{dim, T}}
    r::T
    m::T
    v::Vec{dim, T}
    ω::Vec{3, T}
    q::Quaternion{T}
end

Disk(center::Vec{dim}, r::Real) where {dim} = GeometricObject(Disk, [center], r)

center(x::Disk) = @inbounds x[1]
radius(x::Disk) = x.r

# http://hyperphysics.phy-astr.gsu.edu/hbase/tdisc.html
function moment_of_inertia(x::Disk{2})
    r = radius(x)
    x.m * symmetric(@Mat([0 0 0
                          0 0 0
                          0 0 r^2/2]), :U)
end

distance(line::Line, x::Disk) = distance(line, center(x), radius(x))
distance(poly::Polygon, x::Disk) = distance(poly, center(x), radius(x))
