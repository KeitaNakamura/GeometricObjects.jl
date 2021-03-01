mutable struct Sphere{dim, T} <: GeometricObject{dim, T}
    coordinates::Vector{Vec{dim, T}}
    r::T
    m::T
    I::Mat{3, 3, T}
    v::Vec{dim, T}
    θ::Vec{3, T}
    ω::Vec{3, T}
end

Sphere(center::Vec, r::Real) = GeometricObject(Sphere, [center], r)

center(x::Sphere) = @inbounds x.coordinates[1]
radius(x::Sphere) = x.r

function moment_of_inertia(x::Sphere)
    R = radius(x)
    I = 2R^2 / 5
    x.m * @Mat [I 0 0
                0 I 0
                0 0 I]
end

distance(line::Line, x::Sphere) = distance(line, center(x), radius(x))
distance(poly::Polygon, x::Sphere) = distance(poly, center(x), radius(x))


mutable struct Disk{dim, T} <: GeometricObject{dim, T}
    coordinates::Vector{Vec{dim, T}}
    r::T
    m::T
    I::Mat{3, 3, T}
    v::Vec{dim, T}
    θ::Vec{3, T}
    ω::Vec{3, T}
end

Disk(center::Vec{dim}, r::Real) where {dim} = GeometricObject(Disk, [center], r)

center(x::Disk) = @inbounds x.coordinates[1]
radius(x::Disk) = x.r

# http://hyperphysics.phy-astr.gsu.edu/hbase/tdisc.html
function moment_of_inertia(x::Disk{2})
    R = radius(x)
    x.m * @Mat [0 0 0
                0 0 0
                0 0 R^2/2]
end

distance(line::Line, x::Disk) = distance(line, center(x), radius(x))
distance(poly::Polygon, x::Disk) = distance(poly, center(x), radius(x))
