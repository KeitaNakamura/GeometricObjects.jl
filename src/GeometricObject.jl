abstract type GeometricObject{dim, T, S} end

mutable struct GeometricObject2D{T, S <: Geometry{2, T}} <: GeometricObject{2, T, S}
    geometry::S
    m::T
    v::Vec{2, T}
    ω::T
end

mutable struct GeometricObject3D{T, S <: Geometry{3, T}} <: GeometricObject{3, T, S}
    geometry::S
    m::T
    v::Vec{3, T}
    ω::Vec{3, T}
end

function GeometricObject(geometry::Geometry{2, T}; m = one(T), v = zero(Vec{2, T}), ω = zero(T)) where {T}
    GeometricObject2D(geometry, m, v, ω)
end
function GeometricObject(geometry::Geometry{3, T}; m = one(T), v = zero(Vec{3, T}), ω = zero(Vec{3, T})) where {T}
    GeometricObject3D(geometry, m, v, ω)
end

geometry(x::GeometricObject) = x.geometry
coordinates(x::GeometricObject) = coordinates(geometry(x))
coordinates(x::GeometricObject, i::Int) = (@_propagate_inbounds_meta; coordinates(geometry(x), i))
num_coordinates(x::GeometricObject) = num_coordinates(geometry(x))
quaternion(x::GeometricObject) = quaternion(geometry(x))
moment_of_inertia(x::GeometricObject) = x.m * moment_of_inertia(geometry(x))

function inv_moment_of_inertia(I::SymmetricSecondOrderTensor{3, T}) where {T}
    V, P = eigen(I)
    A⁻¹ = SymmetricSecondOrderTensor{3}((i,j) -> @inbounds i == j ? (abs(V[i]) < sqrt(eps(T)) ? V[i] : inv(V[i])) : zero(T))
    symmetric(P ⋅ A⁻¹ ⋅ inv(P), :U)
end
inv_moment_of_inertia(I::Real) = inv(I)

"""
    translate!(object::GeometricObject, u::Vec)

Translate `object` by the displacement `u`.
"""
function translate!(obj::GeometricObject, u::Vec)
    translate!(geometry(obj), u)
    obj
end

"""
    rotate!(object::GeometricObject, θ::Vec)

Rotate `object` by the angle vector `θ`.
`normalize(θ)` and `norm(θ)` should represent the rotation axis and the angle (radian), respectively.
"""
function rotate!(obj::GeometricObject, θ::Union{Vec, Real})
    rotate!(geometry(obj), θ)
    obj
end

"""
    update_geometry!(::GeometricObject, Δt)

Update geometry of object by timestep `Δt`.
Current linear and angular velocities of object are used in the calculation.
"""
function update_geometry!(obj::GeometricObject, Δt::Real)
    # v and ω need to be updated in advance
    Δx = obj.v * Δt
    Δθ = obj.ω * Δt
    translate!(obj, Δx)
    rotate!(obj, Δθ)
    obj
end

"""
    apply_force!(object::GeometricObject, F, τ, Δt)

Apply linear force `F` and torque (moment of force) `τ` to `object` with timestep `Δt`.
This only updates the linear velocity `object.v` and the angular velocity `object.ω`.
See also [`update_geometry!`](@ref).
"""
function apply_force!(obj::GeometricObject{dim}, F::Vec{dim}, τ::Union{Vec{dim}, Real}, Δt::Real) where {dim}
    isinf(obj.m) && return obj
    m = obj.m
    v = obj.v
    ω = obj.ω
    I = moment_of_inertia(obj)
    I⁻¹ = inv_moment_of_inertia(I)

    # update velocities
    obj.v = v + F / m * Δt
    obj.ω = ω + I⁻¹ ⋅ (τ - _cross(ω, I⋅ω)) * Δt

    obj
end
_cross(x, y) = cross(x, y)
_cross(x::Real, y::Real) = zero(promote_type(typeof(x), typeof(y))) # for 2D case

Broadcast.broadcastable(geo::GeometricObject) = (geo,)

function Base.show(io::IO, mime::MIME"text/plain", x::GeometricObject)
    print(io, typeof(x), ":\n")
    buf = IOBuffer()
    show(buf, MIME("text/plain"), geometry(x))
    strings = map(line -> "  " * line, eachline(IOBuffer(take!(buf))))
    print(io, join(strings, "\n"), "\n")
    print(io, "  Mass: ", x.m, "\n")
    print(io, "  Velocity: ", x.v, "\n")
    print(io, "  Angular velocity: ", x.ω)
end
