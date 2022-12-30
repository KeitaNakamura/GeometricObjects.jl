mutable struct GeometricObject{dim, T, S <: Geometry{dim, T}}
    geometry::S
    m::T
    v::Vec{dim, T}
    ω::Vec{3, T}
end

function GeometricObject(geometry::Geometry{dim, T}; m = one(T), v = zero(Vec{dim, T}), ω = zero(Vec{3, T})) where {dim, T}
    GeometricObject(geometry, m, v, ω)
end

geometry(x::GeometricObject) = x.geometry
coordinates(x::GeometricObject) = coordinates(geometry(x))
coordinates(x::GeometricObject, i::Int) = (@_propagate_inbounds_meta; coordinates(geometry(x), i))
quaternion(x::GeometricObject) = quaternion(geometry(x))
attitude(x::GeometricObject) = attitude(geometry(x))

function inv_moment_of_inertia(I::SymmetricSecondOrderTensor{3, T}, ::Val{3}) where {T}
    V, P = eigen(I)
    A⁻¹ = SymmetricSecondOrderTensor{3}((i,j) -> @inbounds i == j ? (abs(V[i]) < sqrt(eps(T)) ? V[i] : inv(V[i])) : zero(T))
    symmetric(P ⋅ A⁻¹ ⋅ inv(P), :U)
end

function inv_moment_of_inertia(I::SymmetricSecondOrderTensor{3, T}, ::Val{2}) where {T}
    symmetric(@Mat([0 0 0
                    0 0 0
                    0 0 inv(I[3,3])]), :U)
end

function translate!(obj::GeometricObject, u::Vec)
    obj.geometry = translate(geometry(obj), u)
    obj
end

function rotate!(obj::GeometricObject, θ::Union{Vec, Real})
    obj.geometry = rotate(geometry(obj), θ)
    obj
end

"""
    update_geometry!(::GeometricObject, dt)

Update geometry of object by time increment `dt`.
Current linear and angular velocities of object are used in the calculation.
"""
function update_geometry! end

function update_geometry!(obj::GeometricObject{3}, dt::Real)
    # v and ω need to be updated in advance
    dx = obj.v * dt
    dθ = obj.ω * dt
    translate!(obj, dx)
    rotate!(obj, dθ)
    obj
end

function update_geometry!(obj::GeometricObject{2}, dt::Real)
    # v and ω need to be updated in advance
    dx = obj.v * dt
    dθ = obj.ω * dt
    translate!(obj, dx)
    rotate!(obj, dθ[3])
    obj
end

# 3D
function apply_force!(obj::GeometricObject{3}, F::Vec{3}, τ::Vec{3}, dt::Real)
    m = obj.m
    v = obj.v
    ω = obj.ω
    I = moment_of_inertia(geometry(obj))
    I⁻¹ = inv_moment_of_inertia(I, Val(3))

    # update velocities first
    obj.v = v + F / m * dt
    obj.ω = ω + I⁻¹ ⋅ (τ - ω × (I ⋅ ω)) * dt

    update_geometry!(obj, dt)
    obj
end

# 2D
function apply_force!(obj::GeometricObject{2}, F::Vec{2}, τ::Real, dt::Real)
    m = obj.m
    v = obj.v
    ω = obj.ω[3]
    I = moment_of_inertia(geometry(obj))[3,3]

    # update velocities first
    obj.v = v + F / m * dt
    obj.ω = Vec(0, 0, ω + inv(I)*τ*dt)

    update_geometry!(obj, dt)
    obj
end

Base.in(x::Vec, obj::GeometricObject) = x in geometry(obj)

for IterType in (:Tuple, :AbstractArray)
    @eval function Base.findall(pred::Base.Fix2{typeof(in), <: GeometricObject}, iter::$IterType)
        findall(x -> in(x, pred.x), iter)
    end
end

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
