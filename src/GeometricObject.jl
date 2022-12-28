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

function velocityat(obj::GeometricObject{3}, x::Vec{3})
    r = x - centroid(obj)
    obj.v + obj.ω × r
end

function velocityat(obj::GeometricObject{2}, x::Vec{2})
    r = x - centroid(obj)
    @inbounds begin
        v3 = obj.ω × Vec(r[1], r[2], 0)
        obj.v + Vec(v3[1], v3[2])
    end
end

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

function update!(obj::GeometricObject{3}, dt::Real)
    # v and ω need to be updated in advance
    dx = obj.v * dt
    dθ = obj.ω * dt
    translate!(obj, dx)
    rotate!(obj, dθ)
    obj
end

function update!(obj::GeometricObject{2}, dt::Real)
    # v and ω need to be updated in advance
    dx = obj.v * dt
    dθ = obj.ω * dt
    translate!(obj, dx)
    rotate!(obj, dθ[3])
    obj
end

function update!(obj::GeometricObject{3}, F::Vec{3}, τ::Vec{3}, dt::Real)
    m = obj.m
    v = obj.v
    ω = obj.ω
    I = moment_of_inertia(geometry(obj))
    I⁻¹ = inv_moment_of_inertia(I, Val(3))

    # update two velocities first
    obj.v = v + F / m * dt
    obj.ω = ω + I⁻¹ ⋅ (τ - ω × (I ⋅ ω)) * dt

    update!(obj, dt)
    obj
end

# 2D
function update!(obj::GeometricObject{2}, F::Vec{2}, τ::Vec{3}, dt::Real)
    m = obj.m
    v = obj.v
    ω = obj.ω
    I = moment_of_inertia(geometry(obj))

    # update two velocities first
    obj.v = v + F / m * dt
    @inbounds obj.ω = Vec(0, 0, ω[3] + inv(I[3,3]) * τ[3] * dt)

    update!(obj, dt)
    obj
end

function compute_force_moment(obj::GeometricObject{dim, T}, Fᵢ::AbstractArray{Vec{dim, T}}, xᵢ::AbstractArray{Vec{dim, T}}) where {dim, T}
    promote_shape(Fᵢ, xᵢ)
    xc = centroid(geometry(obj))
    F = sum(Fᵢ)
    M = sum((x - xc) × F for (F, x) in zip(Fᵢ, xᵢ); init = zero(Vec{3, T}))
    F, M
end

function update!(obj::GeometricObject{dim, T}, Fᵢ::AbstractArray{<: Vec}, xᵢ::AbstractArray{<: Vec}, dt::Real; body_force_per_unit_mass::Vec = zero(Vec{dim, T})) where {dim, T}
    F, M = compute_force_moment(geometry(obj), Fᵢ, xᵢ)
    update!(obj, F + obj.m*body_force_per_unit_mass, M, dt)
    obj
end

Base.in(x::Vec, obj::GeometricObject) = x in coordinates(obj)
function Base.in(obj::GeometricObject, x::Vec)
    throw(ArgumentError("`in(obj, x)` is invalid, use `in(x, obj)` instead"))
end

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
