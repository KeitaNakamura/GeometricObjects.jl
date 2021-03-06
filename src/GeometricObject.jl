mutable struct GeometricObject{dim, T, S <: Shape{dim, T}} <: AbstractVector{Vec{dim, T}}
    shape::S
    m::T
    v::Vec{dim, T}
    ω::Vec{3, T}
end

function GeometricObject(shape::Shape{dim, T}) where {dim, T}
    m = one(T)
    v = zero(Vec{dim, T})
    ω = zero(Vec{3, T})
    GeometricObject(shape, m, v, ω)
end

coordinates(x::GeometricObject) = coordinates(x.shape)

for f in (:centroid, :centered, :area, :translate!, :rotate!, :distance) # call the same function of `Shape`
    @eval $f(x::GeometricObject, args...; kwargs...) = $f(x.shape, args...; kwargs...)
end

Base.size(x::GeometricObject) = size(coordinates(x))

@inline function Base.getindex(x::GeometricObject, i::Int)
    @boundscheck checkbounds(x, i)
    @inbounds coordinates(x)[i]
end

@inline function Base.setindex!(x::GeometricObject, v, i::Int)
    @boundscheck checkbounds(x, i)
    @inbounds coordinates(x)[i] = v
    x
end

moment_of_inertia(x::GeometricObject) = x.m * moment_of_inertia(x.shape)

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

function update!(obj::GeometricObject, dt::Real)
    # v and ω need to be updated in advance
    dx = obj.v * dt
    dθ = obj.ω * dt
    translate!(obj.shape, dx)
    rotate!(obj.shape, dθ)
    obj
end

function update!(obj::GeometricObject{3}, F::Vec{3}, τ::Vec{3}, dt::Real)
    m = obj.m
    v = obj.v
    ω = obj.ω
    I = moment_of_inertia(obj)
    I⁻¹ = inv_moment_of_inertia(I, Val(3))

    # update two velocities first
    obj.v = v + F / m * dt
    obj.ω = ω + I⁻¹ ⋅ (τ - ω × (I ⋅ ω)) * dt

    update!(obj, dt)
    obj
end

function update!(obj::GeometricObject{2}, F::Vec{2}, τ::Vec{3}, dt::Real)
    m = obj.m
    v = obj.v
    ω = obj.ω
    I = moment_of_inertia(obj)

    # update two velocities first
    obj.v = v + F / m * dt
    @inbounds obj.ω = Vec(0, 0, ω[3] + inv(I[3,3]) * τ[3] * dt)

    update!(obj, dt)
    obj
end

function compute_force_moment(obj::GeometricObject, Fᵢ::AbstractArray{<: Vec}, xᵢ::AbstractArray{<: Vec})
    promote_shape(Fᵢ, xᵢ)
    xc = centroid(obj)
    F = sum(Fᵢ)
    M = sum((x - xc) × F for (F, x) in zip(Fᵢ, xᵢ))
    F, M
end

function update!(obj::GeometricObject{dim, T}, Fᵢ::AbstractArray{<: Vec}, xᵢ::AbstractArray{<: Vec}, dt::Real; body_force_per_unit_mass::Vec = zero(Vec{dim, T})) where {dim, T}
    F, M = compute_force_moment(obj, Fᵢ, xᵢ)
    update!(obj, F + obj.m*body_force_per_unit_mass, M, dt)
    obj
end

Base.in(x::Vec, obj::GeometricObject) = x in obj.shape
function Base.in(obj::GeometricObject, x::Vec)
    throw(ArgumentError("`in(obj, x)` is invalid, use `in(x, obj)` instead"))
end
