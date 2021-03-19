abstract type GeometricObject{dim, T} <: AbstractVector{Vec{dim, T}} end

function GeometricObject(::Type{Obj}, coordinates::AbstractVector{<: Vec{dim, T}}, args...) where {Obj <: GeometricObject, dim, T}
    m = one(T)
    v = zero(Vec{dim, T})
    ω = zero(Vec{3, T})
    q = quaternion(T, 0, Vec(0,0,1))
    Obj(coordinates, args..., m, v, ω, q)
end

coordinates(x::GeometricObject) = x.coordinates

centered(x::GeometricObject) = x .- centroid(x)

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

moment_of_inertia(x::GeometricObject) = throw(ArgumentError("$(typeof(x)) is not supported yet."))

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

function update_position!(obj::GeometricObject{dim}, dt::Real) where {dim}
    ω = obj.ω
    q = obj.q

    # v and ω need to be updated in advance
    # https://www.ashwinnarayan.com/post/how-to-integrate-quaternions/
    xc = centroid(obj) + obj.v * dt
    dq = exp(Quaternion(ω*dt/2))
    r = centered(obj)

    @inbounds for i in eachindex(obj)
        obj[i] = xc + ToVec(Val(dim), (dq * r[i] / dq).vector)
    end
    obj.q = dq * q

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

    update_position!(obj, dt)
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

    update_position!(obj, dt)
    obj
end

function translate!(obj::GeometricObject, u::Vec)
    for i in eachindex(obj)
        @inbounds obj[i] = obj[i] + u
    end
    obj
end

function rotate!(obj::GeometricObject, R::Mat)
    xc = centroid(obj)
    for i in eachindex(obj)
        @inbounds obj[i] = rotate(obj[i] - xc, R) + xc
    end
    obj
end
