@inline Vec3(x::Vec{2, T}) where {T} = @inbounds Vec(x[1], x[2], zero(T))
@inline Vec2(x::Vec{3}) = @inbounds Vec(x[1], x[2])
