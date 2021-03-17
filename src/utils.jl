@inline ToVec(::Val{3}, x::Vec{3}) = x
@inline ToVec(::Val{3}, x::Vec{2, T}) where {T} = @inbounds Vec(x[1], x[2], zero(T))

@inline ToVec(::Val{2}, x::Vec{3}) = @inbounds Vec(x[1], x[2])
@inline ToVec(::Val{2}, x::Vec{2}) = x
