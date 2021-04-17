function unzip(v::AbstractVector{<: Vec{dim, T}}) where {dim, T}
    xyz = ntuple(d -> Vector{T}(undef, length(v)), Val(dim))
    for i in eachindex(v)
        for d in 1:dim
            @inbounds xyz[d][i] = v[i][d]
        end
    end
    xyz
end

@recipe function f(line::AbstractLine)
    seriestype := :path
    aspect_ratio := :equal
    unzip(line)
end

@recipe function f(poly::Polygon)
    seriestype := :shape
    aspect_ratio := :equal
    unzip(poly)
end

@recipe function f(sphere::Union{AbstractSphere{2}, AbstractCircle{2}})
    seriestype := :shape
    aspect_ratio := :equal
    xc = centroid(sphere)
    r = radius(sphere)
    θ -> r*cos(θ) + xc[1], θ -> r*sin(θ) + xc[2], 0, 2π
end
