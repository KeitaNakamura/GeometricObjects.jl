function unzip(v::AbstractVector{<: Vec{dim, T}}) where {dim, T}
    xyz = ntuple(d -> Vector{T}(undef, length(v)), Val(dim))
    for i in eachindex(v)
        for d in 1:dim
            @inbounds xyz[d][i] = v[i][d]
        end
    end
    xyz
end

@recipe function f(line::Line)
    seriestype := :path
    aspect_ratio := :equal
    unzip(coordinates(line))
end

@recipe function f(poly::Polygon)
    seriestype := :geometry
    aspect_ratio := :equal
    unzip(coordinates(poly))
end

@recipe function f(sphere::Union{Sphere{2}, Circle{2}})
    seriestype := :geometry
    aspect_ratio := :equal
    xc = centroid(sphere)
    r = radius(sphere)
    θ -> r*cos(θ) + xc[1], θ -> r*sin(θ) + xc[2], 0, 2π
end
