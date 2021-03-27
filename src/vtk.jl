function vtk_format(x::AbstractVector{<: Vec{dim, T}}) where {dim, T}
    n = length(x)
    v = reinterpret(T, Array(x))
    out = zeros(T, (dim == 2 ? 3 : dim), n)
    out[1:dim, :] .= reshape(v, dim, n)
    out
end

function WriteVTK.vtk_grid(vtk::AbstractString, line::Line)
    coords = vtk_format(line)
    cells = [MeshCell(VTKCellTypes.VTK_LINE, collect(1:size(coords, 2)))]
    vtk_grid(vtk, coords, cells)
end

function WriteVTK.vtk_grid(vtk::AbstractString, poly::Polygon)
    coords = vtk_format(poly)
    cells = [MeshCell(VTKCellTypes.VTK_POLYGON, collect(1:size(coords, 2)))]
    vtk_grid(vtk, coords, cells)
end

function WriteVTK.vtk_grid(vtk::AbstractString, circle::Circle)
    θᵢ = LinRange(0, 2π, 360)
    coords = vtk_format([centroid(circle) + radius(circle) * Vec(cos(θ), sin(θ)) for θ in θᵢ])
    cells = [MeshCell(VTKCellTypes.VTK_POLYGON, collect(1:size(coords, 2)))]
    vtk_grid(vtk, coords, cells)
end

function WriteVTK.vtk_grid(vtk::AbstractString, sphere::Sphere{2})
    vtk_grid(vtk, Circle(centroid(sphere), radius(sphere)))
end
