function vtk_format(x::AbstractVector{Vec{dim, T}}) where {dim, T}
    reshape(reinterpret(T, x), (dim, length(x)))
end

function WriteVTK.vtk_grid(vtk::AbstractString, line::Line)
    coords = vtk_format(coordinates(line))
    cells = [MeshCell(VTKCellTypes.VTK_LINE, collect(1:size(coords, 2)))]
    vtk_grid(vtk, coords, cells)
end

function WriteVTK.vtk_grid(vtk::AbstractString, poly::Polygon)
    coords = vtk_format(coordinates(poly))
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

WriteVTK.vtk_grid(vtk::AbstractString, obj::GeometricObject) = vtk_grid(vtk, geometry(obj))
