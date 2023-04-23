function vtk_format(x::AbstractVector{Vec{dim, T}}) where {dim, T}
    reinterpret(reshape, T, x)
end
function vtk_format(x::SVector{<: Any, Vec{dim, T}}) where {dim, T}
    reshape(vcat(x...), (dim, length(x)))
end

function WriteVTK.vtk_grid(vtk::AbstractString, line::Line; kwargs...)
    coords = vtk_format(coordinates(line))
    cells = [MeshCell(VTKCellTypes.VTK_LINE, 1:size(coords, 2))]
    vtk_grid(vtk, coords, cells; kwargs...)
end

function WriteVTK.vtk_grid(vtk::AbstractString, poly::Polygon; kwargs...)
    coords = vtk_format(coordinates(poly))
    cells = [MeshCell(VTKCellTypes.VTK_POLYGON, 1:size(coords, 2))]
    vtk_grid(vtk, coords, cells; kwargs...)
end

function WriteVTK.vtk_grid(vtk::AbstractString, circle::Circle; kwargs...)
    θᵢ = LinRange(0, 2π, 360)
    coords = vtk_format([centroid(circle) + radius(circle) * Vec(cos(θ), sin(θ)) for θ in θᵢ])
    cells = [MeshCell(VTKCellTypes.VTK_POLYGON, 1:size(coords, 2))]
    vtk_grid(vtk, coords, cells; kwargs...)
end

function WriteVTK.vtk_grid(vtk::AbstractString, sphere::Sphere{2}; kwargs...)
    vtk_grid(vtk, Circle(centroid(sphere), radius(sphere)); kwargs...)
end

WriteVTK.vtk_grid(vtk::AbstractString, obj::GeometricObject; kwargs...) = vtk_grid(vtk, geometry(obj); kwargs...)
