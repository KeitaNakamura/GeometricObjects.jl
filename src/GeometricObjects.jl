module GeometricObjects

using Reexport
@reexport using Tensorial
@reexport using WriteVTK
using RecipesBase

using Base: @_propagate_inbounds_meta

export
    GeometricObject,
    centroid,
    area,
    distance,
    moment_of_inertia,
    velocityat,
    translate!,
    rotate!,
# Line
    Line,
# Polygon
    Polygon,
    Rectangle,
# Sphere/Circle
    Sphere,
    Circle,
    radius,
# VTK
    vtk_grid

include("GeometricObject.jl")
include("Line.jl")
include("Polygon.jl")
include("Sphere.jl")
include("plots.jl")
include("vtk.jl")

end # module
