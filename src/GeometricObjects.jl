module GeometricObjects

using Reexport
@reexport using Tensorial
@reexport using WriteVTK
using StaticArrays
using RecipesBase

using Base: @_propagate_inbounds_meta

export
# Shape
    Shape,
    coordinates,
    centroid,
    area,
    distance,
    moment_of_inertia,
    translate!,
    rotate!,
    enlarge,
# Line
    Line,
# Polygon
    Polygon,
    Rectangle,
# Polyline
    Polyline,
# Sphere/Circle
    Sphere,
    Circle,
    radius,
# GeometricObject,
    GeometricObject,
    velocityat,
# VTK
    vtk_grid

include("Shape.jl")
include("Line.jl")
include("Polygon.jl")
include("Polyline.jl")
include("Sphere.jl")
include("GeometricObject.jl")
include("plots.jl")
include("vtk.jl")

end # module
