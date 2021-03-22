module GeometricObjects

using Reexport
@reexport using Tensorial
using RecipesBase

using Base: @_propagate_inbounds_meta

export
    GeometricObject,
    centroid,
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
    radius

include("GeometricObject.jl")
include("Line.jl")
include("Polygon.jl")
include("Sphere.jl")
include("plots.jl")

end # module
