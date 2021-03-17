module GeometricObjects

using Reexport
@reexport using Tensorial

using Base: @_propagate_inbounds_meta

export
    GeometricObject,
    centroid,
    within,
    distance,
    moment_of_inertia,
# Line
    Line,
# Polygon
    Polygon,
    Rectangle,
# Sphere/Circle
    Sphere,
    Circle,
    radius

include("utils.jl")
include("GeometricObject.jl")
include("Line.jl")
include("Polygon.jl")
include("Sphere.jl")

end # module
