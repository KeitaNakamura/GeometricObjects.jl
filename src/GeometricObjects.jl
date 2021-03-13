module GeometricObjects

using Reexport
@reexport using Tensorial

using Base: @_propagate_inbounds_meta

export
    GeometricObject,
    coordinates,
    center,
    moment_of_inertia,
# Line
    Line,
    perpendicularfoot,
    distance,
    normalunit,
# Polygon
    Polygon,
    Rectangle,
    isinside,
# Sphere/Disk
    Sphere,
    Disk,
    radius

include("utils.jl")
include("GeometricObject.jl")
include("Line.jl")
include("Polygon.jl")
include("Sphere.jl")

end # module
