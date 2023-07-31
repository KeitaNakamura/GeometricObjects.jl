module GeometricObjects

using Reexport
@reexport using Tensorial
@reexport using WriteVTK
using StaticArrays
using RecipesBase
using DelimitedFiles

using Base: @_propagate_inbounds_meta

import Tensorial: quaternion

import LinearAlgebra: norm
export norm

export
# Geometry
    Geometry,
    geometry,
    centroid,
    area,
    volume,
    distance,
    moment_of_inertia,
    translate!,
    rotate!,
    coordinates,
    num_coordinates,
# Line
    Line,
    ison,
# Polygon
    Polygon,
    Rectangle,
    isinside,
# Polyline
    Polyline,
# Sphere/Circle
    Sphere,
    Circle,
    radius,
# FinitePoints
    FinitePoints,
# GeometricObject,
    GeometricObject,
    update_geometry!,
    apply_force!,
    velocityat,
# VTK
    vtk_grid

include("Geometry.jl")
include("Line.jl")
include("Polygon.jl")
include("Polyline.jl")
include("Sphere.jl")
include("FinitePoints.jl")
include("GeometricObject.jl")
include("plots.jl")
include("vtk.jl")

end # module
