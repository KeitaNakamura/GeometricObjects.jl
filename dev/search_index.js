var documenterSearchIndex = {"docs":
[{"location":"#GeometricObjects","page":"Home","title":"GeometricObjects","text":"","category":"section"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"pkg> add https://github.com/KeitaNakamura/GeometricObjects.jl.git","category":"page"},{"location":"#Types-and-Functions","page":"Home","title":"Types and Functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Modules = [GeometricObjects]\nOrder   = [:type, :function]","category":"page"},{"location":"#GeometricObjects.Line","page":"Home","title":"GeometricObjects.Line","text":"Line(a::Vec, b::Vec)\nLine(a::Vec => b::Vec)\n\n\n\n\n\n","category":"type"},{"location":"#Base.intersect-Tuple{Line{2}, Line{2}}","page":"Home","title":"Base.intersect","text":"intersect(::Line, ::Line; [extended = (false, false)])\n\nFind intersection point from two lines. Return nothing if not found.\n\n\n\n\n\n","category":"method"},{"location":"#Base.intersect-Union{Tuple{T}, Tuple{dim}, Tuple{Polygon{dim, T}, Line}} where {dim, T}","page":"Home","title":"Base.intersect","text":"intersect(::Polygon, ::Line; [extended = false])\n\nFind the closest intersection point from line to polygon. Return nothing if not found.\n\n\n\n\n\n","category":"method"},{"location":"#Base.intersect-Union{Tuple{T}, Tuple{dim}, Tuple{Polyline{dim, T}, Line}} where {dim, T}","page":"Home","title":"Base.intersect","text":"intersect(::Polyline, ::Line; [extended = false])\n\nFind the closest intersection point from line to polyline. Return nothing if not found.\n\n\n\n\n\n","category":"method"},{"location":"#GeometricObjects.apply_force!-Union{Tuple{dim}, Tuple{GeometricObject{dim}, Vec{dim}, Union{Real, Vec{dim}}, Real}} where dim","page":"Home","title":"GeometricObjects.apply_force!","text":"apply_force!(object::GeometricObject, F, τ, dt)\n\nApply linear force F and torque (moment of force) τ to object with timestep dt.\n\n\n\n\n\n","category":"method"},{"location":"#GeometricObjects.distance-Tuple{Line, Vec}","page":"Home","title":"GeometricObjects.distance","text":"distance(::Line, x::Vec)\ndistance(::Line, x::Vec, threshold::Real)\n\nReturn the distance vector from x to perpendicular foot. When threshold is given, check the contact between line and point x, and return nothing if contact is not detected.\n\njulia> line = Line(@Vec[0.0, 0.0] => @Vec[1.0, 1.0])\nLine{2, Float64}:\n  Coordinates: [[0.0, 0.0], [1.0, 1.0]]\n  Attitude: [1.0, 0.0]\n\njulia> distance(line, @Vec[1.0, 0.0])\n2-element Vec{2, Float64}:\n -0.5\n  0.5\n\njulia> distance(line, @Vec[1.0, 0.0], 1.0)\n2-element Vec{2, Float64}:\n -0.5\n  0.5\n\njulia> distance(line, @Vec[1.0, 0.0], 0.5) === nothing\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#GeometricObjects.distance-Tuple{Polygon{2}, Vec{2}, Real}","page":"Home","title":"GeometricObjects.distance","text":"distance(::Polygon, x::Vec, threshold::Real)\n\nCheck if the polygon is in contact with x with threshold distance, and return the minimum distance from x to the polygon when contact is detected, otherewise return nothing.\n\n\n\n\n\n","category":"method"},{"location":"#GeometricObjects.distance-Union{Tuple{dim}, Tuple{Sphere{dim}, Vec{dim}}} where dim","page":"Home","title":"GeometricObjects.distance","text":"distance(::Sphere, x::Vec)\ndistance(::Sphere, x::Vec, threshold::Real)\n\nReturn the distance vector from x to perpendicular foot on surface of sphere. When threshold is given, check the contact between line and point x, and return nothing if contact is not detected.\n\njulia> sphere = Sphere(Vec(1.0,1.0), 1.0)\nSphere{2, Float64}:\n  Coordinates: [[1.0, 1.0]]\n  Attitude: [1.0, 0.0]\n\njulia> d = distance(sphere, Vec(0.0,0.0,0.0))\n3-element Vec{3, Float64}:\n 0.4226497308103742\n 0.4226497308103742\n 0.4226497308103742\n\njulia> norm(d) ≈ √3 - 1\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#GeometricObjects.isinside-Tuple{Vec{2}, Polygon{2}}","page":"Home","title":"GeometricObjects.isinside","text":"isinside(x::Vec, ::Polygon; include_bounds = true)\n\nCheck if a point isinside a polygon.\n\n\n\n\n\n","category":"method"},{"location":"#GeometricObjects.isinside-Union{Tuple{dim}, Tuple{Vec{dim}, Sphere{dim}}} where dim","page":"Home","title":"GeometricObjects.isinside","text":"isinside(x::Vec, ::Sphere; include_bounds = true)\n\nCheck if a point isinside a sphere.\n\njulia> sphere = Sphere(Vec(1.0,1.0,1.0), 1.0)\nSphere{3, Float64}:\n  Coordinates: [[1.0, 1.0, 1.0]]\n  Attitude: [1.0, 0.0, 0.0]\n\njulia> isinside(Vec(0.5, 0.5, 0.5), sphere)\ntrue\n\njulia> isinside(Vec(0.0, 0.0, 0.0), sphere)\nfalse\n\n\n\n\n\n","category":"method"},{"location":"#GeometricObjects.ison-Tuple{Vec, Line}","page":"Home","title":"GeometricObjects.ison","text":"ison(x::Vec, line::Line)\n\nCheck if a point ison line.\n\nExamples\n\njulia> line = Line(@Vec[0.0, 0.0], @Vec[2.0, 2.0])\nLine{2, Float64}:\n  Coordinates: [[0.0, 0.0], [2.0, 2.0]]\n  Attitude: [1.0, 0.0]\n\njulia> ison(@Vec[1.0, 1.0], line)\ntrue\n\njulia> ison(@Vec[1.0, 0.0], line)\nfalse\n\n\n\n\n\n","category":"method"},{"location":"#GeometricObjects.perpendicularfoot-Tuple{Line, Vec}","page":"Home","title":"GeometricObjects.perpendicularfoot","text":"GeometricObjects.perpendicularfoot(::Line, x::Vec)\n\nReturn the position of perpendicular foot.\n\n\n\n\n\n","category":"method"},{"location":"#GeometricObjects.rotate!-Tuple{GeometricObject, Union{Real, Vec}}","page":"Home","title":"GeometricObjects.rotate!","text":"rotate!(object::GeometricObject, θ::Vec)\n\nRotate object by the angle vector θ. normalize(θ) and norm(θ) should represent the rotation axis and the angle (radian), respectively.\n\n\n\n\n\n","category":"method"},{"location":"#GeometricObjects.translate!-Tuple{GeometricObject, Vec}","page":"Home","title":"GeometricObjects.translate!","text":"translate!(object::GeometricObject, u::Vec)\n\nTranslate object by the displacement u.\n\n\n\n\n\n","category":"method"},{"location":"#GeometricObjects.translate-Tuple{Geometry, Vec}","page":"Home","title":"GeometricObjects.translate","text":"translate(g::Geometry, u::Vec)\n\nTranslate g by the displacement u.\n\n\n\n\n\n","category":"method"},{"location":"#GeometricObjects.update_geometry!-Tuple{GeometricObject, Real}","page":"Home","title":"GeometricObjects.update_geometry!","text":"update_geometry!(::GeometricObject, dt)\n\nUpdate geometry of object by timestep dt. Current linear and angular velocities of object are used in the calculation.\n\n\n\n\n\n","category":"method"},{"location":"#Tensorial.rotate","page":"Home","title":"Tensorial.rotate","text":"rotate(g::Geometry{2}, θ::Real)\nrotate(g::Geometry{3}, θ::Vec)\n\nRotate g by the angle θ. In 3D, normalize(θ) and norm(θ) should represent the rotation axis and the angle (radian), respectively.\n\n\n\n\n\n","category":"function"}]
}
