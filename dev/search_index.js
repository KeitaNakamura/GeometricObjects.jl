var documenterSearchIndex = {"docs":
[{"location":"#GeometricObjects","page":"Home","title":"GeometricObjects","text":"","category":"section"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"pkg> add https://github.com/KeitaNakamura/GeometricObjects.jl.git","category":"page"},{"location":"#Types-and-Functions","page":"Home","title":"Types and Functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Modules = [GeometricObjects]\nOrder   = [:type, :function]","category":"page"},{"location":"#GeometricObjects.Line","page":"Home","title":"GeometricObjects.Line","text":"Line(a::Vec, b::Vec)\nLine(a::Vec => b::Vec)\n\n\n\n\n\n","category":"type"},{"location":"#Base.in-Tuple{Vec{2, T} where T, Polygon{2, T, V} where {T, V<:AbstractArray{Vec{2, T}, 1}}}","page":"Home","title":"Base.in","text":"in(x::Vec, ::Polygon; include_bounds = true)\n\nCheck if x is in a polygon.\n\n\n\n\n\n","category":"method"},{"location":"#Base.in-Tuple{Vec{dim, T} where {dim, T}, GeometricObjects.AbstractLine}","page":"Home","title":"Base.in","text":"in(x::Vec, line::AbstractLine)\n\nCheck if x is in line.\n\nExamples\n\njulia> line = Line(@Vec[0.0, 0.0], @Vec[2.0, 2.0])\n2-element Line{2,Float64}:\n [0.0, 0.0]\n [2.0, 2.0]\n\njulia> @Vec[1.0, 1.0] in line\ntrue\n\njulia> @Vec[1.0, 0.0] in line\nfalse\n\n\n\n\n\n","category":"method"},{"location":"#Base.in-Union{Tuple{dim}, Tuple{Vec{dim, T} where T, Sphere{dim, T} where T}} where dim","page":"Home","title":"Base.in","text":"in(x::Vec, ::Sphere; include_bounds = true)\n\nCheck if x is in a sphere.\n\njulia> sphere = Sphere(Vec(1.0,1.0), 1.0)\n1-element Sphere{2,Float64}:\n [1.0, 1.0]\n\njulia> Vec(0.5, 0.5) in sphere\ntrue\n\njulia> Vec(0.0, 0.0) in sphere\nfalse\n\n\n\n\n\n","category":"method"},{"location":"#GeometricObjects.distance-Tuple{GeometricObjects.AbstractLine, Vec{dim, T} where {dim, T}}","page":"Home","title":"GeometricObjects.distance","text":"distance(::AbstractLine, x::Vec)\ndistance(::AbstractLine, x::Vec, threshold::Real)\n\nCompute the distance vector from x to perpendicular foot. When threshold is given, check the contact between line and point x, and return nothing if contact is not detected. Note that if the perpendicular foot does not lie on the line, contact detection is performed using distance between x and vertices of line.\n\njulia> line = Line(@Vec[0.0, 0.0] => @Vec[1.0, 1.0])\n2-element Line{2,Float64}:\n [0.0, 0.0]\n [1.0, 1.0]\n\njulia> distance(line, @Vec[1.0, 0.0])\n2-element Tensor{Tuple{2},Float64,1,2}:\n -0.5\n  0.5\n\njulia> distance(line, @Vec[1.0, 0.0], 1.0)\n2-element Tensor{Tuple{2},Float64,1,2}:\n -0.5\n  0.5\n\njulia> distance(line, @Vec[1.0, 0.0], 0.5) === nothing\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#GeometricObjects.distance-Union{Tuple{dim}, Tuple{Sphere{dim, T} where T, Vec{dim, T} where T}} where dim","page":"Home","title":"GeometricObjects.distance","text":"distance(::Sphere, x::Vec)\ndistance(::Sphere, x::Vec, threshold::Real)\n\nCompute the distance vector from x to perpendicular foot on surface of sphere. When threshold is given, check the contact between line and point x, and return nothing if contact is not detected.\n\njulia> sphere = Sphere(Vec(1.0,1.0), 1.0)\n1-element Sphere{2,Float64}:\n [1.0, 1.0]\n\njulia> d = distance(sphere, Vec(0.0,0.0))\n2-element Tensor{Tuple{2},Float64,1,2}:\n 0.29289321881345254\n 0.29289321881345254\n\njulia> norm(d) ≈ √2 - 1\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#GeometricObjects.perpendicularfoot-Tuple{GeometricObjects.AbstractLine, Vec{dim, T} where {dim, T}}","page":"Home","title":"GeometricObjects.perpendicularfoot","text":"GeometricObjects.perpendicularfoot(::AbstractLine, x::Vec)\n\nCompute the position of perpendicular foot.\n\n\n\n\n\n","category":"method"}]
}
