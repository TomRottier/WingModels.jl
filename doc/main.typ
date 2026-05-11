= WingModels.jl

WingModels.jl is a Julia package for generating and analysing wing geometries. Wing models can be specified through combinations of different aerofoils and planforms and the resulting model can be exported as a 3D mesh object. Geometric and aerodynamic properties of the wing models can also be calculated.


== Implementation
A wing can be defined with two independent parts:

- Aerofoil distribution: specifies the cross-section of the wing along a spanwise axis.
- Planform: specifies the transformation (translation + scale) of the aerofoil along the spanwise axis.

A coordinate system ${x,y,z}$ with origin at the root quarter chord defines the wing geometry, $x$ points in the chordwise direction with positive towards the leading edge, $y$ points in the spanwise direction with positive towards the wingtip, and $z$ perpendicular to $x \, y$ following a right-handed system.

An aerofoil distribution thus defines a function $z = f(x)$ for the continuous curve describing the wing cross-section at this spanwise position $y$. The planform similarly defines a function $x = g(y)$ for the continuous curve describing the wing planform.


// z = f(x) and x = g(y) plots for aerofoil and planform


The aerofoil distribution is given relative to the local chord length (leading edge: x = 0 trailing edge x = 1). The planform distribution is given relative to the wingspan (wing root y = 0 wing tip y = 1). The planform scales and translates the aerofoil along the span to produce a 3D wing scaled by the wingspan
