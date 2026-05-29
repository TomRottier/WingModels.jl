# WingModels.jl

## Installation
To install this package run:
```julia
    import Pkg; Pkg.add("https://github.com/TomRottier/WingModels.jl")
```

## Background
A wing geometry can be defined in two parts: 
- an aerofoil distribution: defines the aerofoil shape at each spanwise location. 
- a planform: defines a spanwise chord distribution and a quarter chord line.
    
The whole wing is then given by scaling the aerofoil at each spanwise location by the chord length and translating it onto the quarter chord line. By conforming to this specification arbitrary wing geometries can be created by specifying the aerofoil distribution and the planform.

All aerofoil points are normalised by the chord length and spanwise points are normalised by the span length. This means that chordwise aerofoil points return 0.0 for the leading edge and 1.0 for the trailing edge. Similarly in the spanwise direction, 0.0 is the wing root and 1.0 is the wing tip. To produce a wing with correct dimensions these points must be scaled correctly before outputting to a file.

## Usage
This package provides the tools to generate a wing geometry and output into various formats. Some common aerofoil and planform parameterisations are included but you can define your own as well, see below.

First import the package:

```julia
    using WingModels
```

### Create simple rectangular wing with NACA0012 aerofoil
First, create the aerofoil distribution, for this example we are using the NACA0012 aerofoil which remains the same along the span. This is included with this package as the `NACA00XX` type: 

```julia
    af = NACA00XX(12) # a constant aerofoil distribution using a NACA0012 aerofoil
```
the points defining this aerofoil at a specific spanwise location can be returned with the following function:
```julia
    aerofoil(0.0, af; n=50) # returns 100 points describing this aerofoil at the wing root (0.0)
    aerofoil(1.0, af; n=50) # returns 100 points describing this aerofoil at the wing tip (1.0)
    aerofoil(0.5, af; n=50) # returns 100 points describing this aerofoil at mid-span (0.5)
```
as this aerofoil distribution is constant along the span, these three function calls return the same points. To return a specific chordwise coordinate use the following function instead:
```julia
    aerofoil_height(0.0, 0.0, af; upper=true) # returns the height at the leading edge on the upper surface
    aerofoil_height(0.5, 0.0, af; upper=false) # returns the height at mid-chord on the lower surface
```

Next we create the wing planform, which is just a simple, rectangular wing with aspect ratio of 4:
```julia
    pl = RectangularPlanform(0.25) # the argument is the (constant) chord length normalised by the span length, which for a constant chord planform is just the inverse of the aspect ratio
```
the leading and trailing edges at a spanwise location can be queried with:
```julia
    leading_edge(0.0, pl) # returns chordwise coordinate of leading edge at wing root (0.0)
    trailing_edge(1.0, pl) # returns chordwise coordinate of trailing edge at wing tip (1.0)
```
and the chord length at spanwise location can be queried with:
```julia
    chord(0.5, pl) # returns chord length at midspan
    chord(0.125, pl) # returns the chord length 12.5% along the the span
```
which is of course the same anywhere along the span for a rectangular planform.

The whole wing is can then be created:
```julia
    w = Wing(pl, af)
```
to generate a vector of points on the surface of this wing:
```julia
    wing_pts = wing(w, 0.0, 1.0; nchord=100, nspan=50) # create wing from root (0.0) to tip (1.0)
```
the keyword arguments `nchord` and `nspan` specify the number of points per aerofoil and the number of aerofoils along the span, respectively. 

The wing, aerofoil, and planform can be easily visualised. For example using the `Makie` plotting library:
```julia
    using GLMakie

    # root aerofoil
    pts = Point3.(aerofoil(0.0, af; nchord=100))
    lines(pts)

    # wing planform
    pts = Point3.(planform(pl, n=100))
    lines(pts)

    # full wing
    pts = Point3.(wing(w, 0.0, 1.0; nchord=100, nspan=50))
    lines(pts)

```

The wing can be exported as either a list of points in a text file or as an STL:
```julia
    write_pts("wing.txt", w; nchord=100, nspan=50, sf=100, delim=" ") # export to txt using space as delimiter
    write_stl("wing.stl", w; nchord=100, nspan=50, sf=100) # export stl 
```
where the keyword argument `sf` defines the scale factor. To produce a wing with a 100mm span, set this argument to 100, this will scale each point by a factor of 100 (so e.g. now the wing tip is located at 100 units in the spanwise direction).

Other functions that may be of use: `area`, `mean_chord`, `aspect_ratio`. These do as their name suggests albeit in scaled, non-dimensional units. Check their docstrings for conversions.

### Predefined aerofoil distributions and planforms
- Aerofoil distributions: NACA 00XX series, NACA 4 digit series
- Planform distributions: rectangular, elliptical, trapezoidal, empirical

An empirical planform is just defined from a vector of points describing its leading and trailing edges and can be used for example if digitising a planform from an image. See its documentation for specifics.

### Avian wings from Liu et al., 2006
Also included are the aerofoil distribution and planform `LiuAerofoil` and `LiuPlanform` which are based on the paper *Liu et al. 2006 Avian wing geometry and kinematics. AIAA Journal*. These are parameterisations of avian wings taken from scans of a seagull, merganser, teal, and owl wing. The data provided in the paper to replicate the aerofoil distributions is provided here for convenience (except for the owl). See the paper for details on the wing parameterisation.

### Using your own wing parameterisations
To create your own planform type you must create a new subtype of `AbstractPlanform` and then define the following methods:

```julia
    # define your own planform subtype
    struct YourPlanform <: AbstractPlanform
        data # data needed for your planform
    end

    # add methods for your planform, see docstrings for requirements
    Wings.quarter_chord(y, pl::YourPlanform) = ...
    Wings.chord(y, pl::YourPlanform) = ...
```

Create your own aerofoil distribution in a similar way
```julia
    # define your own aerofoil type
    struct YourAerofoil <:AbstractAerofoil
        data # data needed for your aerofoils
    end

    # add methods for you aerofoils, see docstrings for requirements
    Wings.aerofoil_height(x, y, af::YourAerofoil; upper) = ...
```

The full wing can then be used just as before:
```julia
    af = YourAerofoil(aerofoil_data...)
    pl = YourPlanform(planform_data...)
    w = Wing(pl, af)

    wing_pts = wing(w, 0.0, 1.0; nchord=100, nspan=50)
    write_stl("custom_wing.stl", w; nchord=50, nspan=20, sf=200)
```