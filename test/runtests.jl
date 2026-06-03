using Test
using WingModels

@testset "all" verbose = true begin

    @testset "interface" verbose = true begin
        for af in [
            NACA4(0.1, 0.2, 0.3),
            NACA00(1),
            RectangularAerofoil(0.1),
            seagull.aerofoil,
        ]
            @test length(aerofoil_height(0.0, 0.0, af; upper=true)) == 1
            @test length(aerofoil(0.0, af)[1]) == 2
            @test length(aerofoil(0.0, af; n=24)) == 48
        end

        for pl in [
            TrapezoidalPlanform(0.5, 0.4, 0.3, 0.5),
            EllipticalPlanform(0.3, 0.5),
            RectangularPlanform(0.6),
            TriangularPlanform(0.5, 0.4),
            seagull.planform
        ]
            @test length(chord(rand(), pl)) == 1
            @test length(quarter_chord(rand(), pl)) == 1
            @test length(leading_edge(rand(), pl)) == 1
            @test length(trailing_edge(rand(), pl)) == 1
            @test length(planform(pl)[1]) == 2
            @test length(planform(pl; n=47)) == 94

        end
    end

    @testset "wing transformation" verbose = true begin
        af = RectangularAerofoil(0.2)
        pl = RectangularPlanform(0.5)

        y = rand()
        af_pts = aerofoil(y, af; n=49)
        af_pts_sc = WingModels.scale_aerofoil(y, af_pts, pl)
        @test af_pts_sc[1][1] == 0.0 # wing leading edge has x = 0
        @test af_pts_sc[49][1] == 0.5 # wing trailing edge has 1 * chord length
        @test af_pts_sc[10][2] == af_pts[10][2] * 0.5
        @test af_pts_sc[70][2] == af_pts[70][2] * 0.5

        af_pts_tr = WingModels.translate_aerofoil(y, af_pts_sc, pl)
        @test af_pts_tr[1][1] == -0.125 # leading edge at -0.25c
        @test af_pts_tr[13][1] == 0.0 # quarter chord at x = 0
        @test af_pts_tr[49][1] == 0.375 # trailing edge at 0.75d
        @test af_pts_tr[10][2] == y
        @test af_pts_tr[15][3] == af_pts_sc[15][2] # translation shouldnt change height

    end


    @testset "wing properties" verbose = true begin
        # area of rectangular wing with = c₀
        for c in 0.1:0.1:1.0
            pl = RectangularPlanform(c)
            @test area(pl) == c
            @test mean_chord(pl) == c
            @test aspect_ratio(pl) == 1 / c
            @test second_moment_of_area(pl) ≈ c / 3
        end

        pl = RectangularPlanform(1.0)
        @test mean_aerodynamic_chord(pl) == 1.0
        @test second_moment_of_area(pl) ≈ 1 / 3
        pl = TriangularPlanform(1.0, π / 4)
        @test mean_chord(pl) ≈ 0.5
        @test mean_aerodynamic_chord(pl) ≈ 2 / 3
        @test second_moment_of_area(pl) ≈ 1 / 12
    end


    @testset "Liu" verbose = true begin
        # generic planform
        generic_planform = LiuPlanform(0.5, 0.0, (0.0, 0.0, 0.0, 0.0, 0.0), 1.0)
        foreach(0:0.1:0.5) do y # 0.5x1 rectangle in this region
            @test chord(y, generic_planform) == 1.0
            @test quarter_chord(y, generic_planform) == 0.0
            @test leading_edge(y, generic_planform) == -0.25
            @test trailing_edge(y, generic_planform) == 0.75
        end
        foreach(0.5:0.1:1.0) do y # bounded by 4x(1-x) parabola in this region
            @test chord(y, generic_planform) ≈ 4y * (1 - y)
            @test quarter_chord(y, generic_planform) ≈ 0.0
            @test leading_edge(y, generic_planform) ≈ -y * (1 - y)
            @test trailing_edge(y, generic_planform) ≈ 3y * (1 - y)
        end

        # generic aerofoil: camber line is parabola from 0 to 1 with max camber of 0.25 at 0.5, constant along span
        generic_aerofoil = LiuAerofoil((1.0, 0.0, 0.0), (-1.0, 0.0, 0.0, 0.0), (1.0, 0.0), (1.0, 0.0))
        foreach(0:0.1:1) do y
            @test WingModels.LiuWings.max_camber(y, generic_aerofoil.zcmax) == 1.0
            @test WingModels.LiuWings.max_thickness(y, generic_aerofoil.ztmax) == 1.0
            @test WingModels.LiuWings.camber(0.0, generic_aerofoil.S, WingModels.LiuWings.max_camber(y, generic_aerofoil.zcmax)) == 0.0
            @test WingModels.LiuWings.camber(0.5, generic_aerofoil.S, WingModels.LiuWings.max_camber(y, generic_aerofoil.zcmax)) == 0.25
            @test WingModels.LiuWings.camber(1.0, generic_aerofoil.S, WingModels.LiuWings.max_camber(y, generic_aerofoil.zcmax)) == 0.0
            @test WingModels.LiuWings.thickness(0.0, generic_aerofoil.A, WingModels.LiuWings.max_thickness(y, generic_aerofoil.ztmax)) == 0.001
            @test WingModels.LiuWings.thickness(0.5, generic_aerofoil.A, WingModels.LiuWings.max_thickness(y, generic_aerofoil.ztmax)) == -(0.5^2 - sqrt(0.5))
            @test WingModels.LiuWings.thickness(1.0, generic_aerofoil.A, WingModels.LiuWings.max_thickness(y, generic_aerofoil.ztmax)) == 0.001
            @test aerofoil_height(0.0, y, generic_aerofoil; upper=false) == -0.001
            @test aerofoil_height(0.5, y, generic_aerofoil; upper=false) == 0.25 - -(0.5^2 - sqrt(0.5))
            @test aerofoil_height(0.5, y, generic_aerofoil; upper=true) == 0.25 + -(0.5^2 - sqrt(0.5))
            @test aerofoil_height(1.0, y, generic_aerofoil; upper=false) == -0.001
        end
    end

    @testset "naca" verbose = true begin
        af = NACA4(2, 4, 12) # maximum camber of 0.02 at 0.04, thickness of 0.12
        # test camber at max camber location
        @test WingModels.naca4_aerofoil_camber_front(af.p / 100, af.m / 100, af.p / 100) == af.m / 100
        @test WingModels.naca4_aerofoil_camber_back(af.p / 100, af.m / 100, af.p / 100) == af.m / 100

        # test is maximum at max camber location
        @test WingModels.naca4_aerofoil_camber_gradient_front(af.p / 100, af.m / 100, af.p / 100) == 0.0
        @test WingModels.naca4_aerofoil_camber_gradient_back(af.p / 100, af.m / 100, af.p / 100) == 0.0

        # test values taken from airfoiltools.com: NACA2412, n points 20, cosine spacing
        pts = aerofoil(0.0, af; n=11)
        @test pts[1] == [0.0, 0.0]
        @test pts[2] ≈ [0.022051, 0.028152] atol = 1e-6
        @test pts[5] ≈ [0.344680, 0.079180] atol = 1e-6
        @test pts[8] ≈ [0.795047, 0.037760] atol = 1e-6
        @test pts[11] ≈ [1.0, 0.0] atol = 1e-6
        @test pts[15] ≈ [0.792738, -0.014999] atol = 1e-6
        @test pts[18] ≈ [0.346303, -0.039923] atol = 1e-6
        @test pts[21] ≈ [0.026892, -0.023408] atol = 1e-6
    end

    @testset "geometric wings" verbose = true begin
        pl = TrapezoidalPlanform(0.5, 0.0, 1.0, 1.0) # 1.0x1.0 recatngular planform
        @test chord(0.0, pl) == 1.0
        @test chord(0.5, pl) == 1.0
        @test chord(1.0, pl) == 1.0

        pl = TrapezoidalPlanform(0.5, 0.0, 1.0, 0.0) # triangular planform
        @test chord(0.0, pl) == 1.0
        @test chord(0.5, pl) == 0.5
        @test chord(1.0, pl) == 0.0

        pl = TrapezoidalPlanform(0.5, 0.0, 1.0, 0.5) # trapezoidal planform
        @test chord(0.0, pl) == 1.0
        @test chord(0.3, pl) == 1.0
        @test chord(0.75, pl) == 0.5

        pl = EllipticalPlanform(0.5, 0.5) # symmetrical planform
        @test chord(0.0, pl) ≈ 1.0
        @test chord(0.6, pl) ≈ 0.8
        @test chord(0.8, pl) ≈ 0.6
        @test chord(1.0, pl) ≈ 0.0

        pl = EllipticalPlanform(0.3, 0.1) # non symmetrical planform
        @test chord(0.0, pl) ≈ 0.4
        @test chord(0.6, pl) ≈ 0.32
        @test chord(0.8, pl) ≈ 0.24
        @test chord(1.0, pl) ≈ 0.0

        pl = RectangularPlanform(0.3) # rectangular planform
        @test chord(0.0, pl) == 0.3
        @test chord(0.5, pl) == 0.3
        @test chord(1.0, pl) == 0.3
        @test quarter_chord(0.0, pl) == 0.0
        @test quarter_chord(0.5, pl) == 0.0
        @test quarter_chord(1.0, pl) == 0.0

        pl = TriangularPlanform(1.0, π / 4) # triangular planform
        @test chord(0.0, pl) == 1.0
        @test chord(0.5, pl) == 0.5
        @test chord(1.0, pl) == 0.0
        @test quarter_chord(0.0, pl) == 0.0
        @test quarter_chord(1.0, pl) ≈ -0.75

        af = RectangularAerofoil(0.1) # rectangular aerofoil
        @test aerofoil_height(0.0, rand(), af; upper=true) == 0.1
        @test aerofoil_height(0.0, rand(), af; upper=false) == -0.1
        @test aerofoil_height(1.0, rand(), af; upper=true) == 0.1
        @test aerofoil_height(1.0, rand(), af; upper=false) == -0.1
        @test aerofoil_height(0.5, rand(), af; upper=true) == 0.1
        @test aerofoil_height(0.5, rand(), af; upper=false) == -0.1
    end

    @testset "empirical" verbose = true begin
        # dummy aerofoil data
        _af = NACA4(5, 2, 11)
        xs = range(0, 1; step=0.1)
        _af_pts = aerofoil(0.0, _af; n=length(xs))
        _af_pts_upper = [aerofoil_height(x, 0.0, _af; upper=true) for x in xs]
        _af_pts_lower = [aerofoil_height(x, 0.0, _af; upper=false) for x in xs]
        data = hcat(xs, _af_pts_upper, _af_pts_lower)

        data2 = data * 5rand()
        data2[:, 1] .+ 2rand()
        data2_sc = WingModels.scale_data(data2)
        af = EmpiricalAerofoil(data2, 1.0) # no smoothing
        @test WingModels.smoothing_number(1.0, data2, 3) == size(data2, 1)
        @test WingModels.smoothing_number(0.0, data2, 3) == 4
        @test data2_sc ≈ data

        for x in xs
            @test aerofoil_height(x, 0.0, af; upper=true) ≈ aerofoil_height(x, 0.0, _af; upper=true)
            @test aerofoil_height(x, 0.0, af; upper=false) ≈ aerofoil_height(x, 0.0, _af; upper=false)
        end

        # dummy planform data
        _pl = EllipticalPlanform(0.1, 0.3)
        ys = range(0, 1; step=0.1)
        pts = hcat(ys, [leading_edge(y, _pl) for y in ys], [trailing_edge(y, _pl) for y in ys]) # dummy data
        pl = EmpiricalPlanform(pts, 1.0)
        for f in (leading_edge, trailing_edge, chord, quarter_chord)
            @test f(0.0, pl) ≈ f(0.0, _pl)
            @test f(0.5, pl) ≈ f(0.5, _pl)
            @test f(1.0, pl) ≈ f(1.0, _pl)
        end
    end

    @testset "export" verbose = true begin
        sg = WingModels.seagull
        sga = sg.aerofoil
        sgp = sg.planform
        af = LiuAerofoil(sga.S, sga.A, sga.zcmax, sga.ztmax)
        pf = LiuPlanform(0.5, 0.6, sgp.E, sgp.c₀)
        nchord, nspan = 10, 5
        w = wing(pf, af; nchord, nspan)
        conns = WingModels.get_conns(nchord, nspan)
        @test conns[1] == (1, 11, 12)
        @test conns[2] == (12, 2, 1)
        @test conns[19] == (10, 20, 11)
        @test conns[20] == (11, 1, 10)
    end
end