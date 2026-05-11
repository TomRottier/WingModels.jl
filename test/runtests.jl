using Test
using WingModels

@testset "all" verbose = true begin

    @testset "generic wings" verbose = true begin
        # using Test planform: constant quarter chord at 0.0
        struct TestPlanform <: AbstractPlanform end
        WingModels.quarter_chord(ξ, p::TestPlanform) = 1.0
        WingModels.chord(ξ, p::TestPlanform) = 0.5

        ηs = range(0, 1, length=10)
        ξs = range(0, 1, length=5)
        pl = TestPlanform()
        foreach(ξs) do ξ
            # rectangular aerofoil with height 1
            __aerofoil = [map(η -> [η, ξ, 1.0], ηs); map(η -> [η, ξ, 1.0], ηs)]

            # test scaling
            scaled_aerofoil = WingModels.scale_aerofoil(ξ, __aerofoil, pl)
            @test scaled_aerofoil[1][1] == 0.0
            @test scaled_aerofoil[end][1] == 0.5
            @test scaled_aerofoil[1][3] == 0.5
            @test scaled_aerofoil[end][3] == 0.5

            # test translation (combined with scaling)
            translated_aerofoil = WingModels.translate_aerofoil(ξ, scaled_aerofoil, pl)
            @test translated_aerofoil[1][1] == leading_edge(ξ, pl) == 0.875
            @test translated_aerofoil[end][1] == trailing_edge(ξ, pl) == 1.375
        end


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
        foreach(0:0.1:0.5) do ξ # 0.5x1 rectangle in this region
            @test chord(ξ, generic_planform) == 1.0
            @test quarter_chord(ξ, generic_planform) == 0.0
            @test leading_edge(ξ, generic_planform) == -0.25
            @test trailing_edge(ξ, generic_planform) == 0.75
        end
        foreach(0.5:0.1:1.0) do ξ # bounded by 4x(1-x) parabola in this region
            @test chord(ξ, generic_planform) ≈ 4ξ * (1 - ξ)
            @test quarter_chord(ξ, generic_planform) ≈ 0.0
            @test leading_edge(ξ, generic_planform) ≈ -ξ * (1 - ξ)
            @test trailing_edge(ξ, generic_planform) ≈ 3ξ * (1 - ξ)
        end

        # generic aerofoil: camber line is parabola from 0 to 1 with max camber of 0.25 at 0.5, constant along span
        generic_aerofoil = LiuAerofoil((1.0, 0.0, 0.0), (-1.0, 0.0, 0.0, 0.0), (1.0, 0.0), (1.0, 0.0))
        foreach(0:0.1:1) do ξ
            @test WingModels.LiuWings.max_camber(ξ, generic_aerofoil.zcmax) == 1.0
            @test WingModels.LiuWings.max_thickness(ξ, generic_aerofoil.ztmax) == 1.0
            @test WingModels.LiuWings.camber(0.0, generic_aerofoil.S, WingModels.LiuWings.max_camber(ξ, generic_aerofoil.zcmax)) == 0.0
            @test WingModels.LiuWings.camber(0.5, generic_aerofoil.S, WingModels.LiuWings.max_camber(ξ, generic_aerofoil.zcmax)) == 0.25
            @test WingModels.LiuWings.camber(1.0, generic_aerofoil.S, WingModels.LiuWings.max_camber(ξ, generic_aerofoil.zcmax)) == 0.0
            @test WingModels.LiuWings.thickness(0.0, generic_aerofoil.A, WingModels.LiuWings.max_thickness(ξ, generic_aerofoil.ztmax)) == 0.001
            @test WingModels.LiuWings.thickness(0.5, generic_aerofoil.A, WingModels.LiuWings.max_thickness(ξ, generic_aerofoil.ztmax)) == -(0.5^2 - sqrt(0.5))
            @test WingModels.LiuWings.thickness(1.0, generic_aerofoil.A, WingModels.LiuWings.max_thickness(ξ, generic_aerofoil.ztmax)) == 0.001
            @test aerofoil_pt(0.0, ξ, generic_aerofoil; upper=false) == [0.0, ξ, -0.001]
            @test aerofoil_pt(0.5, ξ, generic_aerofoil; upper=false) == [0.5, ξ, 0.25 - -(0.5^2 - sqrt(0.5))]
            @test aerofoil_pt(0.5, ξ, generic_aerofoil; upper=true) == [0.5, ξ, 0.25 + -(0.5^2 - sqrt(0.5))]
            @test aerofoil_pt(1.0, ξ, generic_aerofoil; upper=false) == [1.0, ξ, -0.001]
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

        # test values taken from airfoiltools.com
        n = 100
        pts = range(1, 0, step=-1 / (n ÷ 2))
        idxs = range(8, 50, step=13)
        ηs = pts[idxs]
        @test aerofoil_pt(0.0, 0.0, af; upper=false) == [0.0, 0.0]
        @test aerofoil_pt(ηs[1], 0.0, af; upper=true) ≈ [0.860952, 0.026875] atol = 1e-6
        @test aerofoil_pt(ηs[2], 0.0, af; upper=true) ≈ [0.601010, 0.063237] atol = 1e-6
        @test aerofoil_pt(ηs[3], 0.0, af; upper=true) ≈ [0.339105, 0.079199] atol = 1e-6
        @test aerofoil_pt(ηs[4], 0.0, af; upper=true) ≈ [0.076565, 0.050135] atol = 1e-6
        pts = range(0, 1, step=1 / (n ÷ 2))
        idxs = range(8, 50, step=13)
        ηs = pts[idxs]
        @test aerofoil_pt(ηs[1], 0.0, af; upper=false) ≈ [0.143397, -0.040719] atol = 1e-6
        @test aerofoil_pt(ηs[2], 0.0, af; upper=false) ≈ [0.400000, -0.037998] atol = 1e-6
        @test aerofoil_pt(ηs[3], 0.0, af; upper=false) ≈ [0.658840, -0.023917] atol = 1e-6
        @test aerofoil_pt(ηs[4], 0.0, af; upper=false) ≈ [0.919362, -0.006059] atol = 1e-6

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
        @test aerofoil_pt(0.0, rand(), af; upper=true) == 0.1
        @test aerofoil_pt(0.0, rand(), af; upper=false) == -0.1
        @test aerofoil_pt(1.0, rand(), af; upper=true) == 0.1
        @test aerofoil_pt(1.0, rand(), af; upper=false) == -0.1
        @test aerofoil_pt(0.5, rand(), af; upper=true) == 0.1
        @test aerofoil_pt(0.5, rand(), af; upper=false) == -0.1
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