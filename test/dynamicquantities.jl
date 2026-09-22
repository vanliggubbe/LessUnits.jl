struct CookiesAndMilk{R} <: DQ.AbstractDimensions{R}
    cookies :: R
    milk :: R
end

@testset "DynamicQuantities" begin
    basis = (2DQ.u"m", 3DQ.u"s")
    @test_logs unitless(basis, 6DQ.u"m")
    @test unitless(basis, 6DQ.u"m") == 3
    @test unitless(basis, 4DQ.u"m/s") ≈ 6
    @test unitof(DQ.u"m/s", basis) ≈ (2/3)DQ.u"m/s"
    @test unitof(100DQ.u"m/s", basis) == unitof(DQ.u"m/s", basis)
    @test unitof(DQ.dimension(DQ.u"m"), basis) == 2DQ.u"m"
    @test unitof(DQ.u"m", basis) isa DQ.Quantity
    @test unitless(basis, 5) == 5
    @test unitless(basis, DQ.Quantity(5)) == 5
    @test unitof(DQ.Quantity(5), basis) == 1
    @test unitof(DQ.Quantity(5), basis) isa Number
    @test !(unitof(DQ.Quantity(5), basis) isa DQ.UnionAbstractQuantity)
    @test unitless((DQ.u"m", 2DQ.u"m"), DQ.Quantity(5)) == 5
    @test unitless((), DQ.Quantity(3)) == 3
    @test_throws ArgumentError unitless((), DQ.u"m")

    @testset "physical reference basis" begin
        reference = (2DQ.Constants.e, DQ.Constants.hbar, 2π * DQ.u"GHz", DQ.Constants.k_B)
        flux = DQ.Constants.h / (2DQ.Constants.e)
        @test unitof(DQ.u"C", reference) ≈ 2DQ.Constants.e
        @test unitless(reference, flux) ≈ 2π
        @test unitless((DQ.Constants.e, reference[2:end]...), flux) ≈ π
        @test unitof(DQ.u"K", reference) ≈ DQ.Constants.h * DQ.u"GHz" / DQ.Constants.k_B
        @test LessUnit(reference)(DQ.u"K") > 0
    end

    @testset "dimensions and scales" begin
        @test unitless((DQ.u"cm",), DQ.u"m") ≈ 100
        @test unitof(DQ.u"m", (DQ.Quantity(4, length=2),)) == 2DQ.u"m"
        @test unitof(DQ.u"m", (DQ.Quantity(2048, length=11),)) ≈ 2DQ.u"m"
        @test DQ.dimension(unitof(DQ.u"m", (DQ.Quantity(2048, length=11),))) == DQ.dimension(DQ.u"m")
        rational_target = DQ.Quantity(1, DQ.Dimensions{Rational{Int}}(length=1//11))
        @test DQ.dimension(unitof(rational_target, (DQ.u"m",))) == DQ.dimension(rational_target)
        @test unitless((DQ.Quantity(big"2.0", length=1),), DQ.Quantity(big"6.0", length=1)) isa BigFloat
        @test unitless((DQ.RealQuantity(2.0, length=1),), DQ.RealQuantity(6.0, length=1)) == 3
        @test unitless((DQ.GenericQuantity(2.0, length=1),), DQ.GenericQuantity(6.0, length=1)) == 3
        @test unitless((2DQ.u"m",), (2+4im)DQ.u"m") == 1+2im
        @test_throws ArgumentError unitof(DQ.u"kg", basis)
        @test_throws ArgumentError unitof(DQ.u"s", (DQ.u"Hz", DQ.u"s"))
        @test_throws ArgumentError unitof(DQ.u"m", (DQ.u"m/s",))
        @test_throws ArgumentError unitof(DQ.u"m", (DQ.Quantity(1), DQ.u"m"))
        @test_throws ArgumentError unitof(typeof(DQ.u"m"), basis)
        @test_throws ArgumentError unitof(DQ.Dimensions, basis)
    end

    @testset "symbolic quantities" begin
        symbolic_basis = (2DQ.us"cm", 3DQ.us"s")
        @test unitless(symbolic_basis, DQ.us"m") ≈ 50
        @test unitless(basis, 400DQ.us"cm") ≈ 2
        @test unitof(DQ.us"cm", basis) == 2DQ.u"m"
        @test unitof(DQ.dimension(DQ.us"cm"), basis) == 2DQ.u"m"
        @test unitof(DQ.u"m", symbolic_basis) ≈ 0.02DQ.u"m"
        @test_throws ArgumentError unitof(DQ.u"m", (DQ.us"cm", DQ.us"m"))
    end

    @testset "custom dimension interface" begin
        reference = (DQ.Quantity(2.0, CookiesAndMilk(cookies=1)), DQ.Quantity(3.0, CookiesAndMilk(milk=1)))
        target = DQ.Quantity(4.0, CookiesAndMilk(cookies=1, milk=-1))
        @test unitless(reference, target) ≈ 6
        @test unitof(DQ.dimension(target), reference) ≈ DQ.Quantity(2/3, DQ.dimension(target))
    end

    @testset "LessUnit and broadcasting" begin
        lu = LessUnit(basis)
        @test lu(6DQ.u"m") == 3
        @test LessUnit(basis...)(6DQ.u"m") == 3
        @test lu(DQ.dimension(DQ.u"m")) == 2DQ.u"m"
        @test_throws ArgumentError lu(typeof(DQ.u"m"))
        vector = [2DQ.u"m", 4DQ.u"m", 6DQ.u"m"]
        matrix = [2DQ.u"m" 4DQ.u"m"; 6DQ.u"m" 8DQ.u"m"]
        tuple = (2DQ.u"m", 6DQ.u"s", 4DQ.u"m/s")
        for quantities in (vector, matrix, tuple, typeof(DQ.u"m")[], DQ.QuantityArray(vector))
            @test unitless.(basis, quantities) == map(Base.Fix1(unitless, basis), quantities)
        end
        @test unitless.(Ref(basis), vector) == [1, 2, 3]
        @test 2 .* unitless.(basis, vector .+ 2DQ.u"m") .+ 1 == [5, 7, 9]
        array = DQ.QuantityArray(vector)
        @test 2 .* unitless.(basis, array .+ 2DQ.u"m") .+ 1 == [5, 7, 9]
        destination = zeros(3)
        destination .= unitless.(basis, array)
        @test destination == [1, 2, 3]
        @test unitless.(basis, DQ.QuantityArray(Float64[], DQ.u"m")) == Float64[]
    end
end
