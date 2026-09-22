@testset "Mixed backends" begin
    warning = (:warn, r"Mixing Unitful and DynamicQuantities")
    # A fresh logger for each operation also verifies rate-limited warnings.
    @test (@test_logs warning unitless((2UF.u"m",), 6DQ.u"m")) == 3
    @test (@test_logs warning unitless((2DQ.u"m",), 6UF.u"m")) == 3
    basis = (2UF.u"m", 3DQ.u"s")
    @test (@test_logs warning unitless(basis, 4UF.u"m/s")) ≈ 6
    @test (@test_logs warning unitless(basis, 4DQ.u"m/s")) ≈ 6
    @test (@test_logs warning unitless(basis, 5)) == 5
    @test (@test_logs warning unitof(DQ.u"m", (2UF.u"m",))) == 2DQ.u"m"
    @test (@test_logs warning unitless((1UF.u"cm",), DQ.u"m")) ≈ 100
    @test (@test_logs warning unitless((DQ.u"m",), 1UF.u"cm")) ≈ 0.01
    @test (@test_logs warning unitless((2DQ.us"cm",), 1UF.u"m")) ≈ 50
    @test (@test_logs warning unitless((DQ.u"m",), 180UF.u"°")) ≈ π
    @test (@test_logs warning unitof(UF.NoDims, (DQ.u"m",))) == 1
    @test (@test_logs warning LessUnit(basis)(4UF.u"m/s")) ≈ 6
    for target in (UF.u"cm", UF.𝐋, UF.Length, 100UF.u"cm")
        result = @test_logs warning unitof(target, (2DQ.u"m",))
        @test result == 2DQ.u"m"
        @test result isa DQ.Quantity
    end
    target = 1UF.u"m"^(1//11)
    result = @test_logs warning unitof(target, (DQ.u"m",))
    @test DQ.ulength(result) == 1//11
    @test (@test_logs warning unitless((DQ.u"m",), 2UF.u"m")) == 2
    @test_throws ArgumentError unitof(typeof(DQ.u"m"), (1UF.u"m", 1UF.u"s"))

    quantities = Any[2UF.u"m", 6DQ.u"s", 4UF.u"m/s"]
    @test (@test_logs warning unitless.(basis, quantities)) ≈ [1, 2, 6]
    @test (@test_logs warning unitless.(basis, DQ.QuantityArray([2DQ.u"m", 4DQ.u"m"]))) == [1, 2]
    @test_logs warning @test_throws ArgumentError unitless(basis, 1UF.u"kg")
    @test_logs warning @test_throws ArgumentError unitless((1UF.u"m", DQ.u"m"), DQ.u"m")
end

@testset "Generated backend dispatch" begin
    warning = (:warn, r"Mixing Unitful and DynamicQuantities")

    # Neither basis independently selects DQ: the type-valued target must do so.
    targets = (
        typeof(DQ.u"m"),
        DQ.Quantity,
        DQ.Quantity{Float64},
        DQ.RealQuantity,
        DQ.GenericQuantity,
        typeof(DQ.dimension(DQ.u"m")),
        DQ.Dimensions,
        DQ.AbstractDimensions,
        DQ.UnionAbstractQuantity,
    )
    for target in targets, basis in ((), (2UF.u"m",))
        @test_throws ArgumentError unitof(target, basis)
        @test_throws ArgumentError LessUnit(basis)(target)
    end

    # Preserve valid Unitful type targets with either backend.
    for target in (typeof(1UF.u"m"), UF.Length, typeof(UF.𝐋))
        @test unitof(target, (2UF.u"m",)) == 2UF.u"m"
        result = @test_logs warning unitof(target, (2DQ.u"m",))
        @test result == 2DQ.u"m"
    end

    # Target-dimension and basis-quantity conversions have separate warning sites.
    result = @test_logs warning warning unitof(
        UF.u"m", (2UF.u"m", 3DQ.u"s")
    )
    @test result == 2DQ.u"m"

    # With a UF-only basis, backend selection can differ between elements.
    result = @test_logs warning unitless.(
        (2UF.u"m",), Any[2UF.u"m", 4DQ.u"m"]
    )
    @test result == [1, 2]
end
