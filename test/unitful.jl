@testset "Unitful" begin
    @test_logs unitless((2UF.u"m",), 6UF.u"m")
    @test (@inferred unitless((2UF.u"m",), 6UF.u"m")) == 3
    @test (@inferred unitof(UF.Length, (2UF.u"m",))) == 2UF.u"m"
    @test unitof(UF.Charge, (2UF.u"q", 1UF.u"ħ", 2π * 1UF.u"GHz")) ≈ 2UF.u"q"
    @test unitless((2UF.u"q", 1UF.u"ħ", 2π * 1UF.u"GHz"), 1UF.u"Φ0") ≈ 2π
    @test unitless((1UF.u"q", 1UF.u"ħ", 2π * 1UF.u"GHz"), 1UF.u"Φ0") ≈ π
    @test unitless((2UF.u"q", 1UF.u"cm"), 5) ≈ 5
    @test unitof(UF.Temperature, (2UF.u"q", 1UF.u"ħ", 2π * 1UF.u"GHz", 1UF.u"k")) ≈ 1UF.u"h*GHz/k"
    @test LessUnit(2π * 1UF.u"GHz", 2UF.u"q", 1UF.u"ħ", 1UF.u"k")(1UF.u"K") > 0.0
    @test_throws ArgumentError unitof(UF.𝐋, (2UF.u"q", 1UF.u"ns"))
    @test_throws ArgumentError unitof(UF.𝐓, (1UF.u"GHz", 1UF.u"ns"))
    @test_throws ArgumentError unitof(UF.Charge, (2UF.u"q", 1UF.u"fF", 1UF.u"nH", 1UF.u"MHz"))
    @test unitof(UF.u"s", (2π * 1UF.u"GHz", 2UF.u"q")) ≈ inv(2π) * UF.u"ns"
    @test UF.unit(unitof(UF.u"s", (2π * 1UF.u"GHz", 2UF.u"q"))) == UF.u"s"

    @testset "unitless broadcasting" begin
        basis = (2UF.u"m", 3UF.u"s")
        quantities = [2UF.u"m", 4UF.u"m", 6UF.u"m"]
        convert_quantity = Base.Fix1(unitless, basis)
        @test unitless.(basis, quantities) == map(convert_quantity, quantities)
        @test unitless.(basis, quantities) == [1, 2, 3]
        matrix = [2UF.u"m" 4UF.u"m"; 6UF.u"m" 8UF.u"m"]
        @test unitless.(basis, matrix) == map(convert_quantity, matrix)
        mixed = (2UF.u"m", 6UF.u"s", 4UF.u"m/s")
        @test unitless.(basis, mixed) == map(convert_quantity, mixed)
        empty_quantities = typeof(2UF.u"m")[]
        @test unitless.(basis, empty_quantities) == map(convert_quantity, empty_quantities)
        @test 2 .* unitless.(basis, quantities .+ 2UF.u"m") .+ 1 == [5, 7, 9]
        destination = zeros(3)
        destination .= unitless.(basis, quantities)
        @test destination == [1, 2, 3]
        @test unitless.(Ref(basis), quantities) == [1, 2, 3]
        @test_throws ArgumentError unitless.(basis, [1UF.u"kg"])
        @test_throws ArgumentError unitless.((1UF.u"m", 2UF.u"m"), quantities)
    end
    @test unitless((), 3) == 3
    @test_throws ArgumentError unitless((), 1UF.u"m")
    @test LessUnit(2UF.u"m")(UF.𝐋) == 2UF.u"m"
end
