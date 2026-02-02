using LessUnits
using Test
using Unitful

@testset "LessUnits.jl" begin
    @test unitof(Unitful.Charge, (2u"q", 1u"ħ", 2π * 1u"GHz")) ≈ 2u"q"
    @test unitless((2u"q", 1u"ħ", 2π * 1u"GHz"), 1u"Φ0") ≈ 2π
    @test unitless((1u"q", 1u"ħ", 2π * 1u"GHz"), 1u"Φ0") ≈ π
    @test unitless((2u"q", 1u"cm"), 5) ≈ 5
    @test unitof(Unitful.Temperature, (2u"q", 1u"ħ", 2π * 1u"GHz", 1u"k")) ≈ 1u"h*GHz/k"
    @test begin
        lu = LessUnit(2π * 1u"GHz", 2u"q", 1u"ħ", 1u"k")
        lu(1u"K") > 0.0
    end
    @test_throws ArgumentError unitof(Unitful.𝐋, (2u"q", 1u"ns"))
    @test_throws ArgumentError unitof(Unitful.𝐓, (1u"GHz", 1u"ns"))
    @test_throws ArgumentError unitof(Unitful.Charge, (2u"q", 1u"fF", 1u"nH", 1u"MHz"))
    @test unitof(u"s", (2π * 1u"GHz", 2u"q")) ≈ inv(2π) * u"ns"
    @test unit(unitof(u"s", (2π * 1u"GHz", 2u"q"))) == u"s"
end
