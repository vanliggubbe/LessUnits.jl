using LessUnits
using Test
using Unitful

@testset "LessUnits.jl" begin
    @test unitof(Unitful.Charge, (2u"q", 1u"ħ", 2π * 1u"GHz")) ≈ 2u"q"
    @test unitless(1u"Φ0", (2u"q", 1u"ħ", 2π * 1u"GHz")) ≈ 2π
    @test unitless(1u"Φ0", (1u"q", 1u"ħ", 2π * 1u"GHz")) ≈ π
    @test unitless(5, (2u"q", 1u"cm")) ≈ 5
    @test unitof(Unitful.Temperature, (2u"q", 1u"ħ", 2π * 1u"GHz", 1u"k")) ≈ 1u"h*GHz/k"
    @test begin
        lu = LessUnit(2π * 1u"GHz", 2u"q", 1u"ħ", 1u"k")
        lu(1u"K") > 0.0
    end
    @test_throws ArgumentError unitof(Unitful.𝐋, (2u"q", 1u"ns"))
    @test_throws ArgumentError unitof(Unitful.𝐓, (1u"GHz", 1u"ns"))
    @test_throws ArgumentError unitof(Unitful.Charge, (2u"q", 1u"fF", 1u"nH", 1u"MHz"))
end
