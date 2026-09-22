using LessUnits
using Test
import Unitful as UF
import DynamicQuantities as DQ

@testset "LessUnits.jl" begin
    include("unitful.jl")
    include("dynamicquantities.jl")
    include("mixed.jl")
    @test isempty(Test.detect_ambiguities(LessUnits; recursive=true))
end
