using Test
using GroupTools
using QuantumHamiltonian

@testset "LocalGeneralizedPermutation" begin
    gp0 = GeneralizedPermutation([1,2,3], [Phase(0//1), Phase(0//1), Phase(0//1)])
    gp1 = GeneralizedPermutation([1,2], [Phase(0//1), Phase(1//2)])
    gp1a = GeneralizedPermutation([1,2], [Phase(0//1), Phase(1//2)])
    gp2 = GeneralizedPermutation([2,1], [Phase(0//1), Phase(0//1)])
    
    @test !isidentity(LocalGeneralizedPermutation([gp0, gp1]))
    @test !isidentity(LocalGeneralizedPermutation([gp0, gp2]))
    @test isidentity(LocalGeneralizedPermutation([gp0, gp0]))

    lgp1 = LocalGeneralizedPermutation([gp0, gp1])
    lgp1a = LocalGeneralizedPermutation([gp0, gp1a])
    lgp2 = LocalGeneralizedPermutation([gp0, gp2])

    @test lgp1 == lgp1a
    @test lgp1 != lgp2
    
    h = UInt(0x12345)
    @test hash(lgp1, h) == hash(lgp1a, h)
    @test hash(lgp1, h) != hash(lgp2, h)

    @test inv(lgp1) == LocalGeneralizedPermutation([inv(gp0), inv(gp1)])
end