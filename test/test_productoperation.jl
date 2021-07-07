using Test
using QuantumHamiltonian
using GroupTools
using LatticeTools

@testset "ProductOperation" begin
    site = Site([State("up", 0), State("dn", 1)])
    nsites = 4
    hs = HilbertSpace([site for i in 1:nsites])
    hsr = represent(hs, [(1,)])

    ϕ = Phase(1//3)
    s = SitePermutation([2,3,4,1])
    p = ProductOperation(ϕ, s)
    q = ProductOperation(s, ϕ)
    for b in hsr.basis
        let
            b1, a1 = b, 1
            b2, a2 = symmetry_apply(hs, ϕ, b1, a1)
            @test b2 == b1
            @test isapprox(a2, cis(2π/3))

            b3, a3 = symmetry_apply(hs, s, b2, a2)

            b4, a4 = symmetry_apply(hs, p, b1, a1)
            @test b3 == b4
            @test isapprox(a3, a4)
        end
    end
end

@testset "DirectProductOperation" begin
    site = Site([State("up", 0), State("dn", 1)])
    nsites = 4
    hs = HilbertSpace([site for i in 1:nsites])
    hsr = represent(hs, [(1,)])

    ϕ = Phase(1//3)
    s = SitePermutation([2,3,4,1])
    p = DirectProductOperation(ϕ, s)
    q = DirectProductOperation(s, ϕ)
    for b in hsr.basis
        let
            b1, a1 = b, 1
            b2, a2 = symmetry_apply(hs, ϕ, b1, a1)
            @test b2 == b1
            @test isapprox(a2, cis(2π/3))

            b3, a3 = symmetry_apply(hs, s, b2, a2)

            b4, a4 = symmetry_apply(hs, p, b1, a1)
            @test b3 == b4
            @test isapprox(a3, a4)

            b5, a5 = symmetry_apply(hs, q, b1, a1)
            @test b3 == b5
            @test isapprox(a3, a5)
        end
    end
end