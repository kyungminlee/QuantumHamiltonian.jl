using Test
using QuantumHamiltonian
using GroupTools
using LatticeTools

@testset "SemidirectProductOperation" begin
    site = Site([State("up", 0), State("dn", 1)])
    nsites = 4
    hs = HilbertSpace([site for i in 1:nsites])
    hss = HilbertSpaceSector(hs, 1)
    hsr = represent(hss)

    p = SitePermutation([2,3,4,1])
    ϕ = LocalGeneralizedPermutation([
        GeneralizedPermutation([1,2], [Phase(0//1), Phase(1//7)]),
        GeneralizedPermutation([1,2], [Phase(0//1), Phase(2//7)]),
        GeneralizedPermutation([1,2], [Phase(0//1), Phase(3//7)]),
        GeneralizedPermutation([1,2], [Phase(0//1), Phase(4//7)]),
    ])

    q = SitePermutation([3,1,2,4])
    ψ = LocalGeneralizedPermutation([
        GeneralizedPermutation([1,2], [Phase(0//1), Phase(1//11)]),
        GeneralizedPermutation([1,2], [Phase(0//1), Phase(2//11)]),
        GeneralizedPermutation([1,2], [Phase(0//1), Phase(3//11)]),
        GeneralizedPermutation([1,2], [Phase(0//1), Phase(4//11)]),
    ])

    pϕ = SemidirectProductOperation(p, ϕ)
    qψ = SemidirectProductOperation(q, ψ)
    ipϕ = inv(pϕ)
    pϕqψ = pϕ * qψ
    for b in hsr.basis_list
        let
            b1, a1 = b, 1
            b2, a2 = symmetry_apply(hs, ϕ, b1, a1)
            b3, a3 = symmetry_apply(hs, p, b2, a2)

            bx, ax = symmetry_apply(hs, pϕ, b1, a1)

            @test b3==bx
            @test isapprox(a3, ax)
        end
        let
            b1, a1 = b, 1
            b2, a2 = symmetry_apply(hs, inv(p), b1, a1)
            b3, a3 = symmetry_apply(hs, inv(ϕ), b2, a2)

            bx, ax = symmetry_apply(hs, ipϕ, b1, a1)
            @test b3==bx
            @test isapprox(a3, ax)
        end
        let
            b1, a1 = b, 1
            b2, a2 = symmetry_apply(hs, qψ, b1, a1)
            b3, a3 = symmetry_apply(hs, pϕ, b2, a2)

            bx, ax = symmetry_apply(hs, pϕqψ, b1, a1)

            @test b3==bx
            @test isapprox(a3, ax)
        end
    end
end