using Test
using QuantumHamiltonian

@testset "product/decompose" begin
    spinhalfsite = Site([State("Up", +1), State("Dn", -1)])
    spinonesite = Site([State("S=1", +2), State("S=0", 0), State("S=-1", -2)])

    s1 = [spinhalfsite for i in 1:8]
    s2 = [spinonesite for i in 1:4]
    # s2 = s1
    hs = HilbertSpace(vcat(s1, s2))
    
    @testset "no quantum number" begin
        hsr = represent(hs)

        hsr1 = represent(HilbertSpace(s1))
        hsr2 = represent(HilbertSpace(s2))
        phsr = ProductHilbertSpaceRepresentation(hsr1, hsr2)

        dhsr = decompose(hsr)
        dphsr = decompose(phsr)
        
        bl = get_basis_list(hsr)
        @test bl == get_basis_list(phsr)
        let bl2 = get_basis_list(dhsr)
            @test bl != bl2
            @test sort(bl) == sort(bl2)
        end

        let bl2 = get_basis_list(dphsr)
            @test bl != bl2
            @test sort(bl) == sort(bl2)
        end
    end
    @testset "with quantum number" begin
        hsr = represent(hs, [(0,)])
        bl = get_basis_list(hsr)

        hsr1 = represent(HilbertSpace(s1))
        hsr2 = represent(HilbertSpace(s2))
        phsr = ProductHilbertSpaceRepresentation(hsr1, hsr2)

        dhsr = sectorslice(decompose(hsr), [(0,)])
        dphsr = sectorslice(
            decompose(phsr),
            x -> reduce((z,w) -> z .+ w, x) == (0,)
        )

        let bl2 = get_basis_list(dhsr)
            @test bl == bl2
        end
        let bl2 = get_basis_list(dphsr)
            @test bl != bl2
            @test sort(bl) == sort(bl2)
        end

        # @show @allocated represent(hs, [(0,)])
        # @show @allocated let
        #     hsr1 = represent(HilbertSpace(s1))
        #     hsr2 = hsr1
        #     phsr = ProductHilbertSpaceRepresentation(hsr1, hsr2)
        #     dphsr = sectorslice(
        #         decompose(phsr),
        #         x -> reduce((z,w) -> z .+ w, x) == (0,)
        #     )
        # end
    end
end
