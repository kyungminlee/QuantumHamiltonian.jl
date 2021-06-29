using Test
using QuantumHamiltonian

@testset "DecomposedHilbertSpaceRepresentation" begin
    spinhalfsite = Site([State("Up", +1), State("Dn", -1)])
    spinonesite = Site([State("S=1", +2), State("S=0", 0), State("S=-1", -2)])

    @testset "trivial case" begin
        hs = HilbertSpace([spinhalfsite, spinhalfsite, spinhalfsite, spinonesite, spinonesite])
        hsr = represent(hs)
        dhsr = DecomposedHilbertSpaceRepresentation(hsr)

        @testset "from hilbert space" begin
            @test qntype(dhsr) == Tuple{Int}
            @test tagtype(dhsr) == Tuple{Int}
            @test basespace(dhsr) == hs
            @test numsites(dhsr) == 5
            @test_throws BoundsError get_site(dhsr, 0)
            @test get_site(dhsr, 1) == spinhalfsite
            @test get_site(dhsr, 4) == spinonesite
            @test_throws BoundsError get_site(dhsr, 100)
            @test bitwidth(dhsr) == 7
            @test_throws BoundsError bitwidth(dhsr, 0)
            @test bitwidth(dhsr, 2) == 1
            @test bitwidth(dhsr, 4) == 2
            # @test bitwidth(phsr, 100) == 0
            @test_throws BoundsError bitoffset(dhsr, 0)
            @test bitoffset(dhsr, 1) == 0
            @test bitoffset(dhsr, 2) == 1
            @test bitoffset(dhsr, 3) == 2
            @test bitoffset(dhsr, 4) == 3
            @test bitoffset(dhsr, 5) == 5
            @test bitoffset(dhsr, 6) == 7
            # @test bitoffset(phsr, 100) == 7
            @test get_bitmask(dhsr, 2) == 0b0000010
            @test get_bitmask(dhsr, 4) == 0b0011000
            # @test get_quantum_numbers(phsr) == [(x,) for x in -7:2:7]
            @test get_quantum_number(dhsr, 0b00_00_0_0_0) == (7,)
            @test get_quantum_number(dhsr, 0b00_10_1_0_0) == (1,)
            # @test get_tags(phsr) == [((x,), (y,)) for y in -4:2:4 for x in -3:2:3]
            @test get_tag(dhsr, 0b00_00_0_0_0) == (7,)
            @test get_tag(dhsr, 0b00_10_1_0_0) == (1,)
            @test extract(dhsr, 0b00_00_0_0_0) == CartesianIndex(1,1,1,1,1)
            @test extract(dhsr, 0b00_10_1_0_0) == CartesianIndex(1,1,2,3,1)
            @test compress(dhsr, CartesianIndex(1,1,1,1,1)) == 0b00_00_0_0_0
            @test compress(dhsr, CartesianIndex(1,1,2,3,1)) == 0b00_10_1_0_0
            @test update(dhsr, 0b00_10_1_0_0, 2, 2) == 0b00_10_1_1_0
            @test_throws BoundsError update(dhsr, 0b00_10_1_0_0, 100, 2)
            @test_throws BoundsError update(dhsr, 0b00_10_1_0_0, 2, 100)
            @test get_state_index(dhsr, 0b00_10_1_0_0, 4) == 3
            @test get_state(dhsr, 0b00_10_1_0_0, 4) ==  State("S=-1", -2)
        end

        @test dimension(dhsr) == 3*3*2*2*2

        let bl = get_basis_list(dhsr)
            @test length(bl) == 3*3*2*2*2
            @test Set(bl) == Set(get_basis_list(hsr))
            for (i, b) in enumerate(bl)
                @test get_basis_index_amplitude(dhsr, b) == (index=i, amplitude=1)
                @test get_basis_state(dhsr, i) == b
            end
        end
        @test_throws BoundsError get_basis_index_amplitude(dhsr, 0b11_10_1_1_1)
        @test_throws BoundsError get_basis_index_amplitude(dhsr, 0xFFFFFFFFFFFFFF)
        @test_throws BoundsError get_basis_state(dhsr, -10)
    end

    @testset "nontrivial" begin
        hs = HilbertSpace([spinhalfsite, spinhalfsite, spinhalfsite, spinonesite, spinonesite])
        hss = HilbertSpaceSector(hs, 0)
        hsr = represent(hss)

        dhsr = DecomposedHilbertSpaceRepresentation(hsr)
        @test dimension(dhsr) == dimension(hsr)

        let bl1 = get_basis_list(hsr),
            bl2 = get_basis_list(dhsr)
            @test length(bl1) == dimension(dhsr)
            @test length(bl2) == dimension(dhsr)
            for (i, b) in enumerate(bl1)
                @test get_basis_index_amplitude(dhsr, b) == (index=i, amplitude=1)
                @test get_basis_state(dhsr, i) == b
                @test get_tag(dhsr, b) == (0,)
            end
        end
        @test_throws BoundsError get_basis_index_amplitude(dhsr, 0b11_10_1_1_1)
        @test_throws BoundsError get_basis_index_amplitude(dhsr, 0xFFFFFFFFFFFFFF)
        @test_throws BoundsError get_basis_state(dhsr, -10)
    end
end