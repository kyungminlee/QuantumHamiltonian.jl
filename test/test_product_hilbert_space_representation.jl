using Test
using QuantumHamiltonian

@testset "ProductHilbertSpaceRepresentation" begin
    spinhalfsite = Site([State("Up", +1), State("Dn", -1)])
    spinonesite = Site([State("S=1", +2), State("S=0", 0), State("S=-1", -2)])

    @testset "trivial case" begin
        bighs = HilbertSpace([spinhalfsite, spinhalfsite, spinhalfsite, spinonesite, spinonesite])
        bighsr = represent(bighs)

        hs1 = HilbertSpace([spinhalfsite, spinhalfsite, spinhalfsite])
        hs2 = HilbertSpace([spinonesite, spinonesite])
        phs = ProductHilbertSpace((hs1, hs2))

        hsr1 = represent(hs1)
        hsr2 = represent(hs2)
        phsr = ProductHilbertSpaceRepresentation(hsr1, hsr2)

        @testset "from hilbert space" begin
            @test valtype(phsr) == Bool
            @test scalartype(phsr) == Bool
            @test bintype(phsr) == UInt
            @test qntype(phsr) == Tuple{Int}
            @test qntype(typeof(phsr)) == Tuple{Int}
            @test tagtype(phsr, Val(:QuantumNumberAsTag)) == Tuple{Tuple{Int}, Tuple{Int}}
            @test basespace(phsr) == phs
            @test basespace(typeof(phsr)) == typeof(basespace(phsr))
            @test numsites(phsr) == 5
            @test_throws BoundsError get_site(phsr, 0)
            @test get_site(phsr, 1) == spinhalfsite
            @test get_site(phsr, 4) == spinonesite
            @test_throws BoundsError get_site(phsr, 100)
            @test bitwidth(phsr) == 7
            @test_throws BoundsError bitwidth(phsr, 0)
            @test bitwidth(phsr, 2) == 1
            @test bitwidth(phsr, 4) == 2
            @test bitwidth(phsr, 100) == 0
            @test_throws BoundsError bitoffset(phsr, 0)
            @test bitoffset(phsr, 1) == 0
            @test bitoffset(phsr, 2) == 1
            @test bitoffset(phsr, 3) == 2
            @test bitoffset(phsr, 4) == 3
            @test bitoffset(phsr, 5) == 5
            @test bitoffset(phsr, 6) == 7
            @test bitoffset(phsr, 100) == 7
            @test get_bitmask(phsr, 2) == 0b0000010
            @test get_bitmask(phsr, 4) == 0b0011000
            # @test get_quantum_numbers(phsr) == [(x,) for x in -7:2:7]
            @test get_quantum_number(phsr, 0b00_00_0_0_0) == (7,)
            @test get_quantum_number(phsr, 0b00_10_1_0_0) == (1,)
            # @test get_tags(phsr) == [((x,), (y,)) for y in -4:2:4 for x in -3:2:3]
            @test get_tag(phsr, 0b00_00_0_0_0, Val(:QuantumNumberAsTag)) == ((3,), (4,))
            @test get_tag(phsr, 0b00_10_1_0_0, Val(:QuantumNumberAsTag)) == ((1,), (0,))
            @test extract(phsr, 0b00_00_0_0_0) == CartesianIndex(1,1,1,1,1)
            @test extract(phsr, 0b00_10_1_0_0) == CartesianIndex(1,1,2,3,1)
            @test compress(phsr, CartesianIndex(1,1,1,1,1)) == 0b00_00_0_0_0
            @test compress(phsr, CartesianIndex(1,1,2,3,1)) == 0b00_10_1_0_0
            @test update(phsr, 0b00_10_1_0_0, 2, 2) == 0b00_10_1_1_0
            @test_throws BoundsError update(phsr, 0b00_10_1_0_0, 100, 2)
            @test_throws BoundsError update(phsr, 0b00_10_1_0_0, 2, 100)
            @test get_state_index(phsr, 0b00_10_1_0_0, 4) == 3
            @test get_state(phsr, 0b00_10_1_0_0, 4) ==  State("S=-1", -2)
        end

        @test dimension(phsr) == 3*3*2*2*2

        let bl = get_basis_list(phsr)
            @test length(bl) == 3*3*2*2*2
            @test bl == get_basis_list(bighsr)
            for (i, b) in enumerate(bl)
                @test get_basis_index_amplitude(phsr, b) == (index=i, amplitude=1)
                @test get_basis_state(phsr, i) == b
            end
        end
        @test get_basis_index_amplitude(phsr, 0b11_10_1_1_1) == (index=-1, amplitude=0)
        @test get_basis_index_amplitude(phsr, 0xFFFFFFFFFFFFFF) == (index=-1, amplitude=0)
        @test_throws BoundsError get_basis_state(phsr, -10)
    end

    @testset "nontrivial" begin
        hs1 = HilbertSpace([spinhalfsite, spinhalfsite, spinhalfsite])
        hs2 = HilbertSpace([spinonesite, spinonesite])
        hss1 = HilbertSpaceSector(hs1, 0)
        hss2 = HilbertSpaceSector(hs2, 1)
        hsr1 = represent(hss1)
        hsr2 = represent(hss2)

        phsr = ProductHilbertSpaceRepresentation(hsr1, hsr2)
        @test dimension(phsr) == dimension(hsr1) * dimension(hsr2)
        
        let bl = get_basis_list(phsr)
            @test length(bl) == dimension(phsr)
            bl1 = get_basis_list(hsr1)
            bl2 = get_basis_list(hsr2)
            bw = (bitwidth(hs1), bitwidth(hs2))
            blp = [bitjoin(bw, (b1, b2)) for b2 in bl2 for b1 in bl1]
            @test bl == blp
            for (i, b) in enumerate(bl)
                @test get_basis_index_amplitude(phsr, b) == (index=i, amplitude=1)
                @test get_basis_state(phsr, i) == b
            end
        end
        @test get_basis_index_amplitude(phsr, 0b10_10_0_0_0) == (index=-1, amplitude=0)
        @test get_basis_index_amplitude(phsr, 0xFFFFFFFFFFFFFF) == (index=-1, amplitude=0)
        @test_throws BoundsError get_basis_state(phsr, -10)
    end
end