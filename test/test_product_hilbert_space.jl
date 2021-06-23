using Base: Bool
using Test
using QuantumHamiltonian

@testset "ProductHilbertSpace" begin
    spinhalfsite = Site([State("Up", +1), State("Dn", -1)])
    spinonesite = Site([State("S=1", +2), State("S=0", 0), State("S=-1", -2)])
    hs1 = HilbertSpace([spinhalfsite, spinhalfsite, spinhalfsite])
    hs2 = HilbertSpace([spinonesite, spinonesite])
    phs = ProductHilbertSpace((hs1, hs2))
    @test qntype(phs) == Tuple{Int}
    @test basespace(phs) === phs
    @test numsites(phs) == 5
    @test get_site(phs, 1) == spinhalfsite
    @test get_site(phs, 4) == spinonesite
    @test_throws BoundsError get_site(phs, 100)
    @test bitwidth(phs) == 7
    @test bitwidth(phs, 2) == 1
    @test bitwidth(phs, 4) == 2
    @test get_bitmask(phs, 2) == 0b0000010
    @test get_bitmask(phs, 4) == 0b0011000
    @test get_quantum_numbers(phs) == [(x,) for x in -7:2:7]
    @test get_quantum_number(phs, 0b00_00_0_0_0) == (7,)
    @test get_quantum_number(phs, 0b00_10_1_0_0) == (1,)
    @test get_tags(phs) == [((x,), (y,)) for y in -4:2:4 for x in -3:2:3]
    @test get_tag(phs, 0b00_00_0_0_0) == ((3,), (4,))
    @test get_tag(phs, 0b00_10_1_0_0) == ((1,), (0,))
    @test extract(phs, 0b00_00_0_0_0) == CartesianIndex(1,1,1,1,1)
    @test extract(phs, 0b00_10_1_0_0) == CartesianIndex(1,1,2,3,1)
    @test compress(phs, CartesianIndex(1,1,1,1,1)) == 0b00_00_0_0_0
    @test compress(phs, CartesianIndex(1,1,2,3,1)) == 0b00_10_1_0_0
    @test update(phs, 0b00_10_1_0_0, 2, 2) == 0b00_10_1_1_0
    @test_throws BoundsError update(phs, 0b00_10_1_0_0, 100, 2)
    @test_throws BoundsError update(phs, 0b00_10_1_0_0, 2, 100)
    @test get_state_index(phs, 0b00_10_1_0_0, 4) == 3
    @test get_state(phs, 0b00_10_1_0_0, 4) ==  State("S=-1", -2)
end