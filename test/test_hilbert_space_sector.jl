using Test
using QuantumHamiltonian
using StaticArrays

@testset "HilbertSpaceSector" begin
    @test isimmutable(HilbertSpaceSector)
    @testset "spinhalf" begin
        QN = Tuple{Int}
        up = State("Up", 1)
        dn = State("Dn",(-1,))
        spin_site = Site([up, dn])
        hs = let
            hs = HilbertSpace([spin_site, spin_site, spin_site, spin_site])
            @test HilbertSpaceSector(hs).allowed_quantum_numbers == Set([(-4,), (-2,), (0,), (2,), (4,)])
            @test HilbertSpaceSector(hs, 0).allowed_quantum_numbers == Set([(0,)])
            @test HilbertSpaceSector(hs, [(2,), (4,)]).allowed_quantum_numbers == Set([(2,), (4,)])
            @test HilbertSpaceSector(hs, [(1,)]).allowed_quantum_numbers == Set(Tuple{Int}[])
            
            @test HilbertSpaceSector(hs) == HilbertSpaceSector(hs, [(-4,), (-2,), (0,), (2,), (4,)])
            @test HilbertSpaceSector(hs) != HilbertSpaceSector(hs, (0,))
            
            @test HilbertSpaceSector(hs, [0, 2]) == HilbertSpaceSector(hs, [(0,), (2,)])
            @test HilbertSpaceSector(hs, Set([0, 2])) == HilbertSpaceSector(hs, [(0,), (2,)])
            
            hss = HilbertSpaceSector(hs, (0,))
            
            @test basespace(hss) != hss
            @test basespace(hss) == hs
            hss
        end
        @test qntype(hs) === Tuple{Int}
        @test qntype(typeof(hs)) === Tuple{Int}
        @test get_quantum_numbers(hs) == [(0,)]

        @test tagtype(hs, Val(:QuantumNumberAsTag)) === Tuple{Int}
        @test tagtype(typeof(hs), Val(:QuantumNumberAsTag)) === Tuple{Int}
        @test get_tags(hs, Val(:QuantumNumberAsTag)) == [(0,)]
    
        @test bitwidth(hs) == 4
        
        @test get_bitmask(hs, 1) == 0b0001
        @test get_bitmask(hs, 2) == 0b0010
        @test get_bitmask(hs, 3) == 0b0100
        @test get_bitmask(hs, 4) == 0b1000
        @test_throws BoundsError get_bitmask(hs, 5)
        @test_throws BoundsError get_bitmask(hs,-1)
        
        @test quantum_number_sectors(hs) == [(0,)]
        @test get_quantum_number(hs, 0b0000) == (4,)
        @test get_quantum_number(hs, 0b0001) == (2,)
        @test get_quantum_number(hs, 0b0010) == (2,)
        @test get_quantum_number(hs, 0b0011) == (0,)
        @test get_quantum_number(hs, 0b0100) == (2,)
        @test get_quantum_number(hs, 0b0101) == (0,)
        @test get_quantum_number(hs, 0b0110) == (0,)
        @test get_quantum_number(hs, 0b0111) == (-2,)
        
        @test get_quantum_number(hs, [1,1,1,1]) == (4,)
        @test get_quantum_number(hs, [2,1,1,1]) == (2,)
        @test get_quantum_number(hs, [1,2,1,1]) == (2,)
        @test get_quantum_number(hs, [2,2,1,1]) == (0,)
        @test get_quantum_number(hs, [1,1,2,1]) == (2,)
        @test get_quantum_number(hs, [2,1,2,1]) == (0,)
        @test get_quantum_number(hs, [1,2,2,1]) == (0,)
        @test get_quantum_number(hs, [2,2,2,1]) == (-2,)
    end
    
    @testset "charge-spin" begin
        QN = Tuple{Int, Int}
        em = State("Em", ( 0, 0))  # charge and spin
        up = State("Up", ( 1, 1))
        dn = State("Dn", ( 1,-1))
        spin_site = Site([up, dn])
        site = Site([em, up, dn])
        hs = HilbertSpace([site, site, spin_site, site])
        
        @test qntype(hs) === QN        
        @test qntype(typeof(hs)) === QN

        @test basespace(hs) === hs
        
        @test hs.bitoffsets[end] == 2 + 2 + 1 + 2
        @test bitwidth(hs) == 2 + 2 + 1 + 2
        @test get_bitmask(hs, 1) == 0b0000011
        @test get_bitmask(hs, 2) == 0b0001100
        @test get_bitmask(hs, 3) == 0b0010000
        @test get_bitmask(hs, 4) == 0b1100000
        @test_throws BoundsError get_bitmask(hs, 5)
        @test_throws BoundsError get_bitmask(hs,-1)
        
        sectors = [(1, -1), (1, 1),
        (2, -2), (2, 0), (2, 2),
        (3, -3), (3,-1), (3, 1), (3, 3),
        (4, -4), (4,-2), (4, 0), (4, 2), (4, 4)]
        @test quantum_number_sectors(hs) == sectors
        
        @test get_state_index(hs, 0b0000000, 1) == 1
        @test get_state_index(hs, 0b0000001, 1) == 2
        @test get_state_index(hs, 0b0000010, 1) == 3
        @test get_state_index(hs, 0b0000011, 1) == 4 # no exception thrown here because of specialization
        
        @test get_state_index(hs, 0b0000011, 2) == 1
        @test get_state_index(hs, 0b0000111, 2) == 2
        @test get_state_index(hs, 0b0001011, 2) == 3
        @test get_state_index(hs, 0b0001111, 2) == 4 # no exception thrown here because of specialization
        
        @test get_state_index(hs, 0b0000000, 3) == 1
        @test get_state_index(hs, 0b0010000, 3) == 2
        @test get_state_index(hs, 0b1101111, 3) == 1
        @test get_state_index(hs, 0b1111111, 3) == 2
        
        @test get_state(hs, 0b0000000, 1) == em
        @test get_state(hs, 0b0000001, 1) == up
        @test get_state(hs, 0b0000010, 1) == dn
        @test_throws BoundsError get_state(hs, 0b0000011, 1)
        
        @test get_state(hs, 0b0000011, 2) == em
        @test get_state(hs, 0b0000111, 2) == up
        @test get_state(hs, 0b0001011, 2) == dn
        @test_throws BoundsError get_state(hs, 0b0001111, 2)
        
        @test get_state(hs, 0b0000000, 3) == up
        @test get_state(hs, 0b0010000, 3) == dn
        @test get_state(hs, 0b1101111, 3) == up
        @test get_state(hs, 0b1111111, 3) == dn
        
        
        @test get_quantum_number(hs, 0b0000000) == (1, 1)
        @test get_quantum_number(hs, 0b0000001) == (2, 2)
        @test get_quantum_number(hs, 0b0000010) == (2, 0)
        @test_throws BoundsError get_quantum_number(hs, 0b0000011)
        @test get_quantum_number(hs, 0b0000100) == (2, 2) # em - up - up - em
        @test get_quantum_number(hs, 0b0000101) == (3, 3) # up - up - up - em
        @test get_quantum_number(hs, 0b0000110) == (3, 1) # dn - up - up - em
        @test_throws BoundsError get_quantum_number(hs, 0b0000111)
        
        @test update(hs, 0b1111110, 1, 1) == 0b1111100
        @test update(hs, 0b1111110, 1, 2) == 0b1111101
        @test update(hs, 0b1111110, 1, 3) == 0b1111110
        @test_throws BoundsError update(hs, 0b1111110, 1, 4)
        
        @test extract(hs, 0b0000000) == CartesianIndex(1,1,1,1)
        @test extract(hs, 0b0000001) == CartesianIndex(2,1,1,1)
        @test extract(hs, 0b0000010) == CartesianIndex(3,1,1,1)
        @test_throws BoundsError extract(hs, 0b0000011)
        
        @test extract(hs, 0b0010000) == CartesianIndex(1,1,2,1)
        @test extract(hs, 0b0010001) == CartesianIndex(2,1,2,1)
        @test extract(hs, 0b0010010) == CartesianIndex(3,1,2,1)
        @test_throws BoundsError extract(hs, 0b0010011)
        
        @test compress(hs, CartesianIndex(1,1,2,1)) == 0b0010000
        @test compress(hs, CartesianIndex(2,1,2,1)) == 0b0010001
        @test compress(hs, CartesianIndex(3,1,2,1)) == 0b0010010
        @test_throws BoundsError compress(hs, CartesianIndex(4, 1, 2, 1))
        @test_throws ArgumentError compress(hs, CartesianIndex(1, 2))
    end
end
