using Test
using QuantumHamiltonian
using StaticArrays

@testset "State" begin
    @test isimmutable(State)
    let
        s = State("MyState")
        @test typeof(s) == State{Tuple{}}
        @test s.name == "MyState"
        @test qntype(s) === Tuple{}
        @test qntype(typeof(s)) === Tuple{}
    end
    
    let
        s = State("Up", 1)
        @test typeof(s) == State{Tuple{Int}}
        @test s.name == "Up"
        @test s.quantum_number == (1,)
        @test qntype(s) === Tuple{Int}
        @test qntype(typeof(s)) === Tuple{Int}
    end
    
    let
        s = State("Up", (1, 1))
        @test s.name == "Up"
        @test s.quantum_number == (1, 1)
        @test qntype(s) === Tuple{Int, Int}
        @test qntype(typeof(s)) === Tuple{Int, Int}
    end
end

@testset "Site" begin
    @test isimmutable(Site)
    @testset "spin-half" begin
        up = State("Up", 1)
        dn = State("Dn",-1)
        site1 = Site([up, dn])
        @test qntype(site1) === Tuple{Int}
        @test qntype(typeof(site1)) === Tuple{Int}
        
        @test dimension(site1) == 2
        @test get_state(site1, 0x0000000) == up
        @test get_state(site1, 0x0000001) == dn
        @test quantum_number_sectors(site1) == [(-1,), (1,)]
        @test get_quantum_number(site1, 1) == (1,)
        @test get_quantum_number(site1, 2) == (-1,)
        @test compress(site1, 1) == 0x0000000
        @test compress(site1, 2) == 0x0000001
        @test get_state_index(site1, 0x0000000) == 1
        @test get_state_index(site1, 0x0000001) == 2
        
        @test_throws BoundsError compress(site1, 10)
        @test_throws BoundsError get_state_index(site1, 0x000000F)
        
        @test [i for i in eachindex(site1)] == [1,2]
        @test collect(keys(site1)) == [1,2]
    end
    
    @testset "spin-charge" begin
        QN = Tuple{Int, Int}
        em = State("Em", ( 0, 0))
        up = State("Up", ( 1, 1))
        dn = State("Dn", (-1, 1))
        ud = State("UpDn", ( 0, 2))
        @test qntype(em) === QN
        @test qntype(up) === QN
        @test qntype(dn) === QN
        @test qntype(ud) === QN
        @test_throws MethodError State{QN}("X", 1)
        
        site = Site([em, up, dn])
        @test qntype(site) === QN
        @test qntype(typeof(site)) === QN
        @test_throws MethodError Site([em, up, dn, State("X", 1)])
        @test bitwidth(site) == 2
        @test dimension(site) == 3
        @test get_state(site, 0b0000000) == em
        @test get_state(site, 0b0000001) == up
        @test get_state(site, 0b0000010) == dn
        @test_throws BoundsError get_state(site, 0b0000011) #TODO best exception type?
    end
end

@testset "HilbertSpace" begin
    @test isimmutable(HilbertSpace)
    @testset "spinhalf" begin
        QN = Tuple{Int}
        up = State("Up", 1)
        dn = State("Dn", (-1,))
        spin_site = Site([up, dn])
        hs = HilbertSpace([spin_site, spin_site, spin_site, spin_site])
        hs2 = HilbertSpace{QN}([spin_site, spin_site, spin_site, spin_site])
        @test hs == hs2
        
        @test numsites(hs) == 4
        @test get_site(hs, 2) == spin_site
        @test_throws BoundsError get_site(hs, 100)
        @test qntype(hs) === Tuple{Int}
        @test scalartype(hs) === Bool
        @test valtype(hs) === Bool
        @test qntype(typeof(hs)) === Tuple{Int}
        @test scalartype(typeof(hs)) === Bool
        @test valtype(typeof(hs)) === Bool
        @test basespace(hs) === hs
        
        @test get_bitmask(hs, 1) == 0b0001
        @test get_bitmask(hs, 2) == 0b0010
        @test get_bitmask(hs, 3) == 0b0100
        @test get_bitmask(hs, 4) == 0b1000
        @test_throws BoundsError get_bitmask(hs, 5)
        @test_throws BoundsError get_bitmask(hs,-1)
        
        @test quantum_number_sectors(hs) == [(-4,), (-2,), (0,), (2,), (4,)]
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
        
        @test numsites(hs) == 4
        @test get_site(hs, 2) == site
        @test get_site(hs, 3) == spin_site
        @test_throws BoundsError get_site(hs, 100)
        @test qntype(hs) === QN
        @test scalartype(hs) === Bool
        @test valtype(hs) === Bool
        @test qntype(typeof(hs)) === QN
        @test scalartype(typeof(hs)) === Bool
        @test valtype(typeof(hs)) === Bool
        @test basespace(hs) === hs
                
        @test hs.bitoffsets[end] == 2 + 2 + 1 + 2
        @test bitwidth(hs) == 2 + 2 + 1 + 2
        @test bitwidth(hs, 1) == 2
        @test bitwidth(hs, 2) == 2
        @test bitwidth(hs, 3) == 1
        @test bitwidth(hs, 4) == 2
        @test bitoffset(hs, 1) == 0
        @test bitoffset(hs, 2) == 2
        @test bitoffset(hs, 3) == 4
        @test bitoffset(hs, 4) == 5
        @test get_bitmask(hs) == 0b1111111
        @test get_bitmask(hs, 1) == 0b0000011
        @test get_bitmask(hs, 2) == 0b0001100
        @test get_bitmask(hs, 3) == 0b0010000
        @test get_bitmask(hs, 4) == 0b1100000
        @test_throws BoundsError get_bitmask(hs, 5)
        @test_throws BoundsError get_bitmask(hs,-1)
        
        sectors = [( 1, -1), ( 1,  1),
        ( 2, -2), ( 2,  0), ( 2, 2),
        ( 3, -3), ( 3, -1), ( 3, 1), ( 3, 3),
        ( 4, -4), ( 4, -2), ( 4, 0), ( 4, 2), ( 4, 4)]
        @test quantum_number_sectors(hs) == sectors
        
        @test get_state_index(hs, 0b0000000, 1) == 1
        @test get_state_index(hs, 0b0000001, 1) == 2
        @test get_state_index(hs, 0b0000010, 1) == 3
        @test get_state_index(hs, 0b0000011, 1) == 4 # No exception thrown here
        
        @test get_state_index(hs, 0b0000011, 2) == 1
        @test get_state_index(hs, 0b0000111, 2) == 2
        @test get_state_index(hs, 0b0001011, 2) == 3
        @test get_state_index(hs, 0b0001111, 2) == 4
        
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
        
        for EX in [extract, uncompress]
            @test EX(hs, 0b0000000) == CartesianIndex(1,1,1,1)
            @test EX(hs, 0b0000001) == CartesianIndex(2,1,1,1)
            @test EX(hs, 0b0000010) == CartesianIndex(3,1,1,1)
            @test_throws BoundsError EX(hs, 0b0000011)
            
            @test EX(hs, 0b0010000) == CartesianIndex(1,1,2,1)
            @test EX(hs, 0b0010001) == CartesianIndex(2,1,2,1)
            @test EX(hs, 0b0010010) == CartesianIndex(3,1,2,1)
            @test_throws BoundsError EX(hs, 0b0010011)
        end

        @test compress(hs, CartesianIndex(1,1,2,1)) == 0b0010000
        @test compress(hs, CartesianIndex(2,1,2,1)) == 0b0010001
        @test compress(hs, CartesianIndex(3,1,2,1)) == 0b0010010
        @test_throws BoundsError compress(hs, CartesianIndex(4, 1, 2, 1))
        @test_throws ArgumentError compress(hs, CartesianIndex(1, 2))
    end
end
