using Test
using QuantumHamiltonian

@testset "indexed vector" begin
    elems = ["A", "BB", "CCC"]
    for VT in [DictIndexedVector, SortedIndexedVector]
        v = VT(elems)
        v1 = VT{String}(elems)
        v2 = VT(["A"])
        @test v == v1
        @test v != v2
        @test hash(v) == hash(v1)
        @test hash(v) != hash(v2)
        @test v[1] == "A"
        @test v[2] == "BB"
        @test v[3] == "CCC"
        @test_throws BoundsError v[4]
        @test findindex(v, "A") == 1
        @test findindex(v, "BB") == 2
        @test findindex(v, "CCC") == 3
        @test findindex(v, "DDDD") <= 0
        @test "A" in v
        @test "BB" in v
        @test "CCC" in v
        @test !("DDDD" in v)
        @test [x for x in v] == elems
        @test length(v) == 3
        @test size(v) == (3,)
        @test size(v, 1) == 3
    end
end