using Test
using QuantumHamiltonian

@testset "indexed vector" begin
    elems = ["A", "BB", "CCC"]
    for VT in [DictIndexedVector, SortedIndexedVector]
        v = VT(elems)
        let
            v1 = VT{String}(elems)
            @test v == v1
            @test hash(v) == hash(v1)
            v2 = VT(["A"])
            @test v != v2
            @test hash(v) != hash(v2)
            if VT == DictIndexedVector
                lookup = Dict(v => i for (i, v) in enumerate(elems))
                v3 = VT{String}(elems, lookup)
                @test v1 == v3

                lookup2 = Dict("A" => 1)
                v4 = VT{String}(elems, lookup2)
                @test_throws KeyError QuantumHamiltonian.checkvalid(v4)
            end
        end
        @test eltype(v) == String
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
        @test collect(v) == elems
        @test length(v) == 3
        @test size(v) == (3,)
        @test size(v, 1) == 3
        @test (QuantumHamiltonian.checkvalid(v); true)
    end
end