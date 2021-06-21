using Test
using QuantumHamiltonian

@testset "util" begin

    @testset "make_bitmask" begin
        for bintype in [UInt8, UInt16, UInt32, UInt64, UInt128]
            @test isa(make_bitmask(4, bintype), bintype)
            @test isa(make_bitmask(4, 2, bintype), bintype)
        end
        @test make_bitmask(5) == 0b11111
        @test make_bitmask(5, 3) == 0b11000
    end

    @testset "bitsplit/bitjoin" begin
        bvec = UInt(0b101_1011_010)
        bw = (3,4,3)
        bvecs = bitsplit(bw, bvec)
        bvec2 = bitjoin(bw, bvecs)
        @test bvecs == (UInt(0b010), UInt(0b1011), UInt(0b101))
        @test bvec == bvec2
        @test bitjoin(bw, (UInt(0b1111111), UInt(0b0), UInt(0b1111111))) == UInt(0b1110000111) # test for bit overflow
    end

    @testset "compress" begin
        bitwidth = [2, 3, 1, 2]
        @test compress(bitwidth, [1,2,0,1], UInt) == UInt(0b01_0_010_01)
        @test_throws ArgumentError compress([1,2,3], [1,2,3,4], UInt)
        @test_throws ArgumentError compress([8,8,8,8], [1,1,1,1], UInt8)
        @test_throws ArgumentError compress([1,1,1], [99, 1,1], UInt)
    end

    @testset "merge_vec" begin
        let
            a = Int[1,5,7,8]
            b = Int[2,3,4,6,9]
            c = QuantumHamiltonian.merge_vec(a, b)
            @test c == [1,2,3,4,5,6,7,8,9]
        end

        let
            a = Int[1,2,3]
            b = Int[4,5,6,7]
            c = QuantumHamiltonian.merge_vec(a, b)
            @test c == [1,2,3,4,5,6,7]
        end

        let
            a = Int[4,5,6,7]
            b = Int[1,2,3]
            c = QuantumHamiltonian.merge_vec(a, b)
            @test c == [1,2,3,4,5,6,7]
        end

        let
            a = Int[1,2,3]
            b = Int[4,5,6,7]
            c = QuantumHamiltonian.merge_vec(a, b)
            @test c == [1,2,3,4,5,6,7]
        end

        let
            a = Int[1,2,3,4]
            b = Int[5,6,7]
            c = QuantumHamiltonian.merge_vec(a, b)
            @test c == [1,2,3,4,5,6,7]
        end

        let
            a = Int[1,2,3]
            b = Int[4,5,6,7]
            c = QuantumHamiltonian.merge_vec(a, b)
            @test c == [1,2,3,4,5,6,7]
        end

        let
            a = Int[1,1,1,1]
            b = Int[1,1,1,1,1]
            c = QuantumHamiltonian.merge_vec(a, b)
            @test c == [1,1,1,1,1,1,1,1,1]
        end

    end

    @testset "choptol!" begin
        let
            d = Dict("A" => 0.0, "B" => 1E-9)
            choptol!(d, 1E-12)
            @test d == Dict("B"=>1E-9)
        end
        let
            d = Dict("A" => 0.0, "B" => 1E-9)
            choptol!(d, 1E-6)
            @test d == Dict()
        end
    end

    @testset "IntegerModulo" begin
        T = IntegerModulo{3}
        t0 = T(0)
        t1 = T(1)
        t2 = T(2)
        t3 = T(3)
        @test t0 != t1
        @test t0 == t3
        @test t0 == 0
        @test t3 == 0
        @test t1 == 1
        @test IntegerModulo{3}(0) != IntegerModulo{4}(0)

        @test t1 + t2 == 1 + t2 == t1 + 2 == t0 == 0 == +t0 == -t3
        @test t1 - t2 == 1 - t2 == t1 - 2 == t2 == 2 == +t2 == -t1
        @test t2 * t2 == 2 * t2 == t2 * 2 == t1 == 1 == +t1 == -t2
        @test t3 * t2 == 3 * t2 == t3 * 2 == t0 == 0 == +t0 == -t3
    end

    @testset "tuple" begin
        QH = QuantumHamiltonian
        t1 = (1.0, 2, 3.0 + 4im)
        T1 = typeof(t1)
        @test QH.tupleone(T1)  === (1.0, 1, 1.0 + 0.0im)
        @test QH.tupleone(T1)  ==  (1.0, 1, 1.0 + 0.0im)
        @test QH.tuplezero(T1) === (0.0, 0, 0.0 + 0.0im)
        @test QH.tuplezero(T1) ==  (0.0, 0, 0.0 + 0.0im)
        @test QH.tupleone(T1)  !== (1, 1, 1)
        @test QH.tupleone(T1)  ==  (1, 1, 1)
        @test QH.tuplezero(T1) !== (0, 0, 0)
        @test QH.tuplezero(T1) ==  (0, 0, 0)

        @test QH.tupleone(t1)  === (1.0, 1, 1.0 + 0.0im)
        @test QH.tupleone(t1)  ==  (1.0, 1, 1.0 + 0.0im)
        @test QH.tuplezero(t1) === (0.0, 0, 0.0 + 0.0im)
        @test QH.tuplezero(t1) ==  (0.0, 0, 0.0 + 0.0im)
        @test QH.tupleone(t1)  !== (1, 1, 1)
        @test QH.tupleone(t1)  ==  (1, 1, 1)
        @test QH.tuplezero(t1) !== (0, 0, 0)
        @test QH.tuplezero(t1) ==  (0, 0, 0)

        @test QH.tupleadd((1.0, 2, 3.0 + 4im), (5.0, 6, 7.0 + 8.0im)) === (6.0, 8, 10.0 + 12.0im)
        @test QH.tuplesubtract((1.0, 2, 3.0 + 4im), (5.0, 6, 7.0 + 8.0im)) === (-4.0, -4, -4.0 - 4.0im)

        @test QH.tuplelength((1, 'B', "γ")) == 3
        @test QH.tuplelength(typeof((1,'B', "γ"))) == 3
    end
end


@testset "FrozenSortedArray" begin
    @test_throws ArgumentError FrozenSortedArrayIndex{UInt}(UInt[0x2, 0x1, 0x3])
    @test_throws ArgumentError FrozenSortedArrayIndex{UInt}(UInt[0x1, 0x1, 0x3])
    fsa = FrozenSortedArrayIndex{UInt}(UInt[0x2, 0x4, 0x6])
    @test fsa[0x2] == 1
    @test fsa[0x4] == 2
    @test fsa[0x6] == 3
    @test collect(fsa) == [0x2 =>1, 0x4=>2, 0x6=>3]
    @test collect(k=>v for (k, v) in fsa) == [0x2 =>1, 0x4=>2, 0x6=>3]
    @test length(fsa) == 3
    @test_throws KeyError fsa[0x1]
    @test haskey(fsa, 0x2)
    @test !haskey(fsa, 0x3)
    @test get(fsa, 0x1, -1) == -1
    @test keys(fsa) == UInt[0x2, 0x4, 0x6]

    @test eltype(fsa) === Pair{UInt, Int}
    @test eltype(typeof(fsa)) === Pair{UInt, Int}
    @test [k for (k, v) in fsa] == UInt[0x2, 0x4, 0x6]
    @test [v for (k, v) in fsa] == UInt[0x1, 0x2, 0x3]
end
