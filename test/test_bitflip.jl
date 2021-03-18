using Test
using QuantumHamiltonian

@testset "GlobalBitFlip" begin
    b0 = GlobalBitFlip()
    b1 = GlobalBitFlip(true)
    
    @test b0 == GlobalBitFlip(false)
    @test b0 != GlobalBitFlip(true)
    @test b1 != GlobalBitFlip(false)
    @test b1 == GlobalBitFlip(true)
    
    @test b0 * b0 == b0
    @test b0 * b1 == b1
    @test b1 * b0 == b1
    @test b1 * b1 == b0

    @test isidentity(b0)
    @test !isidentity(b1)

    @testset "HilbertSpace" begin
        hs, pauli = QuantumHamiltonian.Toolkit.spin_half_system(4)
        p = SitePermutation([2,3,4,1])
        pb0 = p * b0
        pb1 = p * b1
        b0p = b0 * p
        b1p = b1 * p
        @test pb0 == p * b0
        @test pb1 == p * b1
        @test b0p == b0 * p
        @test b1p == b1 * p

        let bvec = 0b0001
            @test symmetry_apply(hs, b0, 0b0001)[1] == 0b0001
            @test symmetry_apply(hs, b1, 0b0001)[1] == 0b1110
            @test symmetry_apply(hs, pb0, 0b0001)[1] == 0b0010
            @test symmetry_apply(hs, pb1, 0b0001)[1] == 0b1101
        end

        hsr = represent(hs)
        for bvec in hsr.basis_list
            bvec_p, sgn_p = symmetry_apply(hs, p, bvec)
            bvec_p_b0, sgn_p_b0 = symmetry_apply(hs, b0, bvec_p)
            bvec_p_b1, sgn_p_b1 = symmetry_apply(hs, b1, bvec_p)
            bvec_pb0, sgn_pb0 = symmetry_apply(hs, pb0, bvec)
            bvec_pb1, sgn_pb1 = symmetry_apply(hs, pb1, bvec)
            bvec_b0p, sgn_b0p = symmetry_apply(hs, pb0, bvec)
            bvec_b1p, sgn_b1p = symmetry_apply(hs, pb1, bvec)
            @test bvec_p_b0 == bvec_pb0
            @test bvec_p_b1 == bvec_pb1
            @test sgn_p_b0 == sgn_pb0
            @test sgn_p_b1 == sgn_pb1
            @test sgn_pb0 == sgn_b0p
            @test sgn_pb1 == sgn_b1p
        end
    end
end