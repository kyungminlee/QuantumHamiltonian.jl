using Test
using QuantumHamiltonian

@testset "state_representation" begin
    n_sites = 4
    up = State("Up", 1)
    dn = State("Dn",-1)
    spin_site = Site([up, dn])
    hs = HilbertSpace([spin_site for i in 1:n_sites])
    hsr = represent(hs)

    ψ = SparseState{ComplexF64, UInt}(0x0 => 1.0, 0x2 => 2.0+3.0im, 0x1000=>4.0+5.0im)
    (ψ_rep, err, err_sq) = represent(hsr, ψ)
    let
        ψ_rep2 = zeros(ComplexF64, 2^n_sites)
        ψ_rep2[1] = 1.0
        ψ_rep2[3] = 2.0 + 3.0im
        @test ψ_rep == ψ_rep2
    end
    @test err == 4.0 + 5.0im && err_sq == abs2(4.0+5.0im)

    ψ2 = SparseState{ComplexF64, UInt}(0x0 => 1.0, 0x2 => 2.0+3.0im)
    ψ3 = SparseState(hsr, ψ_rep)
    @test ψ2 == ψ3

    @test_throws DimensionMismatch SparseState(hsr, zeros(Float64, 999))
end
