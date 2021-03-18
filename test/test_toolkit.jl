using Test
using QuantumHamiltonian

@testset "Toolkit" begin
  @testset "SpinHalf" begin
    @testset "small" begin
      n_sites = 4
      up = State("Up", 1)
      dn = State("Dn",-1)
      spin_site = Site([up, dn])

      hs1 = HilbertSpace([spin_site for i in 1:n_sites])
      (hs2, pauli) = QuantumHamiltonian.Toolkit.spin_half_system(n_sites)

      @test hs1 == hs2
      for i_site in 1:n_sites
        @test pauli(i_site, :x) == pure_operator(hs1, i_site, 1, 2, 1, UInt) + pure_operator(hs1, i_site, 2, 1, 1, UInt)
        @test pauli(i_site, :y) == pure_operator(hs1, i_site, 1, 2, -im, UInt) + pure_operator(hs1, i_site, 2, 1, +im, UInt)
        @test pauli(i_site, :z) == pure_operator(hs1, i_site, 1, 1, 1, UInt) + pure_operator(hs1, i_site, 2, 2, -1, UInt)
        @test pauli(i_site, :+) == pure_operator(hs1, i_site, 1, 2, 1, UInt)
        @test pauli(i_site, :-) == pure_operator(hs1, i_site, 2, 1, 1, UInt)
      end
      @test_throws ArgumentError pauli(1, :unknown)
    end

    @testset "binary type" begin
      (hs1, pauli1) = QuantumHamiltonian.Toolkit.spin_half_system(63)
      @test_throws ArgumentError QuantumHamiltonian.Toolkit.spin_half_system(128)
      (hs2, pauli2) = QuantumHamiltonian.Toolkit.spin_half_system(128, UInt128)
      @test bintype(pauli1(1, :x)) == UInt64
      @test bintype(pauli2(1, :x)) == UInt128
    end
  end

  @testset "SpinSystem" begin
    using QuantumHamiltonian.Toolkit: spin_system
    @testset "Constructor" begin
      spin_system(33, 1//2)
      @test_throws ArgumentError spin_system(33, 1)
      spin_system(33, 1, UInt128)
    end
    r2, r3, r6 = sqrt(2), sqrt(3), sqrt(6)

    spin_matrices = Dict(
      1//2 => Dict(
          :x => [0 1; 1 0] ./2,
          :y => [0 1; -1 0] ./(2im),
          :z => [1 0; 0 -1] ./2,
          :+ => [0 1; 0 0],
          :- => [0 0; 1 0],
        ),
      1//1 => Dict(
          :x => [0 1 0; 1 0 1; 0 1 0] ./ r2,
          :y => [0 1 0; -1 0 1; 0 -1 0] ./ (r2*1im),
          :z => [1 0 0; 0 0 0; 0 0 -1],
          :+ => [0 1 0; 0 0 1; 0 0 0] .* r2,
          :- => [0 0 0; 1 0 0; 0 1 0] .* r2,
        ),
      3//2 => Dict(
          :x => [0 r3 0 0; r3 0 2 0; 0 2 0 r3; 0 0 r3 0] ./ 2,
          :y => [0 r3 0 0;-r3 0 2 0; 0 -2 0 r3; 0 0 -r3 0] ./ (2im),
          :z => [3/2 0 0 0; 0 1/2 0 0; 0 0 -1/2 0; 0 0 0 -3/2],
          :+ => [0 r3 0 0; 0 0 2 0; 0 0 0 r3; 0 0 0 0],
          :- => [0 0 0 0; r3 0 0 0; 0 2 0 0; 0 0 r3 0],
        ),
      2//1 => Dict(
          :x => [0 2 0 0 0; 2 0 r6 0 0; 0 r6 0 r6 0; 0 0 r6 0 2; 0 0 0 2 0] ./ 2,
          :y => [0 2 0 0 0; -2 0 r6 0 0; 0 -r6 0 r6 0; 0 0 -r6 0 2; 0 0 0 -2 0] ./ (2im),
          :z => [2 0 0 0 0; 0 1 0 0 0; 0 0 0 0 0; 0 0 0 -1 0; 0 0 0 0 -2],
          :+ => [0 2 0 0 0; 0 0 r6 0 0; 0 0 0 r6 0; 0 0 0 0 2; 0 0 0 0 0],
          :- => [0 0 0 0 0; 2 0 0 0 0; 0 r6 0 0 0; 0 0 r6 0 0; 0 0 0 2 0],
        ),
    )

    @testset "SingleSite" begin
      @test_throws ArgumentError QuantumHamiltonian.Toolkit.spin_system(1, -3)
      @test_throws ArgumentError QuantumHamiltonian.Toolkit.spin_system(1, 1//3)
      let
        (hs, spin) = QuantumHamiltonian.Toolkit.spin_system(1, 1//2)
        @test_throws ArgumentError spin(1, :q)
      end
      for S in Any[1//2, 1, 1//1, 3//2, 2//1, 2]
        (hs, spin) = QuantumHamiltonian.Toolkit.spin_system(1, S)
        hsr = represent(hs)
        for μ in [:x, :y, :z, :+, :-]
          @test isapprox(spin_matrices[Rational(S)][μ], Matrix(represent(hsr, spin(1, μ))))
        end
      end # for S
    end # testset SingleSite

    @testset "MultipleSiteCommutation" begin
      levi_civita = zeros(Int, (3,3,3))
      levi_civita[1,2,3] = 1
      levi_civita[2,1,3] = -1
      levi_civita[2,3,1] = 1
      levi_civita[1,3,2] = -1
      levi_civita[3,1,2] = 1
      levi_civita[3,2,1] = -1
      component_index = Dict(:x => 1, :y => 2, :z => 3)
      for S in [1//2, 1, 1//1, 3//2, 2//1, 2]
        n = 4
        (hs, spin) = QuantumHamiltonian.Toolkit.spin_system(n, S)
        @test all(let
                    op1 = simplify(spin(i, μ) * spin(j, ν) - spin(j, ν) * spin(i, μ))
                    if i == j
                      op2 = im * sum( levi_civita[iμ, iν, iρ] * spin(i, ρ) for (iρ, ρ) in enumerate([:x, :y, :z]) )
                    else
                      op2 = NullOperator()
                    end
                    simplify(op1 - op2) == NullOperator()
                  end
                    for i in 1:n for j in i:n
                    for (iμ, μ) in enumerate([:x, :y, :z]), (iν, ν) in enumerate([:x, :y, :z]))
      end # for S
    end # testset MultipleSite
  end # testset SpinSystem

  @testset "product_state" begin
    (hs, pauli) = QuantumHamiltonian.Toolkit.spin_half_system(4)
    @test_throws ArgumentError QuantumHamiltonian.Toolkit.product_state(hs, [[1.0, 0.0], [0.0, 0.0]]) # too few
    @test_throws ArgumentError QuantumHamiltonian.Toolkit.product_state(hs, [[1.0, 0.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]) # too many
    @test_throws ArgumentError QuantumHamiltonian.Toolkit.product_state(hs, [[1.0, 0.0], [0.0, 1.0], [0.0, 1.0, 0.0], [0.0, 1.0]])

    local_states = [[1.0, 0.0], [1.0 + 10.0im, 2.0], [1.0, 0.0], [0.0, 1.0]]
    @testset "HilbertSpace" begin
      psi1 = QuantumHamiltonian.Toolkit.product_state(hs, local_states)
      psi2 = SparseState{ComplexF64, UInt}(0b1000 => 1.0 + 10.0im, 0b1010 => 2.0)
      @test isapprox(psi1, psi2; atol=1E-8)
    end

    @testset "HilbertSpaceRepresentation" begin
      hsr = represent(hs)
      psi1 = QuantumHamiltonian.Toolkit.product_state(hsr, local_states)
      psi2 = zeros(ComplexF64, 16)
      psi2[9] = 1.0 + 10.0im
      psi2[11] = 2.0
      @test isapprox(psi1, psi2; atol=1E-8)

      @test_throws ArgumentError QuantumHamiltonian.Toolkit.product_state(hsr, [[1.0, 0.0], [0.0, 0.0]]) # too few
      @test_throws ArgumentError QuantumHamiltonian.Toolkit.product_state(hsr, [[1.0, 0.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]) # too many
      @test_throws ArgumentError QuantumHamiltonian.Toolkit.product_state(hsr, [[1.0, 0.0], [0.0, 1.0], [0.0, 1.0, 0.0], [0.0, 1.0]])
    end
    @testset "sector-rep" begin
      hssr = represent(HilbertSpaceSector(hs, 0))
      # 0011, 0101, 0110, 1001, 1010, 1100
      psi1 = QuantumHamiltonian.Toolkit.product_state(hssr, local_states)
      psi2 = zeros(ComplexF64, 6)
      psi2[5] = 2.0
      @test isapprox(psi1, psi2; atol=1E-8)
    end
  end
end
