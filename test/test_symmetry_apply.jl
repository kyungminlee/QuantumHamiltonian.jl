using Test
using QuantumHamiltonian

using LatticeTools
using QuantumHamiltonian.Toolkit: pauli_matrix

@testset "symmetry_apply" begin

    unitcell = makeunitcell(1.0; SiteType=String)
    addsite!(unitcell, "Spin", FractCoord([0], [0.0]))
    lattice = make_lattice(unitcell, 4)

    # Test State and Site
    up = State("Up", 1)
    dn = State("Dn",-1)
    spin_site = Site([up, dn])

    hs = HilbertSpace(repeat([spin_site], 4))
    hsr = represent(hs)

    transop = SitePermutation([2,3,4,1])
    invop = SitePermutation([1,4,3,2])

    @test symmetry_apply(hs, transop, 0b0001) == (0b0010, true)

    nop = NullOperator()
    pop1 = PureOperator{Float64, UInt}(0b1101, 0b0101, 0b1100, 2.0)
    pop2 = PureOperator{Float64, UInt}(0b0010, 0b0000, 0b0010, 3.0)
    sop = pop1 + pop2

    @test symmetry_apply(hs, transop, nop) == nop
    @test symmetry_apply(hs, transop, pop1) == PureOperator{Float64, UInt}(0b1011, 0b1010, 0b1001, 2.0)
    @test symmetry_apply(hs, transop, sop) == SumOperator{Float64, UInt}([
        symmetry_apply(hs, transop, pop1), symmetry_apply(hs, transop, pop2)
    ])

    @test symmetry_apply(hs, invop, nop) == nop
    @test symmetry_apply(hs, invop, pop1) == PureOperator{Float64, UInt}(0b0111, 0b0101, 0b0110, 2.0)
    @test symmetry_apply(hs, invop, sop) == SumOperator{Float64, UInt}([
        symmetry_apply(hs, invop, pop1), symmetry_apply(hs, invop, pop2)
    ])

    σ = Dict( (isite, j) => pauli_matrix(hs, isite, j) for isite in 1:4, j in [:x, :y, :z, :+, :-])
    j1 = simplify(sum(σ[i, j] * σ[mod(i, 4) + 1 , j] for i in 1:4, j in [:x, :y, :z]))

    tsymbed = translation_symmetry_embedding(lattice)
    psymbed = embed(lattice, project(PointSymmetryDatabase.get(2), [1 0 0;]))
    ssymbed = SymmorphicSymmetryEmbedding(tsymbed, psymbed)

    for hsobj in [hs, hsr]
        @test isinvariant(hsobj, transop, nop)
        @test !isinvariant(hsobj, transop, pop1)
        @test !isinvariant(hsobj, transop, sop)
        @test isinvariant(hsobj, transop, j1)

        @test isinvariant(hsobj, invop, nop)
        @test !isinvariant(hsobj, invop, pop1)
        @test !isinvariant(hsobj, invop, sop)
        @test isinvariant(hsobj, invop, j1)

        @test isinvariant(hsobj, tsymbed, nop)
        @test !isinvariant(hsobj, tsymbed, pop1)
        @test !isinvariant(hsobj, tsymbed, sop)
        @test isinvariant(hsobj, tsymbed, j1)

        @test isinvariant(hsobj, psymbed, nop)
        @test !isinvariant(hsobj, psymbed, pop1)
        @test !isinvariant(hsobj, psymbed, sop)
        @test isinvariant(hsobj, psymbed, j1)

        @test isinvariant(hsobj, ssymbed, nop)
        @test !isinvariant(hsobj, ssymbed, pop1)
        @test !isinvariant(hsobj, ssymbed, sop)
        @test isinvariant(hsobj, ssymbed, j1)
    end
end
