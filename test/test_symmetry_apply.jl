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
    hss = HilbertSpaceSector(hs, 0)
    transop = SitePermutation([2,3,4,1])
    invop = SitePermutation([1,4,3,2])

    @test symmetry_apply(hs, transop, 0b0001) == (0b0010, true)
    @test symmetry_apply(hss, transop, 0b0101) == (0b1010, true)

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

    @test isinvariant(hs, transop, nop)
    @test isinvariant(HilbertSpaceSector(hs, 0), transop, nop)
    @test !isinvariant(hs, transop, pop1)
    @test !isinvariant(HilbertSpaceSector(hs, 0), transop, pop1)
    @test !isinvariant(hs, transop, sop)
    @test !isinvariant(HilbertSpaceSector(hs, 0), transop, sop)
    @test isinvariant(hs, transop, j1)
    @test isinvariant(HilbertSpaceSector(hs, 0), transop, j1)

    @test isinvariant(hs, invop, nop)
    @test isinvariant(HilbertSpaceSector(hs, 0), invop, nop)
    @test !isinvariant(hs, invop, pop1)
    @test !isinvariant(HilbertSpaceSector(hs, 0), invop, pop1)
    @test !isinvariant(hs, invop, sop)
    @test !isinvariant(HilbertSpaceSector(hs, 0), invop, sop)
    @test isinvariant(hs, invop, j1)
    @test isinvariant(HilbertSpaceSector(hs, 0), invop, j1)

    tsymbed = translation_symmetry_embedding(lattice)

    @test isinvariant(hs, tsymbed, nop)
    @test isinvariant(HilbertSpaceSector(hs, 0), tsymbed, nop)
    @test !isinvariant(hs, tsymbed, pop1)
    @test !isinvariant(HilbertSpaceSector(hs, 0), tsymbed, pop1)
    @test !isinvariant(hs, tsymbed, sop)
    @test !isinvariant(HilbertSpaceSector(hs, 0), tsymbed, sop)
    @test isinvariant(hs, tsymbed, j1)
    @test isinvariant(HilbertSpaceSector(hs, 0), tsymbed, j1)

    psymbed = embed(lattice, project(PointSymmetryDatabase.get(2), [1 0 0;]))

    @test isinvariant(hs, psymbed, nop)
    @test isinvariant(HilbertSpaceSector(hs, 0), psymbed, nop)
    @test !isinvariant(hs, psymbed, pop1)
    @test !isinvariant(HilbertSpaceSector(hs, 0), psymbed, pop1)
    @test !isinvariant(hs, psymbed, sop)
    @test !isinvariant(HilbertSpaceSector(hs, 0), psymbed, sop)
    @test isinvariant(hs, psymbed, j1)
    @test isinvariant(HilbertSpaceSector(hs, 0), psymbed, j1)

    ssymbed = tsymbed ⋊ psymbed
    @test isinvariant(hs, ssymbed, nop)
    @test isinvariant(HilbertSpaceSector(hs, 0), ssymbed, nop)
    @test !isinvariant(hs, ssymbed, pop1)
    @test !isinvariant(HilbertSpaceSector(hs, 0), ssymbed, pop1)
    @test !isinvariant(hs, ssymbed, sop)
    @test !isinvariant(HilbertSpaceSector(hs, 0), ssymbed, sop)
    @test isinvariant(hs, ssymbed, j1)
    @test isinvariant(HilbertSpaceSector(hs, 0), ssymbed, j1)

end
