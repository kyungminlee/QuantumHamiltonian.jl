using Test
using Suppressor
using QuantumHamiltonian

using QuantumHamiltonian.Toolkit: pauli_matrix

@testset "prettyprintln" begin
  @testset "spinhalf" begin

    n_sites = 4
    unitcell = make_unitcell(1.0; SiteType=String)
    addsite!(unitcell, "Spin", FractCoord([0], [0.0]))
    lattice = make_lattice(unitcell, n_sites)
    tsym = TranslationSymmetry(lattice)
    tsymbed = embed(lattice, tsym)

    up = State("Up",  1)
    dn = State("Dn", -1)
    spin_site = Site([up, dn])
    hs = HilbertSpace([spin_site, spin_site, spin_site, spin_site])

    σ(i::Integer, j::Symbol) = pauli_matrix(hs, i, j)

    buf = IOBuffer()

    val = NullOperator()
    result1 = @capture_out prettyprintln(val)
    prettyprintln(buf, val)
    result2 = String(take!(buf))
    @test result1 == result2
    @test result1 == "NullOperator\n"

    val = σ(1, :+)
    result1 = @capture_out prettyprintln(val)
    prettyprintln(buf, val)
    result2 = String(take!(buf))
    @test result1 == result2
    @test result1 == join(["PureOperator",
                           "| M: 0000000000000000000000000000000000000000000000000000000000000001",
                           "| R: 0000000000000000000000000000000000000000000000000000000000000000",
                           "| C: 0000000000000000000000000000000000000000000000000000000000000001",
                           "| A: 1", ""], "\n")
    val = σ(2, :x)
    result1 = @capture_out prettyprintln(val)
    prettyprintln(buf, val)
    result2 = String(take!(buf))
    @test result1 == result2
    @test result1 == join(["SumOperator",
                           "| PureOperator",
                           "| | M: 0000000000000000000000000000000000000000000000000000000000000010",
                           "| | R: 0000000000000000000000000000000000000000000000000000000000000000",
                           "| | C: 0000000000000000000000000000000000000000000000000000000000000010",
                           "| | A: 1",
                           "| PureOperator",
                           "| | M: 0000000000000000000000000000000000000000000000000000000000000010",
                           "| | R: 0000000000000000000000000000000000000000000000000000000000000010",
                           "| | C: 0000000000000000000000000000000000000000000000000000000000000000",
                           "| | A: 1",
                           ""], "\n")

    hsr = represent(hs)
    result = @capture_out prettyprintln(hsr)
    @test result == join(["HilbertSpaceRepresentation",
                          "| 0000",
                          "| 0001",
                          "| 0010",
                          "| 0011",
                          "| 0100",
                          "| 0101",
                          "| 0110",
                          "| 0111",
                          "| 1000",
                          "| 1001",
                          "| 1010",
                          "| 1011",
                          "| 1100",
                          "| 1101",
                          "| 1110",
                          "| 1111",
                          ""], "\n")

    rhsr = symmetry_reduce(hsr, IrrepComponent(tsymbed, 1))
    result = @capture_out prettyprintln(rhsr)
    @test result == join(["ReducedHilbertSpaceRepresentation",
                          "| basis_list",
                          "| | 0000",
                          "| | 0001",
                          "| | 0011",
                          "| | 0101",
                          "| | 0111",
                          "| | 1111",
                          "| basis_mapping",
                          "| | 1: 1, 1.0 - 0.0im",
                          "| | 2: 2, 0.5 + 0.0im",
                          "| | 3: 2, 0.5 - 0.0im",
                          "| | 4: 3, 0.5 + 0.0im",
                          "| | 5: 2, 0.5 - 0.0im",
                          "| | 6: 4, 0.7071067811865475 - 0.0im",
                          "| | 7: 3, 0.5 - 0.0im",
                          "| | 8: 5, 0.5 + 0.0im",
                          "| | 9: 2, 0.5 - 0.0im",
                          "| | 10: 3, 0.5 - 0.0im",
                          "| | 11: 4, 0.7071067811865475 - 0.0im",
                          "| | 12: 5, 0.5 - 0.0im",
                          "| | 13: 3, 0.5 - 0.0im",
                          "| | 14: 5, 0.5 - 0.0im",
                          "| | 15: 5, 0.5 - 0.0im",
                          "| | 16: 6, 1.0 - 0.0im",
                          ""], "\n")

    # val = SparseState{Float64, UInt}(hs, UInt(0b0010)=>0.2, UInt(0b100) => 0.3 )
    # result1 = @capture_out prettyprintln(val)
    # prettyprintln(buf, val)
    # result2 = String(take!(buf))
    # @test result1 == result2
    # @test result1 == join(["SparseState",
    #                        "| 0010 : 0.2",
    #                        "| 0100 : 0.3", ""], "\n")
  end
end
