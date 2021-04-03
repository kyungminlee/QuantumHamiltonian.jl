using Test
using QuantumHamiltonian

using LinearAlgebra
using LatticeTools
using QuantumHamiltonian.Toolkit: pauli_matrix

@testset "symmetry_reduce" begin

  n_sites = 4;

  unitcell = makeunitcell(1.0; SiteType=String)
  addsite!(unitcell, "Spin", FractCoord([0], [0.0]))
  lattice = make_lattice(unitcell, n_sites)

  QN = Tuple{Int}
  up = State("Up", 1);
  dn = State("Dn",-1);
  spin_site = Site([up, dn]);

  hs = HilbertSpace([spin_site for i in 1:n_sites]);
  hss = HilbertSpaceSector(hs, 0)
  hsr = represent(hss);

  #translation_group = TranslationGroup([SitePermutation([2,3,4,1])])
  tsym = FiniteTranslationSymmetry(lattice)
  psym = project(PointSymmetryDatabase.find("-1"), [1 0 0;])

  tsymbed = embed(lattice, tsym)
  psymbed = embed(lattice, psym)
  ssymbed = tsymbed ⋊ psymbed

  for symred in [symmetry_reduce, symmetry_reduce_serial, symmetry_reduce_parallel]
    rhsr = symred(hsr, IrrepComponent(tsymbed, 1))
    @test rhsr.basis_list == UInt[0b0011, 0b0101]
    @test rhsr.parent === hsr

    rhsr = symred(hsr, IrrepComponent(tsymbed, 2))
    @test rhsr.basis_list == UInt[0b0011]
    @test rhsr.parent === hsr

    rhsr = symred(hsr, IrrepComponent(tsymbed, 3))
    @test rhsr.basis_list == UInt[0b0011, 0b0101]
    @test rhsr.parent === hsr

    rhsr = symred(hsr, IrrepComponent(tsymbed, 4))
    @test rhsr.basis_list == UInt[0b0011]
    @test rhsr.parent === hsr

    tol = Base.rtoldefault(Float64)
    count_translation = 0
    for tsic in get_irrep_components(tsymbed)
      rhsr = symred(hsr, tsic)
      count_translation += dimension(rhsr)
      for (i_p, b) in enumerate(hsr.basis_list)
        if b in rhsr.basis_list
          @test 1 <= rhsr.basis_mapping_index[i_p] <= dimension(rhsr)
          @test rhsr.basis_list[rhsr.basis_mapping_index[i_p]] == b
          @test isapprox(imag(rhsr.basis_mapping_amplitude[i_p]), 0; atol=tol)
        end
      end
      for (i_p, i_r) in enumerate(rhsr.basis_mapping_index)
        if i_r == -1
          @test ! (rhsr.parent.basis_list[i_p] in rhsr.basis_list)
        end
      end

      @test_throws DimensionMismatch symmetry_unreduce(rhsr, zeros(ComplexF64, dimension(rhsr)+1))
      @test_throws DimensionMismatch symmetry_reduce(rhsr, zeros(ComplexF64, dimension(hsr)+1))
      @test_throws DimensionMismatch symmetry_reduce_serial(rhsr, zeros(ComplexF64, dimension(hsr)+1))
      @test_throws DimensionMismatch symmetry_reduce_parallel(rhsr, zeros(ComplexF64, dimension(hsr)+1))

      sv = rand(ComplexF64, dimension(rhsr))
      lv = symmetry_unreduce(rhsr, sv)

      sv2 = symmetry_reduce(rhsr, lv)
      @test isapprox(sv, sv2; atol=tol)
      sv3 = symmetry_reduce_serial(rhsr, lv)
      @test isapprox(sv, sv3; atol=tol)
      sv4 = symmetry_reduce_parallel(rhsr, lv)
      @test isapprox(sv, sv4; atol=tol)

      sv5 = zero(sv4)
      @test sv5 === symmetry_reduce!(sv5, rhsr, lv)
      @test isapprox(sv, sv5; atol=tol)
      sv6 = zero(sv4)
      @test sv6 === symmetry_reduce_serial!(sv6, rhsr, lv)
      @test isapprox(sv, sv6; atol=tol)
      sv7 = zero(sv4)
      @test sv7 === symmetry_reduce_parallel!(sv7, rhsr, lv)
      @test isapprox(sv, sv7; atol=tol)

      @test_throws DimensionMismatch symmetry_reduce!(zero(sv4), rhsr, Float64[1.0])
      @test_throws DimensionMismatch symmetry_reduce_serial!(zero(sv4), rhsr, Float64[1.0])
      @test_throws DimensionMismatch symmetry_reduce_parallel!(zero(sv4), rhsr, Float64[1.0])

      @test_throws DimensionMismatch symmetry_reduce!(zeros(ComplexF64, dimension(rhsr)+1), rhsr, lv)
      @test_throws DimensionMismatch symmetry_reduce_serial!(zeros(ComplexF64, dimension(rhsr)+1), rhsr, lv)
      @test_throws DimensionMismatch symmetry_reduce_parallel!(zeros(ComplexF64, dimension(rhsr)+1), rhsr, lv)
    end

    #                               0011,   0101,   0110,   1001,   1010,   1100
    #                               1001    0101    1100    0011    1010    0110
    rhsr = symred(hsr, IrrepComponent(psymbed, 1))
    @test rhsr.basis_list == UInt[0b0011, 0b0101, 0b0110, 0b1010]
    rhsr = symred(hsr, IrrepComponent(psymbed, 2))
    @test rhsr.basis_list == UInt[0b0011, 0b0110]

    count_point = 0
    for psic in get_irrep_components(psymbed)
      rhsr = symred(hsr, psic)
      count_point += dimension(rhsr)
      for (i_p, b) in enumerate(hsr.basis_list)
        if b in rhsr.basis_list
          @test 1 <= rhsr.basis_mapping_index[i_p] <= dimension(rhsr)
          @test rhsr.basis_list[rhsr.basis_mapping_index[i_p]] == b
          @test isapprox(imag(rhsr.basis_mapping_amplitude[i_p]), 0; atol=tol)
        end
      end
      for (i_p, i_r) in enumerate(rhsr.basis_mapping_index)
        if i_r == -1
          @test ! (rhsr.parent.basis_list[i_p] in rhsr.basis_list)
        end
      end
    end # for psic

    # space symmetry
    let ssic = first(get_irrep_components(ssymbed))
      rhsr = symred(hsr, ssic)
      @test rhsr.basis_list == UInt[0b0011, 0b0101]
    end
    let ssic = SymmorphicIrrepComponent(IrrepComponent(tsymbed, 1, 1),
                                        IrrepComponent(psymbed, 1, 1))
      rhsr = symred(hsr, ssic)
      @test rhsr.basis_list == UInt[0b0011, 0b0101]
    end
    let ssic = SymmorphicIrrepComponent(IrrepComponent(tsymbed, 1, 1),
                                        IrrepComponent(psymbed, 2, 1))
      rhsr = symred(hsr, ssic)
      @test rhsr.basis_list == UInt[]
    end

    count_space = 0
    for ssic in get_irrep_components(ssymbed)
      rhsr = symred(hsr, ssic)
      count_space += dimension(rhsr)
      for (i_p, b) in enumerate(hsr.basis_list)
        if b in rhsr.basis_list
          @test 1 <= rhsr.basis_mapping_index[i_p] <= dimension(rhsr)
          @test rhsr.basis_list[rhsr.basis_mapping_index[i_p]] == b
          @test isapprox(imag(rhsr.basis_mapping_amplitude[i_p]), 0; atol=tol)
        end
      end
      for (i_p, i_r) in enumerate(rhsr.basis_mapping_index)
        if i_r == -1
          @test ! (rhsr.parent.basis_list[i_p] in rhsr.basis_list)
        end
      end
    end # for ssic

    @test count_translation == count_point == count_space
  end

  @testset "convention" begin
    # want  |ψ(k)⟩ = ∑ exp(+ikx) |ψ(x)⟩

    n_sites = 7;
    unitcell = makeunitcell(1.0; SiteType=String)
    addsite!(unitcell, "Spin", FractCoord([0], [0.0]))
    lattice = make_lattice(unitcell, n_sites)

    QN = Tuple{Int}
    up = State("Up", 1);
    dn = State("Dn",-1);
    spin_site = Site([up, dn]);

    hs = HilbertSpace([spin_site for i in 1:n_sites]);

    let
      hss = HilbertSpaceSector(hs, 5)
      hsr = represent(hss)
      p = SitePermutation([2,3,4,5,6,7,1])
      #translation_group = TranslationGroup(p)
      @test symmetry_apply(hs, p, 0b0000001) == (0b0000010, true)

      tsym = FiniteTranslationSymmetry(lattice)
      tsymbed = embed(lattice, tsym)
      tsic = IrrepComponent(tsymbed, 2, 1)
      rhsr = symmetry_reduce(hsr, tsic)
      @test hsr.basis_list == UInt[0b0000001, 0b0000010, 0b0000100, 0b0001000, 0b0010000, 0b0100000, 0b1000000]
      @test rhsr.basis_list == UInt[0b0000001]
      ψk = symmetry_unreduce(rhsr, [1.0])
      @test isapprox(ψk, [cis(2π * i/n_sites) / sqrt(n_sites) for i in 0:(n_sites-1)]; atol=Base.rtoldefault(Float64))
    end

    let
      hss = HilbertSpaceSector(hs, -5)
      hsr = represent(hss)
      p = SitePermutation([2,3,4,5,6,7,1])
      @test symmetry_apply(hs, p, 0b0000001) == (0b0000010, true)

      tsym = FiniteTranslationSymmetry(lattice)
      tsymbed = embed(lattice, tsym)
      tsic = IrrepComponent(tsymbed, 2, 1)
      rhsr = symmetry_reduce(hsr, tsic)

      # opposite order
      @test hsr.basis_list == UInt[0b0111111, 0b1011111, 0b1101111, 0b1110111, 0b1111011, 0b1111101, 0b1111110]
      @test rhsr.basis_list == UInt[0b0111111]
      ψk = symmetry_unreduce(rhsr, [1.0])
      @test isapprox(ψk, [cis(-2π * i/n_sites) / sqrt(n_sites) for i in 0:(n_sites-1)]; atol=Base.rtoldefault(Float64))
    end

  end

  # σ = Dict( (isite, j) => pauli_matrix(hs, isite, j) for isite in 1:n_sites, j in [:x, :y, :z, :+, :-]);
  # j1 = simplify( sum(σ[(i, j)] * σ[( mod(i, n_sites) + 1 , j)] for i in 1:n_sites, j in [:x, :y, :z]) );
end
