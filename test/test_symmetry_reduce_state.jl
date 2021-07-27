using QuantumHamiltonian
using LatticeTools
using Test

@testset "symmetry_reduce (basis)" begin
    n = 6
    uc = makeunitcell(1.0)
    addsite!(uc, "A", FractCoord([0,], [0.0]))
    lattice = makelattice(uc, n)
    
    tsym = FiniteTranslationSymmetry(lattice)
    tsymbed = embed(lattice, tsym)
    tsym_irrep = [[t[1,1] for t in tr] for tr in irreps(tsymbed)]
    
    hs, pauli = QuantumHamiltonian.Toolkit.spin_half_system(n)
    hsr = represent(hs)
    phsr = let
        site = get_site(hs, 1)
        hsr3 = represent(HilbertSpace([site, site, site]))
        ProductHilbertSpaceRepresentation(hsr3, hsr3)
    end
    phs = basespace(phsr)
    dhsr = decompose(hsr)

    for tir in tsym_irrep
        rhsr = symmetry_reduce(hsr, tsymbed.elements, tir)
        rphsr = symmetry_reduce(phsr, tsymbed.elements, tir)
        rdhsr = symmetry_reduce(dhsr, tsymbed.elements, tir)
        for b in get_basis_iterator(hsr)
            (i1, a1) = get_basis_index_amplitude(rhsr, b)
            (i2, a2) = get_basis_index_amplitude(rphsr, b)
            (i3, a3) = get_basis_index_amplitude(rdhsr, b) # since symmetry operation conserves quantum number, order is conserved within each sector

            (b1p, a1p) = symmetry_reduce(hs, tsymbed.elements, tir, b)
            (i1p, x) = get_basis_index_amplitude(rhsr, b1p)
            @test imag(x) ≈ 0

            (b2p, a2p) = symmetry_reduce(phs, tsymbed.elements, tir, b)
            (i2p, x) = get_basis_index_amplitude(rphsr, b2p)
            @test imag(x) ≈ 0

            (i3p, x) = get_basis_index_amplitude(rdhsr, b1p)
            @test imag(x) ≈ 0

            @test i1 == i2 == i1p == i2p
            @test i3 == i3p
            @test a1 ≈ a2 ≈ a3 ≈ a1p ≈ a2p

            if i1 > 0
                b1 = get_basis_state(rhsr, i1)
                @test b1 == b1p == b2p
            else
                @test a1 ≈ a2 ≈ a3 ≈ 0
            end
            
        end
    end
end
