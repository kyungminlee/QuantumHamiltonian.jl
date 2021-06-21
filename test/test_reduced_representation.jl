using Test
using QuantumHamiltonian

using Suppressor
using LinearAlgebra
using SparseArrays

using LatticeTools

using QuantumHamiltonian.Toolkit: pauli_matrix

@testset "RedRep4" begin
    tol = Base.rtoldefault(Float64)
    
    n = 4
    unitcell = makeunitcell(1.0; SiteType=String)
    addsite!(unitcell, "Spin", FractCoord([0], [0.0]))
    lattice = make_lattice(unitcell, n)
    
    # Test State and Site
    up = State("Up", 1)
    dn = State("Dn",-1)
    spin_site = Site([up, dn])
    
    hs = HilbertSpace(repeat([spin_site], n))
    σ = Dict( (isite, j) => pauli_matrix(hs, isite, j) for isite in 1:n, j in [:x, :y, :z, :+, :-])
    j1 = sum(σ[i, j] * σ[mod(i, n) + 1 , j] for i in 1:n, j in [:x, :y, :z])
    j2 = sum(σ[i, j] * σ[mod(i+1, n) + 1 , j] for i in 1:n, j in [:x, :y, :z])
    
    hsr = represent(HilbertSpaceSector(hs, 0))
    tsym = FiniteTranslationSymmetry(lattice)
    tsymbed = embed(lattice, tsym)
    psym = project(PointSymmetryDatabase.get3d(2), [1 0 0;])
    psymbed = embed(lattice, psym)
    # translation_group = TranslationGroup([Permutation([2,3,4,1])])
    # @test isinvariant(hs, translation_group, j1)
    # @test isinvariant(HilbertSpaceSector(hs, 0), translation_group, j1)
    
    @testset "RHSR-old implementation" begin
        #rhsr = symmetry_reduce(hsr, translation_group, [0//1])
        for symred in [symmetry_reduce, symmetry_reduce_serial, symmetry_reduce_parallel]
            rhsr = symred(hsr, IrrepComponent(tsymbed, 2, 1))
            @test scalartype(rhsr) === ComplexF64
            @test scalartype(typeof(rhsr)) === ComplexF64
            @test valtype(rhsr) === ComplexF64
            @test valtype(typeof(rhsr)) === ComplexF64
            @test bintype(rhsr) === UInt
            @test bintype(typeof(rhsr)) === UInt
            @test dimension(rhsr) <= 2^n
            @test bitwidth(rhsr) == n
            show(devnull, MIME("text/plain"), rhsr)  # make sure it doesn't crash
            
            rhsr = symred(hsr, IrrepComponent(psymbed, 2, 1))
            @test scalartype(rhsr) === ComplexF64
            @test scalartype(typeof(rhsr)) === ComplexF64
            @test valtype(rhsr) === ComplexF64
            @test valtype(typeof(rhsr)) === ComplexF64
            @test bintype(rhsr) === UInt
            @test bintype(typeof(rhsr)) === UInt
            @test dimension(rhsr) <= 2^n
            @test bitwidth(rhsr) == n
            show(devnull, MIME("text/plain"), rhsr)  # make sure it doesn't crash
        end
    end
    
    @testset "RHSR" begin
        # symops_and_amplitudes = collect(get_irrep_iterator(IrrepComponent(tsymbed, 2, 1)))
        symops = elements(tsymbed)
        amplitudes = [m[1,1] for m in irrep(tsymbed, 2)]
        let rhsr = symmetry_reduce(hsr, symops, amplitudes)
            @test get_basis_state(rhsr, 1) == 0b0011
            @test_throws BoundsError get_basis_state(rhsr, 2)
            let out = get_basis_index_amplitude(rhsr, 0b0011)
                @test out.index == 1
                @test isapprox(out.amplitude, 0.5)
            end
            let out = get_basis_index_amplitude(rhsr, 0b1100)
                @test out.index == 1
                @test isapprox(out.amplitude, -0.5)
            end
            let out = get_basis_index_amplitude(rhsr, 0b0101)
                @test out.index < 0
                @test iszero(out.amplitude)
            end
            let out = get_basis_index_amplitude(rhsr, 0b0000)
                @test out.index < 0
                @test iszero(out.amplitude)
            end
        end
        for symred in [symmetry_reduce, symmetry_reduce_serial, symmetry_reduce_parallel]
            rhsr = symred(hsr, symops, amplitudes)
            @test scalartype(rhsr) === ComplexF64
            @test scalartype(typeof(rhsr)) === ComplexF64
            @test valtype(rhsr) === ComplexF64
            @test valtype(typeof(rhsr)) === ComplexF64
            @test bintype(rhsr) === UInt
            @test bintype(typeof(rhsr)) === UInt
            @test dimension(rhsr) <= 2^n
            @test bitwidth(rhsr) == n
            show(devnull, MIME("text/plain"), rhsr)  # make sure it doesn't crash  
        end
        # symops_and_amplitudes = collect(get_irrep_iterator(IrrepComponent(psymbed, 2, 1)))
        symops = elements(psymbed)
        amplitudes = [m[1,1] for m in irrep(psymbed, 2)]
        for symred in [symmetry_reduce, symmetry_reduce_serial, symmetry_reduce_parallel]
            rhsr = symred(hsr, symops, amplitudes)
            @test scalartype(rhsr) === ComplexF64
            @test scalartype(typeof(rhsr)) === ComplexF64
            @test valtype(rhsr) === ComplexF64
            @test valtype(typeof(rhsr)) === ComplexF64
            @test bintype(rhsr) === UInt
            @test bintype(typeof(rhsr)) === UInt
            @test dimension(rhsr) <= 2^n
            @test bitwidth(rhsr) == n
            show(devnull, MIME("text/plain"), rhsr)  # make sure it doesn't crash

            @test_throws ArgumentError symred(hsr, symops, amplitudes[1:end-1])
            @test_throws ArgumentError symred(hsr, empty(symops), empty(amplitudes))
            @test_throws ArgumentError symred(hsr, symops, [amplitudes[1:end-1]..., 0.2])
            symred(hsr, symops, [1.0, amplitudes[2:end]...])
            @test_throws ArgumentError symred(hsr, symops, [0.0, amplitudes[2:end]...])
        end
    end
    
    j1_rep = represent(hsr, j1)
    j1_mat = Matrix(j1_rep)
    @test hsr.basis == [0b0011, 0b0101, 0b0110, 0b1001, 0b1010, 0b1100]
    
    @testset "ROR" begin
        #rhsr = symmetry_reduce(hsr, translation_group, [0//1])
        rhsr = symmetry_reduce(hsr, IrrepComponent(tsymbed, 1, 1))
        @test dimension(rhsr) == 2
        @test collect(rhsr.basis) == UInt[0b0011, 0b0101]
        dim = dimension(rhsr)
        
        j1_redrep = represent(rhsr, j1)
        j2_redrep = represent(rhsr, j2)
        show(devnull, MIME("text/plain"), j1_redrep) # make sure it doesn't crash
        
        @testset "typetraits" begin
            @test scalartype(j1) === Complex{Int}
            @test valtype(j1) === Complex{Int}
            
            @test scalartype(j1_redrep) === ComplexF64
            @test scalartype(typeof(j1_redrep)) === ComplexF64
            @test valtype(j1_redrep) === ComplexF64
            @test valtype(typeof(j1_redrep)) === ComplexF64
            
            @test bintype(j1_redrep) === UInt
            @test bintype(typeof(j1_redrep)) === UInt
            
            @test spacetype(j1_redrep) === typeof(rhsr)
            @test operatortype(j1_redrep) === typeof(j1)
            @test spacetype(typeof(j1_redrep)) === typeof(rhsr)
            @test operatortype(typeof(j1_redrep)) === typeof(j1)
            @test get_space(j1_redrep) === rhsr
        end
        
        @testset "property" begin
            @test bitwidth(j1_redrep) == bitwidth(rhsr)
            @test bitwidth(j2_redrep) == bitwidth(rhsr)
            @test dimension(j1_redrep) == dimension(rhsr)
            @test dimension(j2_redrep) == dimension(rhsr)
        end
        
        @testset "sym" begin
            @test ishermitian(represent(hsr, σ[2, :x]))
            @test !ishermitian(represent(hsr, σ[2, :+]))
        end
        
        @testset "unary operator" begin
            @test +j1_redrep == j1_redrep
            @test -j1_redrep == represent(rhsr, -j1)
        end
        
        @testset "binary operator" begin
            @test simplify(j1_redrep + j2_redrep) == represent(rhsr, simplify(j1 + j2))
            @test simplify(j1_redrep - j2_redrep) == represent(rhsr, simplify(j1 - j2))
            @test simplify(j1_redrep * j2_redrep) == represent(rhsr, simplify(j1 * j2))
            
            @test simplify(j1_redrep * 0) == represent(rhsr, NullOperator())
            @test simplify(0 * j1_redrep) == represent(rhsr, NullOperator())
            @test simplify(j1_redrep * 2) == represent(rhsr, simplify(j1 * 2))
            @test simplify(2 * j1_redrep) == represent(rhsr, simplify(2 * j1))
            @test simplify(j1_redrep / 2) == represent(rhsr, simplify(j1 / 2))
            @test simplify(2 \ j1_redrep) == represent(rhsr, simplify(2 \ j1))
        end
        
        let
            psis = [normalize([1.0, 0.0, 1.0, 1.0, 0.0, 1.0]), normalize([0.0, 1.0, 0.0, 0.0, 1.0, 0.0])]
            H = zeros(ComplexF64, (2,2))
            for i in 1:2, j in 1:2
                H[i,j] = psis[i] ⋅ (j1_mat * psis[j])
            end
            @test isapprox(Matrix(j1_redrep), H; atol=tol)
            
            @testset "get_row_iterator" begin
                rowvec = Vector{ComplexF64}(undef, dimension(rhsr))
                for irow_r in 1:dimension(rhsr)
                    fill!(rowvec, zero(ComplexF64))
                    err = zero(ComplexF64)
                    for (icol_r, ampl) in get_row_iterator(j1_redrep, irow_r)
                        if 1 <= icol_r <= dimension(rhsr)
                            rowvec[icol_r] += ampl
                        else
                            err += ampl
                        end
                    end
                    @test isapprox(err, 0; atol=tol) # is this necessarily true?
                    @test isapprox(rowvec, H[irow_r, :]; atol=tol)
                end
            end
            
            @testset "get_column_iterator" begin
                colvec = zeros(ComplexF64, dimension(rhsr))
                for icol_r in 1:dimension(rhsr)
                    colvec[:] .= zero(ComplexF64)
                    
                    colvec[:] .= zero(ComplexF64)
                    err = zero(ComplexF64)
                    for (irow_r, ampl) in get_column_iterator(j1_redrep, icol_r)
                        if 1 <= irow_r <= dimension(rhsr)
                            colvec[irow_r] += ampl
                        else
                            err += ampl
                        end
                    end
                    @test isapprox(err, 0; atol=tol) # is this necessarily true?
                    @test isapprox(colvec, H[:, icol_r]; atol=tol)
                end
            end # testset get_column_iterator
            
            @testset "get_element" begin
                for irow_r in 1:dimension(rhsr)
                    for icol_r in 1:dimension(rhsr)
                        @test isapprox(get_element(j1_redrep, irow_r, icol_r), H[irow_r, icol_r]; atol=tol)
                        @test isapprox(j1_redrep[irow_r, icol_r], H[irow_r, icol_r]; atol=tol)
                    end
                end
            end
            
            @testset "get exceptions" begin
                dim = dimension(rhsr)
                opr = j1_redrep
                @test_throws BoundsError get_row_iterator(opr, 0)
                @test_throws BoundsError get_row_iterator(opr, dim+1)
                @test_throws BoundsError get_column_iterator(opr, 0)
                @test_throws BoundsError get_column_iterator(opr, dim+1)
                @test_throws BoundsError get_row(opr, 0)
                @test_throws BoundsError get_row(opr, dim+1)
                @test_throws BoundsError get_column(opr, 0)
                @test_throws BoundsError get_column(opr, dim+1)
                @test_throws BoundsError get_element(opr, 0, 1)
                @test_throws BoundsError get_element(opr, dim+1, 1)
                @test_throws BoundsError get_element(opr, 1, 0)
                @test_throws BoundsError get_element(opr, 1, dim+1)
                @test_throws BoundsError opr[0, 1]
                @test_throws BoundsError opr[dim+1, 1]
                @test_throws BoundsError opr[1, 0]
                @test_throws BoundsError opr[1, dim+1]
            end
        end
    end # testset ROR
    
    
    @testset "ROR-ALL" begin
        j1_mat = Matrix(represent(hsr, j1))
        @test isapprox(j1_mat, adjoint(j1_mat); atol=tol)
        eigenvalues1 = eigvals(Hermitian(j1_mat))
        eigenvalues2 = Float64[]
        for tsym_irrep_index in 1:num_irreps(tsym)
            rhsr = symmetry_reduce(hsr, IrrepComponent(tsymbed, tsym_irrep_index))
            
            j1_redrep = represent(rhsr, j1)
            j1_redmat = Matrix(j1_redrep)
            @test isapprox(j1_redmat, adjoint(j1_redmat); atol=tol)
            append!(eigenvalues2, eigvals(Hermitian(j1_redmat)))
        end
        sort!(eigenvalues2)
        @test length(eigenvalues1) == length(eigenvalues2)
        @test isapprox(eigenvalues1, eigenvalues2; atol=tol)
    end # testset iterator
    
end # testset RedOp4

@testset "RedOp7" begin
    # testing complex phase
    n = 7
    unitcell = makeunitcell(1.0; SiteType=String)
    addsite!(unitcell, "Spin", FractCoord([0], [0.0]))
    lattice = make_lattice(unitcell, n)
    
    tol = Base.rtoldefault(Float64)
    spin_site = Site([State("Up", 1), State("Dn",-1)])
    
    hs = HilbertSpace(repeat([spin_site], n))
    σ = Dict( (isite, j) => pauli_matrix(hs, isite, j) for isite in 1:n, j in [:x, :y, :z, :+, :-])
    j1 = sum(σ[i, j] * σ[mod(i, n) + 1 , j] for i in 1:n, j in [:x, :y, :z])
    
    hsr = represent(HilbertSpaceSector(hs, 0))
    # translation_group = TranslationGroup([Permutation([2,3,4,5,6,7,1])])
    # @test isinvariant(hs, translation_group, j1)
    # @test isinvariant(HilbertSpaceSector(hs, 0), translation_group, j1)
    
    j1_rep = represent(hsr, j1)
    mat_size = size(j1_rep)
    
    mat0 = Matrix(j1_rep)
    mat1 = zeros(ComplexF64, mat_size)
    mat2 = zeros(ComplexF64, mat_size)
    for i in 1:mat_size[1]
        for (j, v) in get_row_iterator(j1_rep, i)
            if j > 0
                mat1[i,j] += v
            end
        end
    end
    
    for j in 1:mat_size[2]
        for (i, v) in get_column_iterator(j1_rep, j)
            if i > 0
                mat2[i,j] += v
            end
        end
    end
    
    @test isapprox(mat0, mat1; atol=1E-8)
    @test isapprox(mat0, mat2; atol=1E-8)
    
    
end

@testset "ReducedRepresentation-Square" begin
    tol = Base.rtoldefault(Float64)
    shape = [3 0; 0 3]
    unitcell = makeunitcell([1.0 0.0; 0.0 1.0]; SiteType=String)
    addsite!(unitcell, "Spin", FractCoord([0, 0], [0.0, 0.0]))
    lattice = make_lattice(unitcell, shape)
    n = numsite(lattice.supercell)
    
    # Test State and Site
    up = State("Up", 1)
    dn = State("Dn",-1)
    spin_site = Site([up, dn])
    
    hs = HilbertSpace(repeat([spin_site], n))
    σ = Dict( (isite, j) => pauli_matrix(hs, isite, j) for isite in 1:n, j in [:x, :y, :z, :+, :-])
    # j1 = sum(σ[i, j] * σ[mod(i, n) + 1 , j] for i in 1:n, j in [:x, :y, :z])
    # j2 = sum(σ[i, j] * σ[mod(i+1, n) + 1 , j] for i in 1:n, j in [:x, :y, :z])
    
    # hsr = represent(HilbertSpaceSector(hs, 0))
    hsr = represent(hs)
    tsym = FiniteTranslationSymmetry(lattice)
    tsymbed = embed(lattice, tsym)
    psym = PointSymmetryDatabase.get2d(6)
    psymbed = embed(lattice, psym)
    # translation_group = TranslationGroup([Permutation([2,3,4,1])])
    # @test isinvariant(hs, translation_group, j1)
    # @test isinvariant(HilbertSpaceSector(hs, 0), translation_group, j1)
    
    @testset "RHSR" begin
        #rhsr = symmetry_reduce(hsr, translation_group, [0//1])
        for symred in [symmetry_reduce, symmetry_reduce_serial, symmetry_reduce_parallel]
            rhsr = symred(hsr, IrrepComponent(tsymbed, 2, 1))
            @test scalartype(rhsr) === ComplexF64
            @test scalartype(typeof(rhsr)) === ComplexF64
            @test valtype(rhsr) === ComplexF64
            @test valtype(typeof(rhsr)) === ComplexF64
            @test bintype(rhsr) === UInt
            @test bintype(typeof(rhsr)) === UInt
            @test dimension(rhsr) <= 2^n
            @test bitwidth(rhsr) == n
            show(devnull, MIME("text/plain"), rhsr)  # make sure it doesn't crash
            
            rhsr = symred(hsr, IrrepComponent(psymbed, 2, 1))
            @test scalartype(rhsr) === ComplexF64
            @test scalartype(typeof(rhsr)) === ComplexF64
            @test valtype(rhsr) === ComplexF64
            @test valtype(typeof(rhsr)) === ComplexF64
            @test bintype(rhsr) === UInt
            @test bintype(typeof(rhsr)) === UInt
            @test dimension(rhsr) <= 2^n
            @test bitwidth(rhsr) == n
            show(devnull, MIME("text/plain"), rhsr)  # make sure it doesn't crash
        end
    end
    
    @testset "RHSR-new implementation" begin #TODO: rename
        symops_and_amplitudes = collect(get_irrep_iterator(IrrepComponent(tsymbed, 2, 1)))
        for symred in [symmetry_reduce, symmetry_reduce_serial, symmetry_reduce_parallel]
            rhsr = symred(hsr, symops_and_amplitudes)
            @test scalartype(rhsr) === ComplexF64
            @test scalartype(typeof(rhsr)) === ComplexF64
            @test valtype(rhsr) === ComplexF64
            @test valtype(typeof(rhsr)) === ComplexF64
            @test bintype(rhsr) === UInt
            @test bintype(typeof(rhsr)) === UInt
            @test dimension(rhsr) <= 2^n
            @test bitwidth(rhsr) == n
            show(devnull, MIME("text/plain"), rhsr)  # make sure it doesn't crash  
        end
        symops_and_amplitudes = collect(get_irrep_iterator(IrrepComponent(psymbed, 2, 1)))
        for symred in [symmetry_reduce, symmetry_reduce_serial, symmetry_reduce_parallel]
            rhsr = symred(hsr, symops_and_amplitudes)
            @test scalartype(rhsr) === ComplexF64
            @test scalartype(typeof(rhsr)) === ComplexF64
            @test valtype(rhsr) === ComplexF64
            @test valtype(typeof(rhsr)) === ComplexF64
            @test bintype(rhsr) === UInt
            @test bintype(typeof(rhsr)) === UInt
            @test dimension(rhsr) <= 2^n
            @test bitwidth(rhsr) == n
            show(devnull, MIME("text/plain"), rhsr)  # make sure it doesn't crash  
        end
    end
    
    #=
    j1_rep = represent(hsr, j1)
    j1_mat = Matrix(j1_rep)
    @test hsr.basis == [0b0011, 0b0101, 0b0110, 0b1001, 0b1010, 0b1100]
    
    @testset "ROR" begin
        #rhsr = symmetry_reduce(hsr, translation_group, [0//1])
        rhsr = symmetry_reduce(hsr, IrrepComponent(tsymbed, 1, 1))
        @test dimension(rhsr) == 2
        @test rhsr.basis == UInt[0b0011, 0b0101]
        dim = dimension(rhsr)
        
        j1_redrep = represent(rhsr, j1)
        j2_redrep = represent(rhsr, j2)
        show(devnull, MIME("text/plain"), j1_redrep) # make sure it doesn't crash
        
        @testset "typetraits" begin
            @test scalartype(j1) === Complex{Int}
            @test valtype(j1) === Complex{Int}
            
            @test scalartype(j1_redrep) === ComplexF64
            @test scalartype(typeof(j1_redrep)) === ComplexF64
            @test valtype(j1_redrep) === ComplexF64
            @test valtype(typeof(j1_redrep)) === ComplexF64
            
            @test bintype(j1_redrep) === UInt
            @test bintype(typeof(j1_redrep)) === UInt
            
            @test spacetype(j1_redrep) === typeof(rhsr)
            @test operatortype(j1_redrep) === typeof(j1)
            @test spacetype(typeof(j1_redrep)) === typeof(rhsr)
            @test operatortype(typeof(j1_redrep)) === typeof(j1)
            @test get_space(j1_redrep) === rhsr
        end
        
        @testset "property" begin
            @test bitwidth(j1_redrep) == bitwidth(rhsr)
            @test bitwidth(j2_redrep) == bitwidth(rhsr)
            @test dimension(j1_redrep) == dimension(rhsr)
            @test dimension(j2_redrep) == dimension(rhsr)
        end
        
        @testset "sym" begin
            @test ishermitian(represent(hsr, σ[2, :x]))
            @test !ishermitian(represent(hsr, σ[2, :+]))
        end
        
        @testset "unary operator" begin
            @test +j1_redrep == j1_redrep
            @test -j1_redrep == represent(rhsr, -j1)
        end
        
        @testset "binary operator" begin
            @test simplify(j1_redrep + j2_redrep) == represent(rhsr, simplify(j1 + j2))
            @test simplify(j1_redrep - j2_redrep) == represent(rhsr, simplify(j1 - j2))
            @test simplify(j1_redrep * j2_redrep) == represent(rhsr, simplify(j1 * j2))
            
            @test simplify(j1_redrep * 0) == represent(rhsr, NullOperator())
            @test simplify(0 * j1_redrep) == represent(rhsr, NullOperator())
            @test simplify(j1_redrep * 2) == represent(rhsr, simplify(j1 * 2))
            @test simplify(2 * j1_redrep) == represent(rhsr, simplify(2 * j1))
            @test simplify(j1_redrep / 2) == represent(rhsr, simplify(j1 / 2))
            @test simplify(2 \ j1_redrep) == represent(rhsr, simplify(2 \ j1))
        end
        
        let
            psis = [normalize([1.0, 0.0, 1.0, 1.0, 0.0, 1.0]), normalize([0.0, 1.0, 0.0, 0.0, 1.0, 0.0])]
            H = zeros(ComplexF64, (2,2))
            for i in 1:2, j in 1:2
                H[i,j] = psis[i] ⋅ (j1_mat * psis[j])
            end
            @test isapprox(Matrix(j1_redrep), H; atol=tol)
            
            @testset "get_row_iterator" begin
                rowvec = Vector{ComplexF64}(undef, dimension(rhsr))
                for irow_r in 1:dimension(rhsr)
                    fill!(rowvec, zero(ComplexF64))
                    err = zero(ComplexF64)
                    for (icol_r, ampl) in get_row_iterator(j1_redrep, irow_r)
                        if 1 <= icol_r <= dimension(rhsr)
                            rowvec[icol_r] += ampl
                        else
                            err += ampl
                        end
                    end
                    @test isapprox(err, 0; atol=tol) # is this necessarily true?
                    @test isapprox(rowvec, H[irow_r, :]; atol=tol)
                end
            end
            
            @testset "get_column_iterator" begin
                colvec = zeros(ComplexF64, dimension(rhsr))
                for icol_r in 1:dimension(rhsr)
                    colvec[:] .= zero(ComplexF64)
                    
                    colvec[:] .= zero(ComplexF64)
                    err = zero(ComplexF64)
                    for (irow_r, ampl) in get_column_iterator(j1_redrep, icol_r)
                        if 1 <= irow_r <= dimension(rhsr)
                            colvec[irow_r] += ampl
                        else
                            err += ampl
                        end
                    end
                    @test isapprox(err, 0; atol=tol) # is this necessarily true?
                    @test isapprox(colvec, H[:, icol_r]; atol=tol)
                end
            end # testset get_column_iterator
            
            @testset "get_element" begin
                for irow_r in 1:dimension(rhsr)
                    for icol_r in 1:dimension(rhsr)
                        @test isapprox(get_element(j1_redrep, irow_r, icol_r), H[irow_r, icol_r]; atol=tol)
                        @test isapprox(j1_redrep[irow_r, icol_r], H[irow_r, icol_r]; atol=tol)
                    end
                end
            end
            
            @testset "get exceptions" begin
                dim = dimension(rhsr)
                opr = j1_redrep
                @test_throws BoundsError get_row_iterator(opr, 0)
                @test_throws BoundsError get_row_iterator(opr, dim+1)
                @test_throws BoundsError get_column_iterator(opr, 0)
                @test_throws BoundsError get_column_iterator(opr, dim+1)
                @test_throws BoundsError get_row(opr, 0)
                @test_throws BoundsError get_row(opr, dim+1)
                @test_throws BoundsError get_column(opr, 0)
                @test_throws BoundsError get_column(opr, dim+1)
                @test_throws BoundsError get_element(opr, 0, 1)
                @test_throws BoundsError get_element(opr, dim+1, 1)
                @test_throws BoundsError get_element(opr, 1, 0)
                @test_throws BoundsError get_element(opr, 1, dim+1)
                @test_throws BoundsError opr[0, 1]
                @test_throws BoundsError opr[dim+1, 1]
                @test_throws BoundsError opr[1, 0]
                @test_throws BoundsError opr[1, dim+1]
            end
        end
    end
    =#
end


@testset "RedOp-nontriv" begin
    # TODO: at an angle
end
