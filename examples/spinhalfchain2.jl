using SparseArrays
using LinearAlgebra
using QuantumHamiltonian
using LatticeTools
using MinimalPerfectHash

## Set up lattice
n_sites = 6;

unitcell = makeunitcell(1.0; SiteType=String)
addsite!(unitcell, "Spin", FractCoord([0], [0.0]))
lattice = make_lattice(unitcell, n_sites)
tsym = FiniteTranslationSymmetry(lattice)
psym = project(PointSymmetryDatabase.get(2), [1 0 0;])  # inversion symmetry

## Setup Hilbert Space
(hs, σ) = QuantumHamiltonian.Toolkit.spin_half_system(n_sites)

## Setup Operators
Sx = sum(σ(i,:x) for i in 1:n_sites)
Sy = sum(σ(i,:y) for i in 1:n_sites)
Sz = sum(σ(i,:z) for i in 1:n_sites)

spin_squared = simplify( Sx^2 + Sy^2 + Sz^2 )
j1 = sum(σ(i, j) * σ(mod(i, n_sites)+1 , j) for i in 1:n_sites for j in [:x, :y, :z]);

## Full Solution
hsr = represent(hs)
j1_rep = represent(hsr, j1)
m = Matrix(j1_rep)
alleigenvalues1 = eigvals(Hermitian(m))
@show alleigenvalues1

## Using translation symmetry only
let
    alleigenvalues2 = Float64[]
    alleigenvalues3 = Float64[]
    for qn in quantum_number_sectors(hs)
        hss = HilbertSpaceSector(hs, qn)
        hssr = represent_dict(hss);

        for tsym_irrep_index in 1:num_irreps(tsym)

            tsic = TranslationSymmetryIrrepComponent(tsym, tsym_irrep_index, 1)

            rhssr = symmetry_reduce_serial(hssr, lattice, tsic)
            rhssr2 = symmetry_reduce_parallel(hssr, lattice, tsic)
            #j1_redrep = represent(rhssr, j1)
            #m =  Matrix(j1_redrep)
            if dimension(rhssr) == 0
                continue
            end
            m = Matrix(represent(rhssr, j1))
            m2 = Matrix(represent(rhssr2, j1))
            append!(alleigenvalues2, eigvals(Hermitian(m)))
            append!(alleigenvalues3, eigvals(Hermitian(m2)))
        end
    end
    sort!(alleigenvalues2)
    sort!(alleigenvalues3)

    @show alleigenvalues2
    @show alleigenvalues3
    @show norm(alleigenvalues1 - alleigenvalues2)
    @show norm(alleigenvalues1 - alleigenvalues3)
end


let
    alleigenvalues2 = Float64[]
    for qn in quantum_number_sectors(hs)
        hss = HilbertSpaceSector(hs, qn)
        hssr = represent_dict(hss);
        for psym_irrep_index in 1:num_irreps(psym)
            @show irrep_dimension(psym, psym_irrep_index)
            for psym_irrep_compo in 1:irrep_dimension(psym, psym_irrep_index)
                psic = PointSymmetryIrrepComponent(psym, psym_irrep_index, psym_irrep_compo)
                rhssr = symmetry_reduce_serial(hssr, lattice, psic)
                #rhssr2 = symmetry_reduce_parallel(hssr, lattice, tsic)
                #j1_redrep = represent(rhssr, j1)
                #m =  Matrix(j1_redrep)
                if dimension(rhssr) == 0
                    continue
                end
                m = Matrix(represent(rhssr, j1))
                #m2 = Matrix(represent(rhssr2, j1))
                append!(alleigenvalues2, eigvals(Hermitian(m)))
                #append!(alleigenvalues3, eigvals(Hermitian(m2)))
            end
        end
    end
    sort!(alleigenvalues2)
    #sort!(alleigenvalues3)

    @show alleigenvalues2
    @show norm(alleigenvalues1 - alleigenvalues2)
    #@show alleigenvalues3
end


## Use both

let
    alleigenvalues2 = Float64[]
    for qn in quantum_number_sectors(hs)
        hss = HilbertSpaceSector(hs, qn)
        hssr = represent_dict(hss);

        for tsym_irrep_index in 1:num_irreps(tsym)
            tsic = TranslationSymmetryIrrepComponent(tsym, tsym_irrep_index, 1)
            psym_little = little_symmetry(tsic, psym)
            @show num_irreps(psym_little)
            for psym_irrep_index in 1:num_irreps(psym_little)
                #@show irrep_dimension(psym_little, psym_irrep_index)
                for psym_irrep_compo in 1:irrep_dimension(psym_little, psym_irrep_index)
                    psic = PointSymmetryIrrepComponent(psym_little, psym_irrep_index, psym_irrep_compo)
                    ssic = SymmorphicSpaceSymmetryIrrepComponent(tsic, psic)
                    rhssr = symmetry_reduce_serial(hssr, lattice, ssic)

                    #rhssr2 = symmetry_reduce_parallel(hssr, lattice, tsic)
                    #j1_redrep = represent(rhssr, j1)
                    #m =  Matrix(j1_redrep)
                    if dimension(rhssr) == 0
                        continue
                    end
                    m = Matrix(represent(rhssr, j1))
                    #m2 = Matrix(represent(rhssr2, j1))
                    append!(alleigenvalues2, eigvals(Hermitian(m)))
                    #append!(alleigenvalues3, eigvals(Hermitian(m2)))
                end
            end # for psym_irrep_index
        end # for tsym_irrep_index
    end
    sort!(alleigenvalues2)
    #sort!(alleigenvalues3)

    @show alleigenvalues2
    @show norm(alleigenvalues1 - alleigenvalues2)
    #@show alleigenvalues3
end
