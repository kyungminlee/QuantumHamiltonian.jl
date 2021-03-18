using SparseArrays
using LinearAlgebra
using QuantumHamiltonian
using LatticeTools
using MinimalPerfectHash

@show  Threads.nthreads()

## Set up lattice
n1, n2 = 4, 4
n_sites = n1 * n2;

unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; SiteType=String)
addsite!(unitcell, "Spin", FractCoord([0, 0], [0.0, 0.0]))
lattice = make_lattice(unitcell, [n1 0; 0 n2])
tsym = TranslationSymmetry(lattice)
psym = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0])  # inversion symmetry

## Setup Hilbert Space
(hs, σ) = QuantumHamiltonian.Toolkit.spin_half_system(n_sites)

## Setup Operators
Sx = sum(σ(i,:x) for i in 1:n_sites)
Sy = sum(σ(i,:y) for i in 1:n_sites)
Sz = sum(σ(i,:z) for i in 1:n_sites)

spin_squared = simplify( Sx^2 + Sy^2 + Sz^2 )
j1 = NullOperator()
sub2ind(i1::Integer, i2::Integer) = (i1-1)%n1 + ((i2-1)%n2)*n1 + 1
ind2sub(i::Integer) = (mod(i, n1), mod(i ÷ n1, n2))
for i1 in 1:n1, i2 in 1:n2
    global j1
    j1 += sum(σ(sub2ind(i1, i2), μ) * σ(sub2ind(i1+1, i2), μ) for μ in [:x, :y, :z])
    j1 += sum(σ(sub2ind(i1, i2), μ) * σ(sub2ind(i1, i2+1), μ) for μ in [:x, :y, :z])
end
j1 = simplify(j1)

#j1 = sum(σ(i, j) * σ(mod(i, n_sites)+1 , j) for i in 1:n_sites for j in [:x, :y, :z]);

## Full Solution

# # Full solution too costly
# hsr = represent(hs)
# j1_rep = represent(hsr, j1)
# m = Matrix(j1_rep)
# alleigenvalues1 = eigvals(Hermitian(m))
# @show length(alleigenvalues1)

alleigenvalues1 = nothing

## Using translation symmetry only
println("## Translation Symmetry Only")
if true
    global alleigenvalues1

    alleigenvalues2 = Float64[]
    alleigenvalues3 = Float64[]
    for qn in quantum_number_sectors(hs)
        hss = HilbertSpaceSector(hs, qn)
        hssr = represent_dict(hss)
        # println("  quantum number: $qn")
        # println("  number of irreps: $(num_irreps(tsym))")
        for tsym_irrep_index in 1:num_irreps(tsym)
            tsic = TranslationSymmetryIrrepComponent(tsym, tsym_irrep_index, 1)
            rhssr = symmetry_reduce_serial(hssr, lattice, tsic)
            rhssr2 = symmetry_reduce_parallel(hssr, lattice, tsic)
            #j1_redrep = represent(rhssr, j1)
            #m =  Matrix(j1_redrep)
            dimension(rhssr) == 0 && continue
            println("- QN: $qn\ttsym: $tsym_irrep_index/$(num_irreps(tsym))\tdimension: $(dimension(rhssr))")
            # println("    momentum: $(tsym.hypercube.coordinates[tsym_irrep_index])")
            # println("    hilbert dimension: $(dimension(rhssr))")
            m = Matrix(represent(rhssr, j1))
            m2 = Matrix(represent(rhssr2, j1))
            append!(alleigenvalues2, eigvals(Hermitian(m)))
            append!(alleigenvalues3, eigvals(Hermitian(m2)))
        end
    end
    sort!(alleigenvalues2)
    sort!(alleigenvalues3)
    alleigenvalues1 = alleigenvalues2

    @show length(alleigenvalues2)
    @show length(alleigenvalues3)
    @show norm(alleigenvalues1 - alleigenvalues2)
    @show norm(alleigenvalues1 - alleigenvalues3)
end



if false
    println("## Point Symmetry Only")
    alleigenvalues2 = Float64[]
    alleigenvalues3 = Float64[]
    for qn in quantum_number_sectors(hs)
        hss = HilbertSpaceSector(hs, qn)
        hssr = represent_dict(hss);
        # println("  quantum number: $qn")
        # println("  number of irreps: $(num_irreps(psym))")
        for psym_irrep_index in 1:num_irreps(psym)
            # println("    psym irrep: $(psym.irreps[psym_irrep_index].name)")
            # println("    psym irrep dimension: $(irrep_dimension(psym, psym_irrep_index))")
            for psym_irrep_compo in 1:irrep_dimension(psym, psym_irrep_index)
                psic = PointSymmetryIrrepComponent(psym, psym_irrep_index, psym_irrep_compo)
                rhssr = symmetry_reduce_serial(hssr, lattice, psic)
                rhssr2 = symmetry_reduce_parallel(hssr, lattice, psic)
                dimension(rhssr) == 0 && continue
                println("- QN: $qn\tpsym: $psym_irrep_index/$(num_irreps(psym))\t$psym_irrep_compo/$(irrep_dimension(psym, psym_irrep_index))\tdimension:$(dimension(rhssr))")
                # println("      psym irrep component: $(psym_irrep_compo)")
                # println("      hilbert dimension: $(dimension(rhssr))")
                #
                # for bvec in rhssr.basis_list
                #     println(string(bvec, base=2, pad=9))
                # end
                # println()
                m = Matrix(represent(rhssr, j1))
                m2 = Matrix(represent(rhssr2, j1))
                append!(alleigenvalues2, eigvals(Hermitian(m)))
                append!(alleigenvalues3, eigvals(Hermitian(m2)))
            end
        end
    end
    sort!(alleigenvalues2)
    sort!(alleigenvalues3)

    @show length(alleigenvalues2)
    @show length(alleigenvalues3)
    @show norm(alleigenvalues1 - alleigenvalues2)
    @show norm(alleigenvalues1 - alleigenvalues3)
end


## Use both

println("## Translation AND Point Symmetry")
let
    alleigenvalues2 = Float64[]
    alleigenvalues3 = Float64[]
    for qn in quantum_number_sectors(hs)
        hss = HilbertSpaceSector(hs, qn)
        hssr = represent_dict(hss);
        # println("  quantum number: $qn")
        # println("  number of tsym_irreps: $(num_irreps(tsym))")
        for tsym_irrep_index in 1:num_irreps(tsym)
            tsic = TranslationSymmetryIrrepComponent(tsym, tsym_irrep_index, 1)
            psym_little = little_symmetry(tsic, psym)
            # println("    momentum: $(tsym.hypercube.coordinates[tsym_irrep_index])")
            # println("    little_group: $(psym_little.hermann_mauguin)")
            # println("    number of irreps: $(num_irreps(psym_little))")
            for psym_irrep_index in 1:num_irreps(psym_little)
                # println("      psym_irrep_index: $psym_irrep_index")
                # println("      irrep dimension: $(irrep_dimension(psym_little, psym_irrep_index))")
                for psym_irrep_compo in 1:irrep_dimension(psym_little, psym_irrep_index)
                    psic = PointSymmetryIrrepComponent(psym_little, psym_irrep_index, psym_irrep_compo)
                    ssic = SymmorphicSpaceSymmetryIrrepComponent(tsic, psic)
                    rhssr = symmetry_reduce_serial(hssr, lattice, ssic)
                    rhssr2 = symmetry_reduce_parallel(hssr, lattice, ssic)
                    dimension(rhssr) == 0 && continue


                    print("- QN: $qn")
                    print("\ttsym: $tsym_irrep_index/$(num_irreps(tsym))")
                    print("\tk: $(tsym.hypercube.coordinates[tsym_irrep_index])")
                    print("\tpsym: $(psym_little.hermann_mauguin)")
                    print("\t$psym_irrep_index/$(num_irreps(psym_little))")
                    print("\tdimension:$(dimension(rhssr))")
                    println()
                    # println("        psym_irrepo_compo: $(psym_irrep_compo)")
                    # println("        hilbert dimension: $(dimension(rhssr))")
                    m = Matrix(represent(rhssr, j1))
                    m2 = Matrix(represent(rhssr2, j1))
                    append!(alleigenvalues2, eigvals(Hermitian(m)))
                    append!(alleigenvalues3, eigvals(Hermitian(m2)))
                end
            end # for psym_irrep_index
        end # for tsym_irrep_index
    end
    sort!(alleigenvalues2)
    sort!(alleigenvalues3)

    @show length(alleigenvalues2)
    @show length(alleigenvalues3)
    @show norm(alleigenvalues1 - alleigenvalues2)
    @show norm(alleigenvalues1 - alleigenvalues3)
    #@show alleigenvalues3
end

let
    println("## Symmorphic Space Symmetry")
    alleigenvalues2 = Float64[]
    alleigenvalues3 = Float64[]
    for qn in quantum_number_sectors(hs)
        hss = HilbertSpaceSector(hs, qn)
        hssr = represent_dict(hss);
        # println("  quantum number: $qn")
        # println("  number of tsym_irreps: $(num_irreps(tsym))")
        for ssic in get_irrep_components(lattice, tsym, psym)
            rhssr = symmetry_reduce_serial(hssr, lattice, ssic)
            rhssr2 = symmetry_reduce_parallel(hssr, lattice, ssic)
            dimension(rhssr) == 0 && continue
            print("- QN: $qn")
            print("\ttsym: $(ssic.translation.irrep_index)/$(num_irreps(tsym))")
            print("\tpsym: $(ssic.point.symmetry.hermann_mauguin)")
            print("\t$(ssic.point.irrep_index)/$(num_irreps(ssic.point.symmetry))")
            print("\tdimension:$(dimension(rhssr))")
            println()
            m = Matrix(represent(rhssr, j1))
            m2 = Matrix(represent(rhssr2, j1))
            append!(alleigenvalues2, eigvals(Hermitian(m)))
            append!(alleigenvalues3, eigvals(Hermitian(m2)))
        end # for tsym_irrep_index
    end
    sort!(alleigenvalues2)
    sort!(alleigenvalues3)

    @show length(alleigenvalues2)
    @show length(alleigenvalues3)
    @show norm(alleigenvalues1 - alleigenvalues2)
    @show norm(alleigenvalues1 - alleigenvalues3)
end
