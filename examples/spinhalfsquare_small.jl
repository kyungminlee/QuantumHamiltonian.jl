using SparseArrays
using LinearAlgebra
using QuantumHamiltonian
using LatticeTools

@show  Threads.nthreads()

## Set up lattice
n1, n2 = 3, 3
n_sites = n1 * n2;

unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; SiteType=String)
addsite!(unitcell, "Spin", FractCoord([0, 0], [0.0, 0.0]))
lattice = make_lattice(unitcell, [n1 0; 0 n2])
tsym = TranslationSymmetry(lattice)
psym = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0])  # inversion symmetry

tsymbed = embed(lattice, tsym)
psymbed = embed(lattice, psym)


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


## Full Solution

hsr = represent(hs)
m = Matrix(represent(hsr, j1))
alleigenvalues1 = eigvals(Hermitian(m))
@show length(alleigenvalues1)


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
        for tsic in get_irrep_components(tsymbed)
        #for tsym_irrep_index in 1:num_irreps(tsym)
            #tsic = TranslationSymmetryIrrepComponent(tsym, tsym_irrep_index, 1)
            rhssr = symmetry_reduce_serial(hssr, tsic)
            rhssr2 = symmetry_reduce_parallel(hssr, tsic)
            #j1_redrep = represent(rhssr, j1)
            #m =  Matrix(j1_redrep)
            dimension(rhssr) == 0 && continue
            # println("- QN: $qn\ttsym: $tsym_irrep_index/$(num_irreps(tsym))\tdimension: $(dimension(rhssr))")
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

    @show length(alleigenvalues2)
    @show length(alleigenvalues3)
    @show norm(alleigenvalues1 - alleigenvalues2)
    @show norm(alleigenvalues1 - alleigenvalues3)
end



if true
    println("## Point Symmetry Only")
    alleigenvalues2 = Float64[]
    alleigenvalues3 = Float64[]
    for qn in quantum_number_sectors(hs)
        hss = HilbertSpaceSector(hs, qn)
        hssr = represent_dict(hss);
        # println("  quantum number: $qn")
        # println("  number of irreps: $(num_irreps(psym))")
        # for psym_irrep_index in 1:num_irreps(psymbed)
        for psic in get_irrep_components(psymbed)
            # println("    psym irrep: $(psym.irreps[psym_irrep_index].name)")
            # println("    psym irrep dimension: $(irrep_dimension(psym, psym_irrep_index))")
            # for psym_irrep_compo in 1:irrep_dimension(psymbed, psym_irrep_index)
                # psic = IrrepComponent(psymbed, psym_irrep_index, psym_irrep_compo)
                rhssr = symmetry_reduce_serial(hssr, psic)
                rhssr2 = symmetry_reduce_parallel(hssr, psic)
                dimension(rhssr) == 0 && continue
                # println("- QN: $qn\tpsym: $psym_irrep_index/$(num_irreps(psym))\t$psym_irrep_compo/$(irrep_dimension(psym, psym_irrep_index))\tdimension:$(dimension(rhssr))")
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
            # end
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
            # # println("    momentum: $(tsym.hypercube.coordinates[tsym_irrep_index])")
            # # println("    little_group: $(psym_little.hermann_mauguin)")
            # # println("    number of irreps: $(num_irreps(psym_little))")
            #     # println("      psym_irrep_index: $psym_irrep_index")
            #     # println("      irrep dimension: $(irrep_dimension(psym_little, psym_irrep_index))")
            for tsic in get_irrep_components(tsymbed)
                psymbed_little = little_symmetry(tsic, psymbed)
                for psic in get_irrep_components(psymbed_little)
                    ssic = SymmorphicIrrepComponent(tsic, psic)
                    rhssr = symmetry_reduce_serial(hssr, ssic)
                    rhssr2 = symmetry_reduce_parallel(hssr, ssic)
                    dimension(rhssr) == 0 && continue
                    # print("- QN: $qn")
                    # print("\ttsym: $tsym_irrep_index/$(num_irreps(tsym))")
                    # print("\tpsym: $(psymbed_little.symmetry.hermann_mauguin)")
                    # print("\t$psym_irrep_index/$(num_irreps(psymbed_little))")
                    # print("\tdimension:$(dimension(rhssr))")
                    # println()
                    # println("        psym_irrepo_compo: $(psym_irrep_compo)")
                    # println("        hilbert dimension: $(dimension(rhssr))")
                    m = Matrix(represent(rhssr, j1))
                    m2 = Matrix(represent(rhssr2, j1))
                    append!(alleigenvalues2, eigvals(Hermitian(m)))
                    append!(alleigenvalues3, eigvals(Hermitian(m2)))
                # end
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

ssymbed = tsymbed ⋊ psymbed

println("## Symmorphic Space Symmetry (S = T ⋊ P)")
let
    alleigenvalues2 = Float64[]
    alleigenvalues3 = Float64[]
    for qn in quantum_number_sectors(hs)
        hss = HilbertSpaceSector(hs, qn)
        hssr = represent_dict(hss);
        # println("  quantum number: $qn")
        # println("  number of tsym_irreps: $(num_irreps(tsym))")
        for ssic in get_irrep_components(ssymbed)
            rhssr = symmetry_reduce_serial(hssr, ssic)
            rhssr2 = symmetry_reduce_parallel(hssr, ssic)
            dimension(rhssr) == 0 && continue
            # print("- QN: $qn")
            # print("\ttsym: $(ssic.normal.irrep_index)/$(num_irreps(tsym))")
            # print("\tpsym: $(symmetry_name(ssic.rest.symmetry))")
            # print("\t$(ssic.rest.irrep_index)/$(num_irreps(ssic.rest.symmetry))")
            # print("\tdimension:$(dimension(rhssr))")
            # println()
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



## Use all three

println("## Translation AND Point AND BitFlip Symmetry")
let
    alleigenvalues2 = Float64[]
    alleigenvalues3 = Float64[]
    for qn in quantum_number_sectors(hs)
        hss = HilbertSpaceSector(hs, qn)
        hssr = represent_dict(hss);
        for tsic in get_irrep_components(tsymbed)
            tsym_irrep = collect(get_irrep_iterator(tsic))
            psymbed_little = little_symmetry(tsic, psymbed)
            for psic in get_irrep_components(psymbed_little)
                psym_irrep = collect(get_irrep_iterator(psic))
                if qn[1] == 0
                    for flip_irrep in [ [(BitFlipSymmetry(false), 1), (BitFlipSymmetry(true), 1)],
                                        [(BitFlipSymmetry(false), 1), (BitFlipSymmetry(true), -1)] ]
                        irrep = make_product_irrep(flip_irrep, tsym_irrep, psym_irrep)                            
                        rhssr = symmetry_reduce_serial(hssr, irrep)
                        rhssr2 = symmetry_reduce_parallel(hssr, irrep)
                        dimension(rhssr) == 0 && continue
                        m = Matrix(represent(rhssr, j1))
                        m2 = Matrix(represent(rhssr2, j1))
                        append!(alleigenvalues2, eigvals(Hermitian(m)))
                        append!(alleigenvalues3, eigvals(Hermitian(m2)))
                    end
                else
                    irrep = make_product_irrep(tsym_irrep, psym_irrep)
                    rhssr = symmetry_reduce_serial(hssr, irrep)
                    rhssr2 = symmetry_reduce_parallel(hssr, irrep)
                    dimension(rhssr) == 0 && continue
                    m = Matrix(represent(rhssr, j1))
                    m2 = Matrix(represent(rhssr2, j1))
                    append!(alleigenvalues2, eigvals(Hermitian(m)))
                    append!(alleigenvalues3, eigvals(Hermitian(m2)))                        
                end
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
