using Test
using SparseArrays
using LinearAlgebra
using QuantumHamiltonian
using LatticeTools
using MinimalPerfectHash
using Printf

@show  Threads.nthreads()

## Set up lattice
n1, n2 = 4, 4
n_sites = n1 * n2;

unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; SiteType=String)
addsite!(unitcell, "Spin", FractCoord([0, 0], [0.0, 0.0]))
lattice = make_lattice(unitcell, [n1 0; 0 n2])
tsym = TranslationSymmetry(lattice)
psym = project(PointSymmetryDatabase.get(13), [1 0 0; 0 1 0])  # inversion symmetry

tsymbed = embed(lattice, tsym)
psymbed = embed(lattice, psym)
ssymbed = tsymbed ⋊ psymbed

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

@show isinvariant(hs, tsymbed, j1)
@show isinvariant(hs, psymbed, j1)

@show isinvariant(hs, ssymbed, j1)

# println("## Symmetry analysis")
# println("number of tsym irreps: $(num_irreps(tsym))")
# for tsic in get_irrep_components(tsymbed)
#     psymbed_little = little_symmetry(tsic, psymbed)
# println("- momentum: $(tsym.fractional_momenta[tsic.irrep_index])")
# println("  little_group: $(symmetry(psymbed_little).hermann_mauguin)")
# println("  number of psym irreps: $(num_irreps(psymbed_little))")
# println("  psym_irrep_components: ",
#               join(["[$(psic.irrep_index), $(psic.irrep_component)]"
#                     for psic in get_irrep_components(psymbed_little)], ", "))
# end


# function foo(ssymbed)
#     parent = Dict()
#     for ssic in get_irrep_components(ssymbed)
#         i1 = (ssic.normal.irrep_index, 0)
#         i2 = (ssic.normal.irrep_index, ssic.rest.irrep_index)
#         parent[i2] = i1
#     end

#     nodes = sort(unique(union(collect(keys(parent)), collect(values(parent)))))
#     node_lookup = Dict(k => "node$i" for (i, k) in enumerate(nodes))

#     for node in nodes
#         if !haskey(parent, node)
#             c = node_lookup[node]
#             println("$c -> node0")
#         end
#     end

#     for (k, v) in parent
#         c = node_lookup[k]
#         p = node_lookup[v]
#         println("$c -> $p")
#     end
# end

# function breakup(ssymbed)

#     foo = Dict()

#     nodes = []
#     for ssic in get_irrep_components(ssymbed)
#         ssic.normal.irrep_index
#         ssic.normal.irrep_component
#         k1 = (ssic.normal.irrep_index, ssic.normal.irrep_component)

#         if k1 ∉ nodes
#             push!(nodes, k1)
#         end
#         if !haskey(foo, k1)
#             foo[k1] = Dict("irrep_index" => k1[1], "irrep_component" => k1[2],
#                            "little_group" => symmetry(ssic.rest.symmetry).hermann_mauguin,
#                            "children"=> [])
#         end
#         push!(foo[k1]["children"],
#                 Dict("irrep_index" => ssic.rest.irrep_index,
#                      "irrep_component" => ssic.rest.irrep_component))


#     end
#     bar = [foo[k] for k in nodes]
#     return bar
# end

# using JSON

# println(JSON.json(breakup(ssymbed)))
# exit()



function myshow(io::IO, ssic::SymmorphicIrrepComponent)
    tsymbed = ssic.normal.symmetry
    psymbed = ssic.rest.symmetry
    tsym = symmetry(tsymbed)
    psym = symmetry(psymbed)

    kf = tsym.fractional_momenta[ssic.normal.irrep_index]
    k = tsymbed.lattice.unitcell.reducedreciprocallatticevectors * kf
    @printf(io, "%8d%8d", ssic.normal.irrep_index, num_irreps(tsym))
    @printf(io, "%16s%8s%8d", k, psym.hermann_mauguin, num_irreps(psym))
    @printf(io, "%8d%8d\n", ssic.rest.irrep_index, ssic.rest.irrep_component)
end

@printf(stdout, "%8s%8s%16s%8s%8s\t%s\n", "ITI", "NTI", "momentum", "LG", "NPI", "IPI/CPI")
for ssic in get_irrep_components(ssymbed)
    myshow(stdout, ssic)
end


# exit()

## Full Solution

# Full solution too costly
if false
    hsr = represent(hs)
    m = Matrix(represent(hsr, j1))
    alleigenvalues1 = eigvals(Hermitian(m))
    @show length(alleigenvalues1)
else
    alleigenvalues1 = nothing
end



## Using translation symmetry only

println("## Translation Symmetry Only")
if true
    global alleigenvalues1

    alleigenvalues2 = Float64[]
    alleigenvalues3 = Float64[]
    for qn in quantum_number_sectors(hs)
        hss = HilbertSpaceSector(hs, qn)
        hssr = represent_dict(hss)
        for tsic in get_irrep_components(tsymbed)
            rhssr = symmetry_reduce_serial(hssr, tsic)
            rhssr2 = symmetry_reduce_parallel(hssr, tsic)
            dimension(rhssr) == 0 && continue
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
        for psic in get_irrep_components(psymbed)
                rhssr = symmetry_reduce_serial(hssr, psic)
                rhssr2 = symmetry_reduce_parallel(hssr, psic)
                dimension(rhssr) == 0 && continue
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

let
    println("## Translation AND Point Symmetry")
    alleigenvalues2 = Float64[]
    alleigenvalues3 = Float64[]
    for qn in quantum_number_sectors(hs)
        hss = HilbertSpaceSector(hs, qn)
        hssr = represent_dict(hss);
        for tsic in get_irrep_components(tsymbed)
            psymbed_little = little_symmetry(tsic, psymbed)
            for psic in get_irrep_components(psymbed_little)
                ssic = SymmorphicIrrepComponent(tsic, psic)
                rhssr = symmetry_reduce_serial(hssr, ssic)
                rhssr2 = symmetry_reduce_parallel(hssr, ssic)
                dimension(rhssr) == 0 && continue
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

ssymbed = tsymbed ⋊ psymbed

println("## Symmorphic Space Symmetry (S = T ⋊ P)")
let
    alleigenvalues2 = Float64[]
    alleigenvalues3 = Float64[]
    for qn in quantum_number_sectors(hs)
        hss = HilbertSpaceSector(hs, qn)
        hssr = represent_dict(hss);
        for ssic in get_irrep_components(ssymbed)
            rhssr = symmetry_reduce_serial(hssr, ssic)
            rhssr2 = symmetry_reduce_parallel(hssr, ssic)
            dimension(rhssr) == 0 && continue
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
