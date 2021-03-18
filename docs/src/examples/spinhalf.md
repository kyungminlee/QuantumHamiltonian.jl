# Example 1: S=1/2 Heisenberg Chain

```@example
using SparseArrays
using LinearAlgebra
using Arpack
using Plots

using LatticeTools
using QuantumHamiltonian

println("# S=1/2 Heisenberg Chain")
n_sites = 8;
(hs, σ) = QuantumHamiltonian.Toolkit.spin_half_system(n_sites)
println("Quantum number sectors (2Sz) : ", quantum_number_sectors(hs))

S = Dict(μ => sum(σ(i, μ) for i in 1:n_sites) for μ in [:x, :y, :z])
spin_squared = simplify(S[:x]^2 + S[:y]^2 + S[:z]^2)
j1 = simplify(sum(σ(i, j) * σ(mod(i, n_sites)+1 , j) for i in 1:n_sites for j in [:x, :y, :z]))
hs_rep = represent(hs);     # Use FrozenSortedArrayIndex{UInt} for basis lookup

println("## All sectors at once")
plt = plot(size=(400, 400))
begin
  j1_rep = represent(hs_rep, j1)
  eigenvalues, eigenvectors = eigs(j1_rep; nev=32, which=:SR, )
  eigenvalues = real.(eigenvalues)
  println("E : ", eigenvalues[1:5])
  scatter!(plt,
           zeros(size(eigenvalues)), eigenvalues,
           markershape=:hline,
           markersize=10,
           markerstrokecolor=:red,
           legend=:none)
end

println("## Sz=0, each momentum sectors")
hs_sector = HilbertSpaceSector(hs, 0)
hs_rep = represent_dict(hs_sector) # Use Dict{UInt, Int} for basis lookup

unitcell = make_unitcell(1.0; SiteType=String)
addsite!(unitcell, "Spin", FractCoord([0], [0.0]))
lattice = make_lattice(unitcell, 8)
tsymbed = translation_symmetry_embedding(lattice)

#translation_group = TranslationGroup([Permutation([ mod(i, n_sites)+1 for i in 1:n_sites])])
ks = symmetry(tsymbed).fractional_momenta

for tsic in get_irrep_components(tsymbed)
  k = ks[tsic.irrep_index]
  hs_redrep = symmetry_reduce(hs_rep, tsic)
  j1_redrep = represent(hs_redrep, j1)
  j1_redrep_dense = Matrix(j1_redrep)   # Make a dense matrix
  eigenvalues = eigvals(Hermitian(j1_redrep_dense))
  println("E(k=", join(string.(k), ","), ") : ", eigenvalues[1:5])
  scatter!(plt,
           ones(size(eigenvalues)).*tsic.irrep_index, eigenvalues,
           markershape=:hline,
           markersize=10,
           markerstrokecolor=:blue,
           legend=:none)
end
xticks!(plt, collect(0:n_sites), ["All", string.(1:n_sites)...])
xlims!(plt, -1, n_sites+1)
xlabel!(plt, "Momentum index")
ylabel!(plt, "Energy")
savefig(plt, "spinchain.svg"); nothing
```

![](spinchain.svg)
