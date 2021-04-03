using SparseArrays
using LinearAlgebra
using QuantumHamiltonian
using LatticeTools
# using MinimalPerfectHash

n_sites = 7;
(hs, σ) = QuantumHamiltonian.Toolkit.spin_half_system(n_sites)

unitcell = make_unitcell(1.0; SiteType=String)
addsite!(unitcell, "Spin", FractCoord([0], [0.0]))
lattice = make_lattice(unitcell, n_sites)
tsym = FiniteTranslationSymmetry(lattice)

Sx = sum(σ(i,:x) for i in 1:n_sites)
Sy = sum(σ(i,:y) for i in 1:n_sites)
Sz = sum(σ(i,:z) for i in 1:n_sites)

spin_squared = simplify( Sx^2 + Sy^2 + Sz^2 )
j1 = sum(σ(i, j) * σ(mod(i, n_sites)+1 , j) for i in 1:n_sites for j in [:x, :y, :z]);
@show j1


hss = HilbertSpaceSector(hs, 5)
hsr = represent_dict(hss);
j1_rep = represent(hsr, j1)

x = rand(Float64, (7, 5))
y = j1_rep * x
j1_mat = Matrix(j1_rep)
y2 = j1_mat * x

@show y - y2