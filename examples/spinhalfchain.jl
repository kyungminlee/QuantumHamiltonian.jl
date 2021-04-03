using SparseArrays
using LinearAlgebra
using QuantumHamiltonian
using LatticeTools
using MinimalPerfectHash

n_sites = 7;
(hs, σ) = QuantumHamiltonian.Toolkit.spin_half_system(n_sites)

unitcell = makeunitcell(1.0; SiteType=String)
addsite!(unitcell, "Spin", FractCoord([0], [0.0]))
lattice = make_lattice(unitcell, n_sites)
tsym = FiniteTranslationSymmetry(lattice)

Sx = sum(σ(i,:x) for i in 1:n_sites)
Sy = sum(σ(i,:y) for i in 1:n_sites)
Sz = sum(σ(i,:z) for i in 1:n_sites)

spin_squared = simplify( Sx^2 + Sy^2 + Sz^2 )
j1 = sum(σ(i, j) * σ(mod(i, n_sites)+1 , j) for i in 1:n_sites for j in [:x, :y, :z]);
@show j1

for qn in quantum_number_sectors(hs)
hss = HilbertSpaceSector(hs, qn)
hsr = represent_dict(hss);

let j1_rep = represent(represent(hs), j1)
    @show eigvals(Hermitian(Matrix(j1_rep)))
end

let
    rhsr = symmetry_reduce_serial(hsr, lattice, tsym, 2)
    @show tsym.orthogonal_coordinates[2]
    j1_redrep = represent(rhsr, j1)
    m =  Matrix(j1_redrep)
    println("m = ")
    display(m)
end
exit()



for tsym_irrep_index in 1:num_irreps(tsym)
    # @show tsym_irrep_index
    rhsr = symmetry_reduce_serial(hsr, lattice, tsym, tsym_irrep_index)
    # @show rhsr
    # @show rhsr.parent
    # @show rhsr.basis_list
    # @show rhsr.basis_mapping_index
    # @show rhsr.basis_mapping_amplitude
    # # println()
    j1_redrep = represent(rhsr, j1)
    # @show typeof(j1_redrep), j1_redrep
    m =  Matrix(j1_redrep)
    @show eigvals(Hermitian(m))
end

psym = project(PointSymmetryDatabase.get(2), [1 0 0;])  # inversion symmetry

@show psym.hermann_mauguin

@show iscompatible(tsym, psym)
for tsym_irrep_index in 1:num_irreps(tsym)
    psym_little = little_symmetry(tsym, tsym_irrep_index, psym)
    @assert iscompatible(tsym, tsym_irrep_index, psym_little)
    @show psym_little.hermann_mauguin


end




#=
translation_group = TranslationGroup([Permutation([ mod(i, n_sites)+1 for i in 1:n_sites])])
ks = translation_group.fractional_momenta
rhsr = symmetry_reduce(hsr, translation_group, ks[1])

j1_redrep = represent(rhsr, j1)
j1_redrep_sparse = sparse(j1_redrep)
=#




#
# using BenchmarkTools
# using Arpack
#
#
#
#
# function simplify_combine(so::SumOperator{S, BR}; tol::Real=Base.rtoldefault(Float64)) where {S, BR}
#   bw = sizeof(BR)*8
#   diagonal_terms = [t for t in so.terms if t.bitrow == t.bitcol]
#   offdiagonal_terms = [t for t in so.terms if t.bitrow != t.bitcol]
#
#   @label simplify_combine_loop_start
#
#   for ib in 0:(bw-1)
#     mask = make_bitmask(ib+1, ib)
#
#     for (i1, t1) in enumerate(diagonal_terms)
#       (t1.bitmask & mask) == 0x0 && continue
#       t1.bitrow != t1.bitcol && continue
#
#       for i2 in (i1+1):length(diagonal_terms)
#         t2 = diagonal_terms[i2]
#         (t2.bitmask & mask) == 0x0 && continue
#         t2.bitrow != t2.bitcol && continue
#
#         if ( ((t1.bitmask & ~mask) == (t2.bitmask & ~mask)) &&
#              ((t1.bitrow  & ~mask) == (t2.bitrow  & ~mask)) &&
#              ((t1.bitrow  &  mask) != (t2.bitrow  &  mask)) )
#
#           if isapprox(t1.amplitude, t2.amplitude; atol=tol)
#             i_small, i_large = i1, i2
#             new_bitmask = t1.bitmask & (~mask)
#             new_bitrow  = t1.bitrow  & (~mask)
#             new_bitcol  = t1.bitcol  & (~mask)
#             new_amplitude = t1.amplitude
#
#             diagonal_terms[i_small] = PureOperator{S, BR}(new_bitmask, new_bitrow, new_bitcol, new_amplitude)
#             deleteat!(diagonal_terms, i_large)
#
#             @goto simplify_combine_loop_start
#             #dirty = true
#             #break
#           else
#             # different amplitude. do nothing
#           end # if amplitude
#         end # if bits
#
#       end # for i2
#     end # for i1
#   end # for ib
#
#   #end # while dirty
#   return simplify(SumOperator{S, BR}([diagonal_terms..., offdiagonal_terms...]))
# end
#
#
# spin_squared_2 = simplify_combine(spin_squared)
#
# # eigs(j1_redrep)
# # @btime eigs(j1_redrep)
# # @timev eigs(j1_redrep)
