using StaticArrays
using DataStructures
using LinearAlgebra
using SparseArrays

include("util.jl")
include("HilbertSpace.jl")
include("Operator.jl")

function make_square_lattice(n1 ::Integer, n2 ::Integer)
  n_sites = n1 * n2
  
  nearest_neighbor_bond_types = [(1,0), (0,1)]
  nearest_neighbor_pairs = Tuple{Int, Int}[]
  for i1 in 1:n1
      for i2 in 1:n2
          i = (i1-1) * n2 + (i2-1) + 1
          for (d1, d2) in nearest_neighbor_bond_types
              j1 = mod(i1 + d1 - 1, n1) + 1
              j2 = mod(i2 + d2 - 1, n2) + 1
              j = (j1-1) * n2 + (j2-1) + 1

              push!(nearest_neighbor_pairs, (i,j))
          end
      end
  end

  second_nearest_neighbor_bond_types = [(1,1), (-1,1)]
  second_nearest_neighbor_pairs = Tuple{Int, Int}[]
  for i1 in 1:n1
      for i2 in 1:n2
          i = (i1-1) * n2 + (i2-1) + 1
          for (d1, d2) in second_nearest_neighbor_bond_types
              j1 = mod(i1 + d1 - 1, n1) + 1
              j2 = mod(i2 + d2 - 1, n2) + 1
              j = (j1-1) * n2 + (j2-1) + 1

              push!(second_nearest_neighbor_pairs, (i,j))
          end
      end
  end

  #=
        o
      / | \
    o - o - o 
      \ | /
        o
  =#
  chiral_triplet_types = [
      (( 1, 0), ( 0, 1)),
      (( 0, 1), (-1, 0)),
      ((-1, 0), ( 0,-1)),
      (( 0,-1), ( 1, 0)),
  ]
  chiral_triplets = Tuple{Int64, Int64, Int64}[]
  for i1 in 1:n1, i2 in 1:n2
    i = (i1-1) * n2 + (i2-1) + 1
    for ((d1, d2), (e1, e2)) in chiral_triplet_types
      j1 = mod(i1 + d1 - 1, n1) + 1
      j2 = mod(i2 + d2 - 1, n2) + 1
      j = (j1-1) * n2 + (j2-1) + 1

      k1 = mod(i1 + e1 - 1, n1) + 1
      k2 = mod(i2 + e2 - 1, n2) + 1
      k = (k1-1) * n2 + (k2-1) + 1
      push!(chiral_triplets, (i,j,k))
    end
  end
  return (nearest_neighbor_pairs, second_nearest_neighbor_pairs, chiral_triplets)
end

#function main()
  QN = Int
  up = State{QN}("Up", 1)
  dn = State{QN}("Dn",-1)
  spinsite = Site{QN}([up, dn])

  (n1, n2) = (2, 2)
  n_sites = n1 * n2
  
  hs = AbstractHilbertSpace{QN}()
  for i in 1:n_sites
    add_site!(hs, spinsite)
  end

  @show hs

  # for basis in hsb.basis_list
  #     @show string(basis, base=2, pad=hs.bitoffsets[end])
  # end
  # for (k, v) in hsb.basis_lookup
  #     @show string(k, base=2), v        
  # end

  PAULI_MATRICES = [ Float64[0 1.0; 1.0 0.0], ComplexF64[0.0 -1.0*im; 1.0*im 0.0], Float64[1.0 0.0; 0.0 -1.0]]

  sigma(i ::Integer, j ::Integer) = KroneckerProductOperator{ComplexF64}(hs, 1.0, Dict(i=>PAULI_MATRICES[j]))
  sigma_plus(i ::Integer) = KroneckerProductOperator{ComplexF64}(hs, 1.0, Dict(i=>[0.0 1.0; 0.0 0.0]))
  sigma_minus(i ::Integer) = KroneckerProductOperator{ComplexF64}(hs, 1.0, Dict(i=>[0.0 0.0; 1.0 0.0]))

  (nearest_neighbor_pairs, second_nearest_neighbor_pairs, chiral_triplets) = make_square_lattice(n1, n2)

  (J1, J2, J3) = (1.0, 1.0, 1.0)

  j1_terms = KroneckerProductOperator{ComplexF64}[]
  for (i,j) in nearest_neighbor_pairs
    push!(j1_terms, 2 * sigma_plus(i) * sigma_minus(j))
    push!(j1_terms, 2 * sigma_minus(i) * sigma_plus(j))
    #push!(j1_terms, J1 * sigma(i,1) * sigma(j,1))
    #push!(j1_terms, J1 * sigma(i,2) * sigma(j,2))
    push!(j1_terms, J1 * sigma(i,3) * sigma(j,3))
  end

  j2_terms = KroneckerProductOperator{ComplexF64}[]
  for (i,j) in second_nearest_neighbor_pairs
    push!(j2_terms, J2 * 2 * sigma_plus(i) * sigma_minus(j))
    push!(j2_terms, J2 * 2 * sigma_minus(i) * sigma_plus(j))
    push!(j2_terms, J2 * sigma(i,3) * sigma(j,3))
  end

  j3_terms = KroneckerProductOperator{ComplexF64}[]
  for (i,j,k) in chiral_triplets
    push!(j3_terms, J3 * sigma(i,1) * sigma(j,2) * sigma(k, 3))
    push!(j3_terms, J3 * sigma(i,2) * sigma(j,3) * sigma(k, 1))
    push!(j3_terms, J3 * sigma(i,3) * sigma(j,1) * sigma(k, 2))
  end

  #for op in j2_terms
  #  @show op.operators
  #end
  #println("------")
  
  # @show apply(hs, 1.0 * sigma(1, 1), UInt64(0x0))
  # @show apply(hs, 1.0 * sigma(2, 1), UInt64(0x0))
  # @show apply(hs, 1.0 * sigma(3, 1), UInt64(0x0))
  # @show apply(hs, 1.0 * sigma(4, 1), UInt64(0x0))

  # @show apply(hs, 1.0 * sigma(1, 1), UInt64(0x1))
  # @show apply(hs, 1.0 * sigma(2, 1), UInt64(0x1))
  # @show apply(hs, 1.0 * sigma(3, 1), UInt64(0x1))
  # @show apply(hs, 1.0 * sigma(4, 1), UInt64(0x1))

  # @show apply(hs, [1.0 * sigma(1, 1), 1.0*sigma(3,1)], UInt64(0x0))

  # return

  let
    chs = concretize(hs)
    H, ε = materialize(chs, j1_terms)
    Hd = Matrix(real.(H))
    println(Hd)
    #@show Matrix(H)
    #@show size(H)
    #(eigenvalues, eigenvectors) = eigs(H; which=:SR)
    #@show sort(eigenvalues)
  end
  return

  sectors = quantum_number_sectors(hs)
  for qn in sectors
    println("------------------------------")
    println("Sector = ", qn)
    chs = concretize(hs, qn)
    # for phi in chs.basis_list
    #   println("phi = ", string(phi, base=2, pad=4))
    #   for isite in 1:4
    #     @show isite, get_state_index(hs, isite, phi)
    #   end
    # end
    H, ε = materialize(chs, j1_terms)
    #@show Matrix(H)
    #continue
    @show size(H)
    @show ε
    if size(H)[1] <= 10
      @show eigvals(Hermitian(Matrix(H)))
    else
      (eigenvalues, eigenvectors) = eigs(H; which=:SR)
      @show sort(real.(eigenvalues))
    end
  end


  # let
  #   chs = concretize(hs)
  #   H, ε = materialize(chs, j1_terms)
  #   @show size(H)
  #   (eigenvalues, eigenvectors) = eigs(H; which=:SR)
  #   @show sort(eigenvalues)

  # end


#end

#main()