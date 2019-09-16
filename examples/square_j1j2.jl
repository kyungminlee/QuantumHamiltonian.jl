using LinearAlgebra
using SparseArrays
using Arpack
using JLD2

using ExactDiagonalization

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

function main()
  QN = Int
  up = State{QN}("Up", 1)
  dn = State{QN}("Dn",-1)
  spinsite = Site{QN}([up, dn])

  (n1, n2) = (4, 4)
  n_sites = n1 * n2
  
  hs = AbstractHilbertSpace{QN}()
  for i in 1:n_sites
    add_site!(hs, spinsite)
  end

  @show hs

  PAULI_MATRICES = [ Float64[0 1.0; 1.0 0.0], ComplexF64[0.0 -1.0*im; 1.0*im 0.0], Float64[1.0 0.0; 0.0 -1.0]]

  sigma(i ::Integer, j ::Integer) = KroneckerProductOperator{ComplexF64}(hs, 1.0, Dict(i=>PAULI_MATRICES[j]))
  sigma_plus(i ::Integer) = KroneckerProductOperator{ComplexF64}(hs, 1.0, Dict(i=>[0.0 1.0; 0.0 0.0]))
  sigma_minus(i ::Integer) = KroneckerProductOperator{ComplexF64}(hs, 1.0, Dict(i=>[0.0 0.0; 1.0 0.0]))

  (nearest_neighbor_pairs, second_nearest_neighbor_pairs, chiral_triplets) = make_square_lattice(n1, n2)

  (J1, J2, J3) = (1.0, 1.0, 1.0)

  J2s = 0:0.2:2
  J3s = 0:0.2:2

  j1_terms = KroneckerProductOperator{ComplexF64}[]
  for (i,j) in nearest_neighbor_pairs
    push!(j1_terms, 2 * sigma_plus(i) * sigma_minus(j))
    push!(j1_terms, 2 * sigma_minus(i) * sigma_plus(j))
    push!(j1_terms, sigma(i,3) * sigma(j,3))
  end

  j2_terms = KroneckerProductOperator{ComplexF64}[]
  for (i,j) in second_nearest_neighbor_pairs
    push!(j2_terms, 2 * sigma_plus(i) * sigma_minus(j))
    push!(j2_terms, 2 * sigma_minus(i) * sigma_plus(j))
    push!(j2_terms, sigma(i,3) * sigma(j,3))
  end

  j3_terms = KroneckerProductOperator{ComplexF64}[]
  for (i,j,k) in chiral_triplets
    #push!(j3_terms, sigma(i,1) * sigma(j,2) * sigma(k, 3))
    #push!(j3_terms, sigma(i,2) * sigma(j,3) * sigma(k, 1))
    #push!(j3_terms, sigma(i,3) * sigma(j,1) * sigma(k, 2))
    push!(j3_terms, 2*im * sigma(i,3)* sigma_plus(j) * sigma_minus(k))
    push!(j3_terms,-2*im * sigma(i,3)* sigma_minus(j) * sigma_plus(k))

    push!(j3_terms, 2*im * sigma(j,3)* sigma_plus(k) * sigma_minus(i))
    push!(j3_terms,-2*im * sigma(j,3)* sigma_minus(k) * sigma_plus(i))
    
    push!(j3_terms, 2*im * sigma(k,3)* sigma_plus(i) * sigma_minus(j))
    push!(j3_terms,-2*im * sigma(k,3)* sigma_minus(i) * sigma_plus(j))
    
  end



  sectors = quantum_number_sectors(hs)
  sectors = [x for x in sectors if x >= 0]

  for qn in sectors
    println("------------------------------")
    println("Sector = ", qn)
    println("Concretizing Hilbert Space")
    flush(stdout)
    
    chs = concretize(hs, Set([qn]))
    println("Materializing Terms")
    flush(stdout)
    j1_sparse, ϵ = materialize_parallel(chs, j1_terms)
    @show "J1", ϵ
    j2_sparse, ϵ = materialize_parallel(chs, j2_terms)
    @show "J2", ϵ
    j3_sparse, ϵ = materialize_parallel(chs, j3_terms)
    @show "J3", ϵ

    spectrum = Dict()
    #spectrum ::Dict{Int, Vector{Float64}} = if size(H)[1] <= 20
    for J2 in J2s, J3 in J3s
      @show (J1, J2, J3)

      H = J1 * j1_sparse + J2 * j2_sparse + J3 * j3_sparse
      println("Diagonlizating Hamiltonian")
      flush(stdout)

      if size(H)[1] <= 20
        spectrum[(J1, J2, J3)] = sort(eigvals(Hermitian(Matrix(H))))
      else
        (eigenvalues, eigenvectors) = eigs(H; which=:SR)
        spectrum[(J1, J2, J3)] = sort(real.(eigenvalues))
      end

    end # for J2, J3
    filename = "spectrum_$(n1)_$(n2)_$(qn).jld2"
    @save filename n1, n2, qn, spectrum
  end # for qn
end

main()
