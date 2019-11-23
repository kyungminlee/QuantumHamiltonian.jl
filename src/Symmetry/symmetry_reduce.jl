export symmetry_reduce, symmetry_reduce_serial, symmetry_reduce_parallel
export symmetry_reduce, symmetry_unreduce

import TightBindingLattice.TranslationGroup

import Dates

function symmetry_reduce(
    hsr ::HilbertSpaceRepresentation{QN, BR, DT},
    trans_group ::TranslationGroup,
    fractional_momentum ::AbstractVector{<:Rational};
    ComplexType::DataType=ComplexF64,
    tol::Real=sqrt(eps(Float64))) where {QN, BR, DT}
  symred = Threads.nthreads() == 1 ? symmetry_reduce_serial : symmetry_reduce_parallel
  return symred(hsr, trans_group, fractional_momentum; ComplexType=ComplexType, tol=tol)
end

"""
    symmetry_reduce_serial(hsr, trans_group, frac_momentum; ComplexType=ComplexF64, tol=sqrt(eps(Float64)))

Symmetry-reduce the HilbertSpaceRepresentation using translation group.

"""
function symmetry_reduce_serial(
    hsr ::HilbertSpaceRepresentation{QN, BR, DT},
    trans_group ::TranslationGroup,
    fractional_momentum ::AbstractVector{<:Rational};
    ComplexType::DataType=ComplexF64,
    tol::Real=sqrt(eps(Float64))) where {QN, BR, DT}

  ik = let
    match(k ::Vector{Rational{Int}}) ::Bool = k == fractional_momentum
    findfirst(match, trans_group.fractional_momenta)
  end
  HSR = HilbertSpaceRepresentation{QN, BR, DT}
  ik === nothing && throw(ArgumentError("fractional momentum $(fractional_momentum) not an irrep of the translation group"))

  phases = conj.(trans_group.character_table[ik, :])
  n_basis = length(hsr.basis_list)

  reduced_basis_list = BR[]
  RepresentativeAmplitudeType = NamedTuple{(:representative,:amplitude), Tuple{Int,ComplexType}}
  #representative_amplitude_list = RepresentativeAmplitudeType[(representative=-1,amplitude=zero(ComplexType)) for i in 1:n_basis]
  representative_amplitude_list = Vector{RepresentativeAmplitudeType}(undef, n_basis)
  fill!(representative_amplitude_list, (representative=-1, amplitude=zero(ComplexType)))

  size_estimate = let
    denom = max(1, length(trans_group.fractional_momenta) - 1)
    n_basis ÷ denom
  end

  sizehint!(reduced_basis_list, size_estimate)
  visited = falses(n_basis)

  for ivec_p in 1:n_basis
    visited[ivec_p] && continue
    bvec = hsr.basis_list[ivec_p]

    compatible = true
    ψ = Dict{BR, ComplexType}()

    for i in 1:length(trans_group.elements)
      t = trans_group.translations[i]
      g = trans_group.elements[i]

      bvec_prime = symmetry_apply(hsr.hilbert_space, g, bvec)

      if bvec_prime < bvec
        compatible = false
        break
      elseif bvec_prime == bvec && !is_compatible(fractional_momentum, t)
        compatible = false
        break
      end

      p = phases[i]
      ψ[bvec_prime] = get(ψ, bvec_prime, zero(ComplexType)) + p
    end
    (!compatible) && continue

    choptol!(ψ, tol)
    #@assert !isempty(ψ)

    ivec_p_primes = Int[hsr.basis_lookup[bvec_prime] for bvec_prime in keys(ψ)]

    #@assert !any(visited[ivec_p_primes])
    visited[ivec_p_primes] .= true

    push!(reduced_basis_list, bvec)
    norm = sqrt(sum((abs2.(values(ψ)))))
    for (bvec_prime, amplitude) in ψ
      ivec_p_prime = hsr.basis_lookup[bvec_prime]
      representative_amplitude_list[ivec_p_prime] = (representative=ivec_p, amplitude=amplitude / norm)
    end
  end

  sort!(reduced_basis_list)

  ItemType = NamedTuple{(:index, :amplitude), Tuple{Int, ComplexType}}
  #basis_mapping = ItemType[(index=-1, amplitude=zero(ComplexType)) for b in hsr.basis_list]
  basis_mapping = Vector{ItemType}(undef, length(hsr.basis_list))
  fill!(basis_mapping, (index=-1, amplitude=zero(ComplexType)))
  for (ivec_r, bvec) in enumerate(reduced_basis_list)
    ivec_p = hsr.basis_lookup[bvec]
    amplitude = representative_amplitude_list[ivec_p].amplitude
    basis_mapping[ivec_p] = (index=ivec_r, amplitude=amplitude)
  end

  for (ivec_p_prime, (ivec_p, amplitude)) in enumerate(representative_amplitude_list)
    (ivec_p == -1) && continue  # not in this irrep
    (ivec_p_prime == ivec_p) && continue # already in the lookup
    ivec_r = basis_mapping[ivec_p].index
    basis_mapping[ivec_p_prime] = (index=ivec_r, amplitude=amplitude)
  end

  return ReducedHilbertSpaceRepresentation{HSR, BR, ComplexType}(hsr, trans_group, reduced_basis_list, basis_mapping)
end


function symmetry_reduce_parallel(
    hsr ::HilbertSpaceRepresentation{QN, BR, DT},
    trans_group ::TranslationGroup,
    fractional_momentum ::AbstractVector{<:Rational};
    ComplexType::DataType=ComplexF64,
    tol::Real=sqrt(eps(Float64))) where {QN, BR, DT}

  HSR = HilbertSpaceRepresentation{QN, BR, DT}
  @debug "BEGIN symmetry_reduce_parallel"
  ik = let
    match(k ::Vector{Rational{Int}}) ::Bool  = k == fractional_momentum
    findfirst(match, trans_group.fractional_momenta)
  end

  if ik === nothing
    throw(ArgumentError("fractional momentum $(fractional_momentum)" *
                        " not an irrep of the translation group"))
  end
  phases = conj.(trans_group.character_table[ik, :])
  n_basis = length(hsr.basis_list)
  @debug "Original Hilbert space dimension: $n_basis"

  visit_lock = Threads.SpinLock()
  nthreads = Threads.nthreads()
  local_progress = zeros(Int, nthreads)
  local_reduced_basis_list = [BR[] for i in 1:nthreads]
  RepresentativeAmplitudeType = NamedTuple{(:representative,:amplitude), Tuple{Int,ComplexType}}
  #representative_amplitude_list = RepresentativeAmplitudeType[(representative=-1,amplitude=zero(ComplexType)) for i in 1:n_basis]
  representative_amplitude_list = Vector{RepresentativeAmplitudeType}(undef, n_basis)
  fill!(representative_amplitude_list, (representative=-1, amplitude=zero(ComplexType)))

  size_estimate = let
    denom = max(1, length(trans_group.fractional_momenta) - 1)
    n_basis ÷ denom
  end
  @debug "Estimate for the reduced Hilbert space dimension: $size_estimate"
  for i in eachindex(local_reduced_basis_list)
    sizehint!(local_reduced_basis_list[i], size_estimate ÷ nthreads + 1)
  end

  # Load balancing (the representatives are the smaller binary numbers)
  reorder = Int[]
  sizehint!(reorder, n_basis)
  nblocks = (n_basis + nthreads - 1) ÷ nthreads
  for i in 1:nthreads, j in 1:nblocks
    k = i + nthreads * (j-1)
    if 1 <= k <= n_basis
      push!(reorder, k)
    end
  end
  #@assert sort(reorder) == 1:n_basis

  #prev_progress_time = Dates.now()

  visited = falses(n_basis)
  @debug "Starting reduction (parallel)"
  Threads.@threads for itemp in 1:n_basis
    ivec_p = reorder[itemp]
    visited[ivec_p] && continue
    id = Threads.threadid()
    bvec = hsr.basis_list[ivec_p]

    compatible = true
    ψ = Dict{BR, ComplexType}()

    for i in eachindex(trans_group.elements)
      t = trans_group.translations[i]
      g = trans_group.elements[i]

      bvec_prime = symmetry_apply(hsr.hilbert_space, g, bvec)

      if bvec_prime < bvec
        compatible = false
        break
      elseif bvec_prime == bvec && !is_compatible(fractional_momentum, t)
        compatible = false
        break
      end

      p = phases[i]
      ψ[bvec_prime] = get(ψ, bvec_prime, zero(ComplexType)) + p
    end
    (!compatible) && continue

    choptol!(ψ, tol)
    #@assert !isempty(ψ)

    ivec_p_primes = Int[hsr.basis_lookup[bvec_prime] for bvec_prime in keys(ψ)]

    lock(visit_lock)
    if any(visited[ivec_p_primes]) # check again for thread-safety
      unlock(visit_lock)
      continue
    else
      visited[ivec_p_primes] .= true
      unlock(visit_lock)
      local_progress[id] += length(ivec_p_primes)
    end

    push!(local_reduced_basis_list[id], bvec)

    norm = sqrt( sum( (abs2.(values(ψ))) ) )
    for (bvec_prime, amplitude) in ψ
      ivec_p_prime = hsr.basis_lookup[bvec_prime]
      representative_amplitude_list[ivec_p_prime] = (representative=ivec_p, amplitude=amplitude / norm)
    end
  end
  @debug "Finished reduction (parallel)"

  @debug "Starting garbage collection"
  GC.gc()
  @debug "Finished garbage collection"

  @debug "Collecting basis list"
  reduced_basis_list = BR[]
  sizehint!(reduced_basis_list, sum(length(x) for x in local_reduced_basis_list))
  while !isempty(local_reduced_basis_list)
    lbl = pop!(local_reduced_basis_list)
    append!(reduced_basis_list, lbl)
  end
  @debug "Starting garbage collection"
  GC.gc()
  @debug "Finished garbage collection"

  @debug "Sorting basis list"
  sort!(reduced_basis_list)

  @debug "Starting garbage collection"
  GC.gc()
  @debug "Finished garbage collection"

  ItemType = NamedTuple{(:index, :amplitude), Tuple{Int, ComplexType}}
  #basis_mapping = ItemType[(index=-1, amplitude=zero(ComplexType)) for i in hsr.basis_list]
  basis_mapping = Vector{ItemType}(undef, length(hsr.basis_list))
  fill!(basis_mapping, (index=-1, amplitude=zero(ComplexType)))

  @debug "Collecting basis lookup (diagonal)"
  Threads.@threads for ivec_r in eachindex(reduced_basis_list)
    bvec = reduced_basis_list[ivec_r]
    ivec_p = hsr.basis_lookup[bvec]
    amplitude = representative_amplitude_list[ivec_p].amplitude
    basis_mapping[ivec_p] = (index=ivec_r, amplitude=amplitude)
  end

  @debug "Starting garbage collection"
  GC.gc()
  @debug "Finished garbage collection"

  @debug "Collecting basis lookup (offdiagonal)"
  Threads.@threads for ivec_p_prime in eachindex(representative_amplitude_list)
    ivec_p, amplitude = representative_amplitude_list[ivec_p_prime]
    (ivec_p <= 0) && continue  # not in this irrep
    (ivec_p_prime == ivec_p) && continue # already in the lookup
    bvec = hsr.basis_list[ivec_p]
    ivec_r = basis_mapping[ivec_p].index
    #@assert ivec_r > 0 "ivec_r $ivec_r <= 0"
    basis_mapping[ivec_p_prime] = (index=ivec_r, amplitude=amplitude)
  end

  @debug "Starting garbage collection"
  GC.gc()
  @debug "Finished garbage collection"

  @debug "END symmetry_reduce_parallel"
  return ReducedHilbertSpaceRepresentation{HSR, BR, ComplexType}(hsr, trans_group, reduced_basis_list, basis_mapping)
end



raw"""
    symmetry_unreduce

```math
\begin{pmatrix} l_1 \\ l_2 \\ l_3 \\ \vdots \\ l_n \end{pmatrix}
=
\begin{pmatrix}
. & \cdots & . \\
. & \cdots & . \\
. & \cdots & . \\
  & \dots & \\
. & \cdots &
\end{pmatrix}
\begin{pmatrix} s_1 \\ \vdots \\ s_m \end{pmatrix}
```
"""
function symmetry_unreduce(
    rhsr::ReducedHilbertSpaceRepresentation{HSR, BR, C},
    small_vector::AbstractVector{Si}
  ) where {HSR, BR, C, Si<:Number}
  if length(small_vector) != dimension(rhsr)
    throw(DimensionMismatch("Dimension of the input vector should match the reduced representation"))
  end
  So = promote_type(C, Si)
  large_vector = zeros(So, dimension(rhsr.parent))
  for (i_p, (i_r, ampl)) in enumerate(rhsr.basis_mapping)
    if i_r > 0
      large_vector[i_p] += ampl * small_vector[i_r]
    end
  end
  return large_vector
end

"""
"""
function symmetry_reduce(
    rhsr::ReducedHilbertSpaceRepresentation{HSR, BR, C},
    large_vector::AbstractVector{Si}
  ) where {HSR, BR, C, Si<:Number}
  if length(large_vector) != dimension(rhsr.parent)
    throw(DimensionMismatch("Dimension of the input vector should match the larger representation"))
  end
  So = promote_type(C, Si)
  small_vector = zeros(So, dimension(rhsr))

  # basis mapping
  # (i_p | i_r | ampl) indicates : U_(p, r) = ampl
  for (i_p, (i_r, ampl)) in enumerate(rhsr.basis_mapping)
    if i_r > 0
      small_vector[i_r] += conj(ampl) * large_vector[i_p]
    end
  end
  return small_vector
end
