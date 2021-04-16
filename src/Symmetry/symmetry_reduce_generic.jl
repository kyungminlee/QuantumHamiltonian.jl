export symmetry_reduce_serial, symmetry_reduce_parallel
export make_symmetrizer

using LatticeTools


function symmetry_reduce(
    hsr::HilbertSpaceRepresentation{QN, BR, DT},
    symops_and_amplitudes::AbstractArray{Tuple{OperationType, InputScalarType}};
    tol::Real=Base.rtoldefault(float(real(InputScalarType)))
) where {QN, BR, DT, OperationType<:AbstractSymmetryOperation, InputScalarType<:Number}
    symred = Threads.nthreads() == 1 ? symmetry_reduce_serial : symmetry_reduce_parallel
    return symred(hsr, symops_and_amplitudes; tol=tol)
end


function symmetry_reduce(
    hsr::HilbertSpaceRepresentation{QN, BR, DT},
    symops::AbstractArray{OperationType},
    amplitudes::AbstractArray{InputScalarType};
    tol::Real=Base.rtoldefault(float(real(InputScalarType)))
) where {QN, BR, DT, OperationType<:AbstractSymmetryOperation, InputScalarType<:Number}
    symred = Threads.nthreads() == 1 ? symmetry_reduce_serial : symmetry_reduce_parallel
    return symred(hsr, symops, amplitudes; tol=tol)
end


"""
    symmetry_reduce_serial(hilbert_space_representation, symops_and_amplitudes; tol=√ϵ)

The irreps have to follow certain order.
"""
function symmetry_reduce_serial(
    hsr::HilbertSpaceRepresentation{QN, BR, DT},
    symops_and_amplitudes::AbstractArray{Tuple{OperationType, InputScalarType}};
    tol::Real=Base.rtoldefault(float(real(InputScalarType)))
) where {QN, BR, DT, OperationType<:AbstractSymmetryOperation, InputScalarType<:Number}
    symops = [x for (x, y) in symops_and_amplitudes]
    amplitudes = [y for (x, y) in symops_and_amplitudes]
    return symmetry_reduce_serial(hsr, symops, amplitudes; tol=tol)
end


function symmetry_reduce_serial(
    hsr::HilbertSpaceRepresentation{QN, BR, DT},
    symops::AbstractArray{OperationType},
    amplitudes::AbstractArray{InputScalarType};
    tol::Real=Base.rtoldefault(float(real(InputScalarType)))
) where {QN, BR, DT, OperationType<:AbstractSymmetryOperation, InputScalarType<:Number}
    if length(symops) != length(amplitudes)
        throw(ArgumentError("number of symmetry operations and number of amplitudes should match ($(length(symops)) != $(length(amplitudes)))"))
    end
    if !all(isapprox(abs(y), one(abs(y))) || isapprox(abs(y), zero(abs(y))) for y in amplitudes)
        throw(ArgumentError("all amplitudes need to have norm 1 or 0"))
    end
    ScalarType = float(InputScalarType)

    n_basis = length(hsr.basis_list)

    basis_mapping_representative = Vector{Int}(undef, n_basis)
    fill!(basis_mapping_representative, -1)
    basis_mapping_amplitude = zeros(ScalarType, n_basis)
    subgroup_size = length(symops)
    size_estimate = n_basis ÷ max(1, subgroup_size - 1)

    inv_norm_cache = [inv(sqrt(float(i))) for i in 1:subgroup_size]

    reduced_basis_list = BR[]
    sizehint!(reduced_basis_list, size_estimate)

    visited = falses(n_basis)

    basis_states = Vector{BR}(undef, subgroup_size)  # recomputed for every bvec
    basis_phases = ones(ScalarType, subgroup_size)
    basis_amplitudes = Dict{BR, ScalarType}()
    sizehint!(basis_amplitudes, subgroup_size + subgroup_size ÷ 2)

    for ivec_p in 1:n_basis
        visited[ivec_p] && continue
        bvec = hsr.basis_list[ivec_p]

        compatible = true
        for i in 2:subgroup_size
            symop, ampl = symops[i], amplitudes[i]
            bvec_prime, sgn = symmetry_apply(hsr.hilbert_space, symop, bvec)
            if bvec_prime < bvec
                compatible = false
                break
            elseif bvec_prime == bvec && !isapprox(conj(ampl) * sgn, one(ScalarType); atol=tol)
                compatible = false
                break
            end
            basis_states[i] = bvec_prime
            basis_phases[i] = conj(ampl) * sgn
        end # for i
        (!compatible) && continue
        basis_states[1] = bvec

        push!(reduced_basis_list, bvec)

        empty!(basis_amplitudes)
        for i in 1:subgroup_size
            bvec_prime = basis_states[i]
            basis_amplitudes[bvec_prime] = basis_phases[i]
        end
        inv_norm = inv_norm_cache[length(basis_amplitudes)]
        for (bvec_prime, amplitude) in basis_amplitudes
            ivec_p_prime = hsr.basis_lookup[bvec_prime]
            visited[ivec_p_prime] = true
            basis_mapping_representative[ivec_p_prime] = ivec_p
            basis_mapping_amplitude[ivec_p_prime] = amplitude * inv_norm
        end
    end # for ivec_p

    basis_mapping_index = Vector{Int}(undef, n_basis)
    fill!(basis_mapping_index, -1)

    for (ivec_r, bvec) in enumerate(reduced_basis_list)
        ivec_p = hsr.basis_lookup[bvec]
        basis_mapping_index[ivec_p] = ivec_r
    end

    for (ivec_p_prime, ivec_p) in enumerate(basis_mapping_representative)
        (ivec_p <= 0) && continue  # not in this irrep
        (ivec_p_prime == ivec_p) && continue  # already in the lookup
        ivec_r = basis_mapping_index[ivec_p]
        basis_mapping_index[ivec_p_prime] = ivec_r
    end

    HSR = HilbertSpaceRepresentation{QN, BR, DT}
    RHSR = ReducedHilbertSpaceRepresentation{HSR, BR, ScalarType}
    return RHSR(hsr, reduced_basis_list, basis_mapping_index, basis_mapping_amplitude)
end


"""
    symmetry_reduce_parallel(hsr, symops_and_amplitudes; tol=√ϵ)

Symmetry-reduce the HilbertSpaceRepresentation using translation group (multi-threaded).

"""
function symmetry_reduce_parallel(
    hsr::HilbertSpaceRepresentation{QN, BR, DT},
    symops_and_amplitudes::AbstractArray{Tuple{OperationType, InputScalarType}};
    tol::Real=Base.rtoldefault(float(real(InputScalarType)))
) where {QN, BR, DT, OperationType<:AbstractSymmetryOperation, InputScalarType<:Number}
    symops = [x for (x, y) in symops_and_amplitudes]
    amplitudes = [y for (x, y) in symops_and_amplitudes]
    return symmetry_reduce_parallel(hsr, symops, amplitudes; tol=tol)
end


function symmetry_reduce_parallel(
    hsr::HilbertSpaceRepresentation{QN, BR, DT},
    symops::AbstractArray{OperationType},
    amplitudes::AbstractArray{InputScalarType};
    tol::Real=Base.rtoldefault(float(real(InputScalarType)))
) where {QN, BR, DT, OperationType<:AbstractSymmetryOperation, InputScalarType<:Number}
    @debug "BEGIN symmetry_reduce_parallel"
    if length(symops) != length(amplitudes)
        throw(ArgumentError("number of symmetry operations and number of amplitudes should match ($(length(symops)) != $(length(amplitudes)))"))
    end
    if !all(isapprox(abs(y), one(abs(y)); atol=tol) || isapprox(abs(y), zero(abs(y)); atol=tol) for y in amplitudes)
        throw(ArgumentError("all amplitudes need to have norm 1"))
    end
    ScalarType = float(InputScalarType)

    n_basis = length(hsr.basis_list)
    @debug "Original Hilbert space dimension: $n_basis"

    # basis_mapping_index and basis_mapping_amplitude contain information about
    # which basis vector of the larger Hilbert space is included
    # in the basis of the smaller Hilbert space with what amplitude.
    basis_mapping_representative = Vector{Int}(undef, n_basis)
    fill!(basis_mapping_representative, -1)
    basis_mapping_amplitude = zeros(ScalarType, n_basis)

    subgroup_size = length(symops)

    size_estimate = n_basis ÷ max(1, subgroup_size - 1)
    @debug "Estimate for the reduced Hilbert space dimension: $size_estimate"

    inv_norm_cache = [inv(sqrt(float(i))) for i in 1:subgroup_size]

    nthreads = Threads.nthreads()
    local_reduced_basis_list = Vector{Vector{BR}}(undef, nthreads)
    for i in 1:nthreads
        local_reduced_basis_list[i] = BR[]
        sizehint!(local_reduced_basis_list[i], size_estimate ÷ nthreads + 1)
    end

    visited = zeros(UInt8, n_basis) # use UInt8 rather than Bool for thread safety

    local_basis_states = Matrix{BR}(undef, (nthreads, subgroup_size))
    local_basis_phases = ones(ScalarType, (nthreads, subgroup_size))
    local_basis_amplitudes = Vector{Dict{BR, ScalarType}}(undef, nthreads)
    for id in 1:nthreads
        local_basis_amplitudes[id] = Dict{BR, ScalarType}()
        sizehint!(local_basis_amplitudes[id], subgroup_size)
    end

    # Load balancing
    #   The representatives are the smaller binary numbers.
    #   Distribute them equally among threads.
    reorder = Int[]
    sizehint!(reorder, n_basis)
    nblocks = (n_basis + nthreads - 1) ÷ nthreads
    for i in 1:nthreads, j in 1:nblocks
        k = i + nthreads * (j-1)
        if 1 <= k <= n_basis
            push!(reorder, k)
        end
    end

    @debug "Starting reduction (parallel)"
    Threads.@threads for itemp in 1:n_basis
        ivec_p = reorder[itemp]
        (visited[ivec_p] != 0x0) && continue

        id = Threads.threadid()
        bvec = hsr.basis_list[ivec_p]

        # A basis binary representation is incompatible with the reduced Hilbert space if
        # (1) it is not the smallest among the its star, or
        # (2) its star is smaller than the representation
        compatible = true
        for i in 2:subgroup_size
            symop, ampl = symops[i], amplitudes[i]
            bvec_prime, sgn = symmetry_apply(hsr.hilbert_space, symop, bvec)
            if bvec_prime < bvec
                compatible = false
                break
            elseif bvec_prime == bvec && !isapprox(conj(ampl) * sgn, one(ScalarType); atol=tol)
                compatible = false
                break
            end
            local_basis_states[id, i] = bvec_prime
            local_basis_phases[id, i] = conj(ampl) * sgn
        end # for i
        (!compatible) && continue
        local_basis_states[id, 1] = bvec

        push!(local_reduced_basis_list[id], bvec)

        empty!(local_basis_amplitudes[id])
        for i in 1:subgroup_size
            bvec_prime = local_basis_states[id, i]
            local_basis_amplitudes[id][bvec_prime] = local_basis_phases[id, i] # Same bvec_prime, same p.
        end
        inv_norm = inv_norm_cache[length(local_basis_amplitudes[id])]        
        for (bvec_prime, amplitude) in local_basis_amplitudes[id]
            ivec_p_prime = hsr.basis_lookup[bvec_prime]
            visited[ivec_p_prime] = 0x1
            basis_mapping_representative[ivec_p_prime] = ivec_p
            basis_mapping_amplitude[ivec_p_prime] = amplitude * inv_norm
        end
    end
    @debug "Finished reduction (parallel)"

    reduced_basis_list = BR[]
    sizehint!(reduced_basis_list, sum(length(x) for x in local_reduced_basis_list))
    while !isempty(local_reduced_basis_list)
        lbl = pop!(local_reduced_basis_list)
        append!(reduced_basis_list, lbl)
    end
    @debug "Collected basis list"

    sort!(reduced_basis_list)
    @debug "Sorted basis list"

    # Basis vectors of the unreduced Hilbert space that are
    # not included in the reduced Hilbert space are marked as -1
    basis_mapping_index = Vector{Int}(undef, n_basis)
    fill!(basis_mapping_index, -1)

    Threads.@threads for ivec_r in eachindex(reduced_basis_list)
        bvec = reduced_basis_list[ivec_r]
        ivec_p = hsr.basis_lookup[bvec]
        basis_mapping_index[ivec_p] = ivec_r
    end
    @debug "Collected basis lookup (diagonal)"

    Threads.@threads for ivec_p_prime in eachindex(basis_mapping_representative)
        ivec_p = basis_mapping_representative[ivec_p_prime]
        (ivec_p <= 0) && continue  # not in this irrep
        (ivec_p_prime == ivec_p) && continue  # already in the lookup
        ivec_r = basis_mapping_index[ivec_p]
        basis_mapping_index[ivec_p_prime] = ivec_r
    end
    @debug "Collected basis lookup (offdiagonal)"

    @debug "END symmetry_reduce_parallel"
    HSR = HilbertSpaceRepresentation{QN, BR, DT}
    RHSR = ReducedHilbertSpaceRepresentation{HSR, BR, ScalarType}
    return RHSR(hsr, reduced_basis_list, basis_mapping_index, basis_mapping_amplitude)
end
