export symmetry_reduce_serial, symmetry_reduce_parallel
export make_symmetrizer

using LatticeTools


"""
    symmetry_reduce(hsr, symops, amplitudes, bvec; tol=√ϵ)

Returns bᵢ => ⟨B|ϕᵢ⟩, i.e., the basis state (represented by bᵢ, and the amount of that basis
state that overlaps with the input.
Returns the same amplitude as the `get_basis_index_amplitude` of the reduced Hilbert space representation

Basis states: |ϕᵢ⟩ with a representative bᵢ
input       : |B⟩
"""
function symmetry_reduce(
    hs::AbstractHilbertSpace,
    symops::AbstractArray{O},
    amplitudes::AbstractArray{SI},
    bvec::BR,
    ::Type{S}=float(SI);
    tol::Real=Base.rtoldefault(float(real(S)))
) where {BR<:Unsigned, O, SI<:Number, S<:Number}
    min_bvec = BR(bvec)
    min_amplitude = one(S)
    subgroup_size = length(symops)
    multiplicity = 1
    for i in 2:subgroup_size
        symop, ampl = symops[i], amplitudes[i]
        bvec_prime, ampl_prime = symmetry_apply(hs, symop, bvec, conj(ampl))
        if bvec_prime < min_bvec
            min_bvec = BR(bvec_prime)
            min_amplitude = S(ampl_prime)
        elseif bvec_prime == bvec
            if !isapprox(ampl_prime, one(S); atol=tol)
                return (representative=typemax(BR), amplitude=zero(S))
            else
                multiplicity += 1
            end
        end
    end
    @assert mod(subgroup_size, multiplicity) == 0
    return (representative=min_bvec, amplitude=conj(min_amplitude) / sqrt(subgroup_size / multiplicity))
end


function symmetry_reduce(
    hsr::AbstractHilbertSpaceRepresentation,
    symops::AbstractArray{OperationType},
    amplitudes::AbstractArray{InputScalarType},
    ::Type{BT}=DictIndexedVector;
    tol::Real=Base.rtoldefault(float(real(InputScalarType)))
) where {OperationType<:AbstractSymmetryOperation, InputScalarType<:Number, BT<:AbstractIndexedVector}
    symred = Threads.nthreads() == 1 ? symmetry_reduce_serial : symmetry_reduce_parallel
    return symred(hsr, symops, amplitudes, BT; tol=tol)
end


function symmetry_reduce_serial(
    hsr::AbstractHilbertSpaceRepresentation{BR, C},
    symops::AbstractArray{OperationType},
    amplitudes::AbstractArray{InputScalarType},
    ::Type{BT}=DictIndexedVector;
    tol::Real=Base.rtoldefault(float(real(InputScalarType)))
) where {BR, C, OperationType<:AbstractSymmetryOperation, InputScalarType<:Number, BT<:AbstractIndexedVector}
    if length(symops) != length(amplitudes)
        throw(ArgumentError("number of symmetry operations and number of amplitudes should match ($(length(symops)) != $(length(amplitudes)))"))
    elseif length(symops) < 1
        throw(ArgumentError("length of symops cannot be less than 1"))
    end
    if !all(isapprox(abs(y), one(abs(y)); atol=tol) || isapprox(abs(y), zero(abs(y)); atol=tol) for y in amplitudes)
        throw(ArgumentError("all amplitudes need to have norm 1 or 0 (monomial representation)"))
    elseif !isapprox(amplitudes[1], one(amplitudes[1]); atol=tol)
        throw(ArgumentError("amplitude of first element (identity) needs to be one"))
    end

    symops, amplitudes = let new_symops = OperationType[], new_amplitudes = InputScalarType[]
        for (x, y) in zip(symops, amplitudes)
            if !isapprox(abs(y), zero(abs(y)); atol=tol)
                push!(new_symops, x)
                push!(new_amplitudes, y)
            end
        end
        new_symops, new_amplitudes
    end

    ScalarType = float(promote_type(InputScalarType, C))

    n_basis = dimension(hsr)

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
        bvec = get_basis_state(hsr, ivec_p)
        compatible = true
        for i in 2:subgroup_size
            symop, ampl = symops[i], amplitudes[i]
            bvec_prime, ampl_prime = symmetry_apply(hsr, symop, bvec, conj(ampl))
            if bvec_prime < bvec
                compatible = false
                break
            elseif bvec_prime == bvec && !isapprox(ampl_prime, one(ScalarType); atol=tol)
                compatible = false
                break
            end
            basis_states[i] = bvec_prime
            basis_phases[i] = ampl_prime
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
            # TODO: check nested reduction
            ivec_p_prime, ampl_p_prime = get_basis_index_amplitude(hsr, bvec_prime)
            visited[ivec_p_prime] = true
            basis_mapping_representative[ivec_p_prime] = ivec_p
            basis_mapping_amplitude[ivec_p_prime] = amplitude * inv_norm * ampl_p_prime
        end
    end # for ivec_p

    basis_mapping_index = Vector{Int}(undef, n_basis)
    fill!(basis_mapping_index, -1)

    for (ivec_r, bvec) in enumerate(reduced_basis_list)
        ivec_p, _ = get_basis_index_amplitude(hsr, bvec)
        basis_mapping_index[ivec_p] = ivec_r
    end

    for (ivec_p_prime, ivec_p) in enumerate(basis_mapping_representative)
        (ivec_p <= 0) && continue  # not in this irrep
        (ivec_p_prime == ivec_p) && continue  # already in the lookup
        ivec_r = basis_mapping_index[ivec_p]
        basis_mapping_index[ivec_p_prime] = ivec_r
    end

    reduced_basis = BT(reduced_basis_list)
    return ReducedHilbertSpaceRepresentation(hsr, reduced_basis, basis_mapping_index, basis_mapping_amplitude)
end


function symmetry_reduce_parallel(
    hsr::AbstractHilbertSpaceRepresentation{BR, C},
    symops::AbstractArray{OperationType},
    amplitudes::AbstractArray{InputScalarType},
    ::Type{BT}=DictIndexedVector;
    tol::Real=Base.rtoldefault(float(real(InputScalarType)))
) where {BR, C, OperationType<:AbstractSymmetryOperation, InputScalarType<:Number, BT<:AbstractIndexedVector}
    @debug "BEGIN symmetry_reduce_parallel"
    if length(symops) != length(amplitudes)
        throw(ArgumentError("number of symmetry operations and number of amplitudes should match ($(length(symops)) != $(length(amplitudes)))"))
    elseif length(symops) < 1
        throw(ArgumentError("length of symops cannot be less than 1"))
    end
    if !all(isapprox(abs(y), one(abs(y)); atol=tol) || isapprox(abs(y), zero(abs(y)); atol=tol) for y in amplitudes)
        throw(ArgumentError("all amplitudes need to have norm 1 or 0 (monomial representation)"))
    elseif !isapprox(amplitudes[1], one(amplitudes[1]); atol=tol)
        throw(ArgumentError("amplitude of first element (identity) needs to be one"))
    end

    symops, amplitudes = let new_symops = OperationType[], new_amplitudes = InputScalarType[]
        for (x, y) in zip(symops, amplitudes)
            if !isapprox(abs(y), zero(abs(y)); atol=tol)
                push!(new_symops, x)
                push!(new_amplitudes, y)
            end
        end
        new_symops, new_amplitudes
    end

    ScalarType = float(promote_type(InputScalarType, C))

    n_basis = dimension(hsr)
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
        bvec = get_basis_state(hsr, ivec_p)

        # A basis binary representation is incompatible with the reduced Hilbert space if
        # (1) it is not the smallest among the its star, or
        # (2) its star is smaller than the representation
        compatible = true
        for i in 2:subgroup_size
            symop, ampl = symops[i], amplitudes[i]
            bvec_prime, ampl_prime = symmetry_apply(hsr, symop, bvec, conj(ampl))
            if bvec_prime < bvec
                compatible = false
                break
            elseif bvec_prime == bvec && !isapprox(ampl_prime, one(ScalarType); atol=tol)
                compatible = false
                break
            end
            local_basis_states[id, i] = bvec_prime
            local_basis_phases[id, i] = ampl_prime
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
            # TODO: check nested reduction
            ivec_p_prime, ampl_p_prime = get_basis_index_amplitude(hsr, bvec_prime)
            visited[ivec_p_prime] = 0x1
            basis_mapping_representative[ivec_p_prime] = ivec_p
            basis_mapping_amplitude[ivec_p_prime] = amplitude * inv_norm * ampl_p_prime
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
        ivec_p, _ = get_basis_index_amplitude(hsr, bvec)
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
    reduced_basis = BT(reduced_basis_list)
    return ReducedHilbertSpaceRepresentation(hsr, reduced_basis, basis_mapping_index, basis_mapping_amplitude)
end
