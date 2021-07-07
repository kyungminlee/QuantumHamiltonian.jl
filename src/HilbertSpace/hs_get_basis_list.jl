
function hs_get_basis_list(
    hs::AbstractHilbertSpace,
    ::Type{BR}=UInt
)::Vector{BR} where {BR<:Unsigned}
    if sizeof(BR) * 8 <= bitwidth(hs)
        throw(ArgumentError("type $(BR) not enough to represent the hilbert space (need $(bitwidth(hs)) bits)"))
    end
    basis_list = BR[]
    for indexarray in keys(hs)
        push!(basis_list, compress(hs, indexarray, BR))
    end
    sort!(basis_list)
    return basis_list
end


"""
    hs_get_basis_list(hs, allowed_quantum_numbers, binary_type=UInt)

Generate a basis for the `HilbertSpaceSector`.
"""
function hs_get_basis_list(
    hs::AbstractHilbertSpace{QN},
    allowed_quantum_numbers::AbstractVector{QN},
    ::Type{BR}=UInt
)::Vector{BR} where {QN, BR<:Unsigned}
    if sizeof(BR) * 8 <= bitwidth(hs)
        throw(ArgumentError("type $(BR) not enough to represent the hilbert space (need $(bitwidth(hs)) bits)"))
    end
    if isempty(intersect(Set(allowed_quantum_numbers), quantum_number_sectors(hs)))
        return BR[]
    end

    n_sites = numsites(hs)

    quantum_numbers = Vector{QN}[
        let site = get_site(hs, isite)
            [
                # state.quantum_number
                #     for state in site.states
                get_quantum_number(site, istate)
                    for istate in 1:dimension(site)
            ]
        end
            for isite in 1:n_sites
    ]

    # `qn_possible[i]` contains all possible quantum numbers in the subspace
    # spanned by sites 1 ⊗ 2 ⊗ ... ⊗ (i-1).
    # `qn_requested[i]` contains quantum numbers that are requested in the subspace
    # spanned by sites 1 ⊗ 2 ⊗ ... ⊗ (i-1).
    # `qn_schedule` is the intersection of these two, and are the ones we need to
    # generate the basis.
    qn_possible::Vector{Vector{QN}} = let qn_possible = Vector{Vector{QN}}(undef, n_sites+1)
        qn_possible[1] = [tuplezero(QN)]
        for i in 1:n_sites
            a = Set{QN}()
            for qi in quantum_numbers[i], qa in qn_possible[i]
                push!(a, qa .+ qi)
            end
            qn_possible[i+1] = sort(collect(a))
        end
        qn_possible
    end

    qn_requested::Vector{Vector{QN}} = let qn_requested = Vector{Vector{QN}}(undef, n_sites+1)
        qn_requested[n_sites+1] = sort(allowed_quantum_numbers)
        for i in n_sites:-1:1
            a = Set{QN}()
            for qi in quantum_numbers[i], qa in qn_requested[i+1]
                push!(a, qa .- qi)
            end
            qn_requested[i] = sort(collect(a))
        end
        qn_requested
    end

    qn_schedule = [intersect(x,y) for (x, y) in zip(qn_requested, qn_possible)]

    sector_basis_list = Dict{QN, Vector{BR}}(tuplezero(QN) => BR[BR(0x0)])
    new_sector_basis_list = Dict{QN, Vector{BR}}()

    #sl = Threads.SpinLock()
    for i in 1:n_sites
        empty!(new_sector_basis_list)
        # @threads
        for q in qn_schedule[i+1]
            new_sector_basis_list_q = BR[]
            for (i_state, q_curr) in enumerate(quantum_numbers[i])
                q_prev::QN = q .- q_curr
                if haskey(sector_basis_list, q_prev)
                    append!(
                        new_sector_basis_list_q,
                        # (s | (represent(hs.sites[i], i_state, BR) << hs.bitoffsets[i])) for s in sector_basis_list[q_prev]
                        (s | (compress(get_site(hs, i), i_state, BR) << bitoffset(hs, i))) for s in sector_basis_list[q_prev]
                    )
                end
            end
            #lock(sl)
            new_sector_basis_list[q] = new_sector_basis_list_q
            #unlock(sl)
        end
        sector_basis_list, new_sector_basis_list = new_sector_basis_list, sector_basis_list
    end

    basis_list = BR[]
    qs = collect(keys(sector_basis_list))
    for q in qs
        states = pop!(sector_basis_list, q)
        sort!(states)
        basis_list = merge_vec(basis_list, states)
        GC.gc()
    end
    @assert issorted(basis_list)
    return basis_list
end
