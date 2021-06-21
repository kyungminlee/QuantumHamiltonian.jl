export HilbertSpaceSector
export scalartype
export qntype
export basespace
export bitwidth


"""
    HilbertSpaceSector{QN}

Hilbert space sector.
"""
struct HilbertSpaceSector{HS<:AbstractHilbertSpace, QN<:Tuple{Vararg{<:AbstractQuantumNumber}}}<:AbstractHilbertSpace{QN}
    parent::HS
    allowed_quantum_numbers::Set{QN}

    """
        HilbertSpaceSector(parent::HilbertSpace{QN}) where QN
    """
    function HilbertSpaceSector(parent::HS) where {HS<:AbstractHilbertSpace}
        QN = qntype(HS)
        sectors = quantum_number_sectors(parent)
        return new{HS, QN}(parent, Set(sectors))
    end

    """
        HilbertSpaceSector(parent::HilbertSpace{QN}, allowed::Integer) where {QN<:Tuple{<:Integer}}
    """
    function HilbertSpaceSector(parent::HS, allowed::Integer) where {HS<:AbstractHilbertSpace{<:Tuple{<:Integer}}}
        QN = qntype(HS)
        sectors = Set{QN}(quantum_number_sectors(parent))
        return new{HS, QN}(parent, intersect(sectors, Set([(allowed,)])))
    end

    """
        HilbertSpaceSector(parent::HilbertSpace{QN}, allowed::QN) where QN
    """
    function HilbertSpaceSector(parent::AbstractHilbertSpace{QN}, allowed::QN) where {QN}
        HS = typeof(parent)
        sectors = Set{QN}(quantum_number_sectors(parent))
        return new{HS, QN}(parent, intersect(sectors, Set([allowed])))
    end

    """
        HilbertSpaceSector(parent::HilbertSpace{QN}, allowed::Union{AbstractSet{<:Integer}, AbstractVector{<:Integer}}) where QN
    """
    function HilbertSpaceSector(
        parent::AbstractHilbertSpace{QN},
        allowed::Union{AbstractSet{<:Integer}, AbstractVector{<:Integer}}
    ) where {QN}
        HS = typeof(parent)
        sectors = Set{QN}(quantum_number_sectors(parent))
        return new{HS, QN}(parent, intersect(sectors, Set((x,) for x in allowed)))
    end

    function HilbertSpaceSector(
        parent::AbstractHilbertSpace{QN},
        allowed::Union{AbstractSet{QN}, AbstractVector{QN}}
    ) where QN
        HS = typeof(parent)
        sectors = Set{QN}(quantum_number_sectors(parent))
        return new{HS, QN}(parent, intersect(sectors, Set(allowed)))
    end
end


"""
    scalartype(arg::Type{HilbertSpaceSector{HS, QN}})

Returns the scalar type of the given hilbert space sector type.
For HilbertSpaceSector{QN}, it is always `Bool`.
"""
scalartype(::Type{HilbertSpaceSector{HS, QN}}) where {HS, QN} = Bool
scalartype(::HilbertSpaceSector{HS, QN}) where {HS, QN} = Bool


"""
    valtype(arg::Type{HilbertSpaceSector{HS, QN}})

Returns the `valtype` (scalar type) of the given hilbert space sector type.
For HilbertSpaceSector{QN}, it is always `Bool`.
"""
Base.valtype(::Type{HilbertSpaceSector{HS, QN}}) where {HS, QN} = Bool
Base.valtype(::HilbertSpaceSector{HS, QN}) where {HS, QN} = Bool


"""
    qntype(arg::Type{HilbertSpaceSector{QN}})

Returns the quantum number type of the given hilbert space sector type.
"""
qntype(::Type{HilbertSpaceSector{HS, QN}}) where {HS, QN} = QN
qntype(::HilbertSpaceSector{HS, QN}) where {HS, QN} = QN


"""
    basespace(hss)

Get the base space of the `HilbertSpaceSector`, which is
its parent `HilbertSpace` (with no quantum number restriction).
"""
basespace(hss::HilbertSpaceSector{HS, QN}) where {HS, QN} = basespace(hss.parent)::HS

#bitwidth(hss::HilbertSpaceSector) = bitwidth(basespace(hss))

function Base.:(==)(lhs::HilbertSpaceSector{HS, Q1}, rhs::HilbertSpaceSector{HS, Q2}) where {HS, Q1, Q2}
    return (
        basespace(lhs) == basespace(rhs) &&
        lhs.allowed_quantum_numbers == rhs.allowed_quantum_numbers
    )
end




"""
    hs_get_basis_list(hss, binary_type=UInt)

Generate a basis for the `HilbertSpaceSector`.
"""
function hs_get_basis_list(hss::HilbertSpaceSector{HS, QN}, binary_type::Type{BR}=UInt)::Vector{BR} where {HS, QN, BR<:Unsigned}
    hs = hss.parent
    if sizeof(BR) * 8 <= bitwidth(hs)
        throw(ArgumentError("type $(BR) not enough to represent the hilbert space (need $(bitwidth(hs)) bits)"))
    end
    if isempty(intersect(hss.allowed_quantum_numbers, quantum_number_sectors(hs)))
        return BR[]
    end

    quantum_numbers = Vector{QN}[QN[state.quantum_number for state in site.states] for site in hs.sites]

    n_sites = length(hs.sites)

    # `qn_possible[i]` contains all possible quantum numbers in the subspace
    # spanned by sites 1 ⊗ 2 ⊗ ... ⊗ (i-1).
    # `qn_requested[i]` contains quantum numbers that are requested in the subspace
    # spanned by sites 1 ⊗ 2 ⊗ ... ⊗ (i-1).
    # `qn_schedule` is the intersection of these two, and are the ones we need to
    # generate the basis.
    qn_possible::Vector{Vector{QN}} = let qn_possible = Vector{Vector{QN}}(undef, n_sites+1)
        qn_possible[1] = [tuplezero(QN)]
        for i in eachindex(hs.sites)
            a = Set{QN}()
            for qi in quantum_numbers[i], qa in qn_possible[i]
                push!(a, qa .+ qi)
            end
            qn_possible[i+1] = sort(collect(a))
        end
        qn_possible
    end

    qn_requested::Vector{Vector{QN}} = let qn_requested = Vector{Vector{QN}}(undef, n_sites+1)
        qn_requested[n_sites+1] = sort(collect(hss.allowed_quantum_numbers))
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
                        (s | (represent(hs.sites[i], i_state, BR) << hs.bitoffsets[i])) for s in sector_basis_list[q_prev]
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
        basis_list = merge_vec(basis_list, states)
        GC.gc()
    end
    @assert issorted(basis_list)
    return basis_list
end




for fname in [
    :bitwidth,
    :get_bitmask,
    :quantum_number_sectors,
    :get_quantum_number,
    :extract,
    :compress,
    :update,
    :get_state_index,
    :get_state,
]
    @eval begin
        """
            $($fname)(hss::HilbertSpaceSector, args...;kwargs...)

        Call `$($fname)` with basespace of `hss`.
        """
        @inline $fname(hss::HilbertSpaceSector, args...;kwargs...) = $fname(hss.parent, args...; kwargs...)
    end
end

