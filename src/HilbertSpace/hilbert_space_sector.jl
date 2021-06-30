export HilbertSpaceSector
export scalartype
export qntype
export basespace
export bitwidth


"""
    HilbertSpaceSector{QN}

Hilbert space sector.
"""
struct HilbertSpaceSector{
    QN<:Tuple{Vararg{<:AbstractQuantumNumber}},
    HS<:AbstractHilbertSpace{QN}
}<:AbstractHilbertSpace{QN}

    parent::HS
    allowed_quantum_numbers::Set{QN}

    """
        HilbertSpaceSector(parent::HilbertSpace)
    """
    function HilbertSpaceSector(parent::HS) where {HS<:AbstractHilbertSpace}
        QN = qntype(HS)
        sectors = quantum_number_sectors(parent)
        return new{QN, HS}(parent, Set(sectors))
    end

    """
        HilbertSpaceSector(parent::HilbertSpace{QN, TT}, allowed::Integer) where {QN<:Tuple{<:Integer}}
    """
    function HilbertSpaceSector(parent::HS, allowed::Integer) where {HS<:AbstractHilbertSpace{<:Tuple{<:Integer}}}
        QN = qntype(HS)
        sectors = Set{QN}(quantum_number_sectors(parent))
        return new{QN, HS}(parent, intersect(sectors, Set([(allowed,)])))
    end

    """
        HilbertSpaceSector(parent::HilbertSpace{QN, TT}, allowed::QN)
    """
    function HilbertSpaceSector(parent::AbstractHilbertSpace{QN}, allowed::QN) where {QN}
        HS = typeof(parent)
        sectors = Set{QN}(quantum_number_sectors(parent))
        return new{QN, HS}(parent, intersect(sectors, Set([allowed])))
    end

    """
        HilbertSpaceSector(parent::HilbertSpace{QN, TT}, allowed::Union{AbstractSet{<:Integer}, AbstractVector{<:Integer}})
    """
    function HilbertSpaceSector(
        parent::AbstractHilbertSpace{QN},
        allowed::Union{AbstractSet{<:Integer}, AbstractVector{<:Integer}}
    ) where {QN}
        HS = typeof(parent)
        sectors = Set{QN}(quantum_number_sectors(parent))
        return new{QN, HS}(parent, intersect(sectors, Set((x,) for x in allowed)))
    end

    function HilbertSpaceSector(
        parent::AbstractHilbertSpace{QN},
        allowed::Union{AbstractSet{QN}, AbstractVector{QN}}
    ) where QN
        HS = typeof(parent)
        sectors = Set{QN}(quantum_number_sectors(parent))
        return new{QN, HS}(parent, intersect(sectors, Set(allowed)))
    end
end


"""
    qntype(arg::Type{HilbertSpaceSector{QN}})

Returns the quantum number type of the given hilbert space sector type.
"""
qntype(::Type{<:HilbertSpaceSector{QN, HS}}) where {QN, HS} = QN
tagtype(::Type{<:HilbertSpaceSector{QN, HS}}) where {QN, HS} = QN

function Base.:(==)(lhs::HilbertSpaceSector{QN, HS}, rhs::HilbertSpaceSector{QN, HS}) where {QN, HS}
    return (
        basespace(lhs) == basespace(rhs) &&
        lhs.allowed_quantum_numbers == rhs.allowed_quantum_numbers
    )
end


"""
    hs_get_basis_list(hss, binary_type=UInt)

Generate a basis for the `HilbertSpaceSector`.
"""
function hs_get_basis_list(hss::HilbertSpaceSector{QN, HS}, ::Type{BR}=UInt)::Vector{BR} where {QN, HS, BR<:Unsigned}
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


get_quantum_numbers(hss::HilbertSpaceSector) = sort(collect(hss.allowed_quantum_numbers))
get_tags(hss::HilbertSpaceSector) = get_quantum_numbers(hss)
quantum_number_sectors(hss::HilbertSpaceSector) = get_quantum_numbers(hss)


for fname in [
    :basespace,
    :numsites,
    :get_site,
    :bitwidth,
    :bitoffset,
    :get_quantum_number,
    :get_tag,
    :extract,
    :compress,
]
    @eval begin
        """
            $($fname)(hss::HilbertSpaceSector, args...;kwargs...)

        Call `$($fname)` with basespace of `hss`.
        """
        @inline $fname(hss::HilbertSpaceSector, args...;kwargs...) = $fname(hss.parent, args...; kwargs...)
    end
end

