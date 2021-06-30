export HilbertSpace
import LatticeTools.numsites
export numsites, get_site
export quantum_number_sectors, get_quantum_numbers, get_quantum_number
export extract, uncompress, compress, update, get_state, get_state_index
export get_bitmask
export get_tag, get_tags
export bitwidth, bitoffset
export scalartype, qntype
export basespace


"""
    HilbertSpace{QN, TT}

Abstract Hilbert space with quantum number type `QN`.

# Examples
```jldoctest
julia> using QuantumHamiltonian

julia> spin_site = Site([State("Up", +1), State("Dn", -1)]);

julia> hs = HilbertSpace([spin_site, spin_site])
HilbertSpace{Tuple{Int64}, Tuple{Int64}}(Site{Tuple{Int64}}[Site{Tuple{Int64}}(State{Tuple{Int64}}[State{Tuple{Int64}}("Up", (1,)), State{Tuple{Int64}}("Dn", (-1,))]), Site{Tuple{Int64}}(State{Tuple{Int64}}[State{Tuple{Int64}}("Up", (1,)), State{Tuple{Int64}}("Dn", (-1,))])], [1, 1], [0, 1, 2])
```
"""
struct HilbertSpace{QN<:Tuple{Vararg{<:AbstractQuantumNumber}}, TagType}<:AbstractHilbertSpace{QN}
    sites::Vector{Site{QN}}
    bitwidths::Vector{Int}
    bitoffsets::Vector{Int}

    function HilbertSpace(sites::AbstractArray{Site{QN}, 1}) where QN
        bitwidths = map(bitwidth, sites)
        bitoffsets = Int[0, cumsum(bitwidths)...]
        new{QN, QN}(sites, bitwidths, bitoffsets)
    end

    function HilbertSpace{QN, TT}(sites::AbstractArray{Site{QN}, 1}) where {QN, TT}
        bitwidths = map(bitwidth, sites)
        bitoffsets = Int[0, cumsum(bitwidths)...]
        new{QN, TT}(sites, bitwidths, bitoffsets)
    end
end


"""
    qntype(::Type{HilbertSpace{QN, TT}})

Returns the quantum number type of the given hilbert space type.
"""
qntype(::Type{HilbertSpace{QN, TT}}) where {QN, TT} = QN


"""
    tagtype(::Type{HilbertSpace{QN, TT}})

Returns the tag type of the given Hilbert space type, which is its qntype.
"""
tagtype(::Type{HilbertSpace{QN, TT}}) where {QN, TT} = TT


"""
    basespace(hs)

Get the base space of the HilbertSpace `hs`, which is itself.
"""
basespace(hs::HilbertSpace) = hs

numsites(hs::HilbertSpace) = length(hs.sites)

get_site(hs::HilbertSpace, i::Integer) = hs.sites[i]


"""
    bitwidth(hs, [isite])
Total number of bits

```jldoctest
julia> using QuantumHamiltonian

julia> spin_site = Site([State("Up", +1), State("Dn", -1)]);

julia> hs = HilbertSpace([spin_site, spin_site, spin_site,]);

julia> bitwidth(hs)
3
```
"""
bitwidth(hs::HilbertSpace) = hs.bitoffsets[end]

bitwidth(hs::HilbertSpace, isite::Integer) = hs.bitwidths[isite]

bitoffset(hs::HilbertSpace, isite::Integer) = hs.bitoffsets[isite]

Base.:(==)(lhs::HilbertSpace, rhs::HilbertSpace) = lhs.sites == rhs.sites


"""
    get_bitmask(hs, [isite, [type]])

Get the bitmask of the site `isite`, optionally in type `type`
"""
# optimization from generic using bitoffsets
function get_bitmask(hs::HilbertSpace, isite::Integer, ::Type{BR}=UInt)::BR where {BR<:Unsigned}
    return make_bitmask(hs.bitoffsets[isite+1], hs.bitoffsets[isite], BR)
end

# function get_bitmask(hs::HilbertSpace, ::Type{BR}=UInt)::BR where {BR<:Unsigned}
#     return make_bitmask(bitwidth(hs), BR)
# end


"""
    get_quantum_numbers(hs)

Return a sorted list of quantum numbers of the hilbert space `hs`.
"""
function get_quantum_numbers(hs::HilbertSpace{QN, TT})::Vector{QN} where {QN, TT}
    qns = Set{QN}([tuplezero(QN)])
    for site in hs.sites
        qns_next = Set{QN}()
        for state in site.states, q in qns
            push!(qns_next, q .+ state.quantum_number)
        end
        qns = qns_next
    end
    return sort(collect(qns))
end

quantum_number_sectors(hs::HilbertSpace) = get_quantum_numbers(hs)

"""
    get_quantum_number(hs, rep)

Get the quantum number of `rep`, which is either a binary representation, or a `CartesianIndex`.
"""
function get_quantum_number(hs::HilbertSpace, binrep::Unsigned)
    return mapreduce(
        identity,
        tupleadd,
        let i = get_state_index(hs, binrep, isite)
            site.states[i].quantum_number
        end
            for (isite, site) in enumerate(hs.sites)
    )
end

function get_quantum_number(hs::HilbertSpace, indexarray::AbstractVector{<:Integer})
    return mapreduce(
        identity,
        tupleadd,
        site.states[indexarray[isite]].quantum_number
            for (isite, site) in enumerate(hs.sites)
    )
end

# TODO: introduce TagStrategy template argument for HilbertSpace
get_tag(hs::HilbertSpace, binrep::Unsigned) = get_quantum_number(hs, binrep)
get_tags(hs::HilbertSpace) = quantum_number_sectors(hs)

"""
    extract(hs, binrep)

Convert binary representation to an array of indices (of states)

# Examples
```jldoctest
julia> using QuantumHamiltonian

julia> spin_site = Site([State("Up", +1), State("Dn", -1)]);

julia> hs = HilbertSpace([spin_site, spin_site]);

julia> extract(hs, 0x03)
CartesianIndex(2, 2)
```
"""
function extract(hs::HilbertSpace, binrep::BR)::CartesianIndex where {BR<:Unsigned}
    out = Int[]
    for (isite, site) in enumerate(hs.sites)
        @inbounds mask = make_bitmask(hs.bitwidths[isite], BR)
        index = Int(binrep & mask) + 1
        @boundscheck if !(1 <= index <= length(site.states))
            throw(BoundsError(1:length(site.states), index))
        end
        push!(out, index)
        binrep = binrep >> hs.bitwidths[isite]
    end
    return CartesianIndex(out...)
end


"""
    uncompress(hs, binrep)

Convert binary representation to an array of indices (of states)

# Examples
```jldoctest
julia> using QuantumHamiltonian

julia> spin_site = Site([State("Up", +1), State("Dn", -1)]);

julia> hs = HilbertSpace([spin_site, spin_site]);

julia> extract(hs, 0x03)
CartesianIndex(2, 2)
```
"""
function uncompress(hs::HilbertSpace, binrep::BR)::CartesianIndex where {BR<:Unsigned}
    out = Int[]
    for (isite, site) in enumerate(hs.sites)
        @inbounds mask = make_bitmask(hs.bitwidths[isite], BR)
        index = Int(binrep & mask) + 1
        @boundscheck if !(1 <= index <= length(site.states))
            throw(BoundsError(1:length(site.states), index))
        end
        push!(out, index)
        binrep = binrep >> hs.bitwidths[isite]
    end
    return CartesianIndex(out...)
end


"""
    compress(hs, indexarray::CartesianIndex, binary_type=UInt)

Convert a cartesian index (a of state) to its binary representation

# Examples
```jldoctest
julia> using QuantumHamiltonian

julia> spin_site = Site([State("Up", +1), State("Dn", -1)]);

julia> hs = HilbertSpace([spin_site, spin_site]);

julia> compress(hs, CartesianIndex(2,2))
0x0000000000000003
```
"""
function compress(
    hs::HilbertSpace,
    indexarray::CartesianIndex,
    ::Type{BR}=UInt
) where {BR<:Unsigned}
    if length(indexarray) != length(hs.sites)
        throw(ArgumentError("length of indexarray should be the number of sites"))
    end
    binrep = zero(BR)
    for (isite, site) in enumerate(hs.sites)
        @boundscheck if !(1 <= indexarray[isite] <= dimension(site))
            throw(BoundsError(1:dimension(site), indexarray[isite]))
        end
        binrep |= (BR(indexarray[isite] - 1) << hs.bitoffsets[isite] )
    end
    return binrep
end


function Base.keys(hs::HilbertSpace)
    return CartesianIndices(((1:length(site.states) for site in hs.sites)...,))
end


function hs_get_basis_list(
    hs::HilbertSpace,
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
    hs::HilbertSpace{QN, TT},
    allowed_quantum_numbers::AbstractVector{QN},
    ::Type{BR}=UInt
)::Vector{BR} where {QN, TT, BR<:Unsigned}
    if sizeof(BR) * 8 <= bitwidth(hs)
        throw(ArgumentError("type $(BR) not enough to represent the hilbert space (need $(bitwidth(hs)) bits)"))
    end
    if isempty(intersect(Set(allowed_quantum_numbers), quantum_number_sectors(hs)))
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
