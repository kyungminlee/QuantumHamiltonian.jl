export HilbertSpace
import LatticeTools.numsites
export numsites, get_site
export quantum_number_sectors, get_quantum_numbers, get_quantum_number
export extract, uncompress, compress, update, get_state, get_state_index
export get_bitmask
export get_tag, get_tags
export bitwidth, bitoffset
export scalartype, qntype, tagtype
export basespace


"""
    HilbertSpace{QN, TS}

Abstract Hilbert space with quantum number type `QN`.

# Examples
```jldoctest
julia> using QuantumHamiltonian

julia> spin_site = Site([State("Up", +1), State("Dn", -1)]);

julia> hs = HilbertSpace([spin_site, spin_site]);
```
"""
struct HilbertSpace{QN<:Tuple{Vararg{<:AbstractQuantumNumber}}}<:AbstractHilbertSpace{QN}
    sites::Vector{Site{QN}}
    bitwidths::Vector{Int}
    bitoffsets::Vector{Int}

    function HilbertSpace(sites::AbstractVector{Site{QN}}) where {QN}
        bitwidths = map(bitwidth, sites)
        bitoffsets = Int[0, cumsum(bitwidths)...]
        new{QN}(sites, bitwidths, bitoffsets)
    end
end


# """
#     qntype(::Type{HilbertSpace{QN}})

# Returns the quantum number type of the given hilbert space type.
# """
# qntype(::Type{HilbertSpace{QN}}) where {QN, TS} = QN


# """
#     tagtype(::Type{HilbertSpace{QN}})

# Returns the tag type of the given Hilbert space type, which is its qntype.
# """
tagtype(::Type{<:HilbertSpace{QN}}, ::Val{:QuantumNumberAsTag}) where {QN} = QN


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
    get_quantum_numbers(hs)

Return a sorted list of quantum numbers of the hilbert space `hs`.
"""
function get_quantum_numbers(hs::HilbertSpace{QN})::Vector{QN} where {QN}
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


get_tag(hs::HilbertSpace{QN}, binrep::Unsigned, ::Val{:QuantumNumberAsTag}) where {QN} = get_quantum_number(hs, binrep)
get_tags(hs::HilbertSpace{QN}, ::Val{:QuantumNumberAsTag}) where {QN} = quantum_number_sectors(hs)


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


# == performance specializations ==


"""
    get_state_index(hs, binrep, isite)

Get the *index of* the local state at site `isite` for the basis state
represented by `binrep`.
"""
# performance specialization
function get_state_index(hs::HilbertSpace, binrep::BR, isite::Integer) where {BR<:Unsigned}
    sitebinrep = (binrep >> bitoffset(hs, isite)) & make_bitmask(bitwidth(hs, isite), BR)
    return Int(sitebinrep) + 1
end

"""
    get_bitmask(hs, [isite, [type]])

Get the bitmask of the site `isite`, optionally in type `type`
"""
# performance specialization
function get_bitmask(hs::HilbertSpace, isite::Integer, ::Type{BR}=UInt)::BR where {BR<:Unsigned}
    return make_bitmask(hs.bitoffsets[isite+1], hs.bitoffsets[isite], BR)
end

# performance specialization
function Base.keys(hs::HilbertSpace)
    return CartesianIndices(((length(site.states) for site in hs.sites)...,))
end

