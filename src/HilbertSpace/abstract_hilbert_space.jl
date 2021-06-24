export AbstractQuantumNumber
export AbstractHilbertSpace

export qntype, tagtype
export get_bitmask
export update, get_state_index, get_state

## TODO: Think about this
AbstractQuantumNumber = Integer

abstract type AbstractHilbertSpace{QN<:Tuple{Vararg{<:AbstractQuantumNumber}}} end

qntype(::T) where {T<:AbstractHilbertSpace} = qntype(T)
tagtype(::T) where {T<:AbstractHilbertSpace} = tagtype(T)


function get_bitmask(hs::AbstractHilbertSpace, ::Type{BR}=UInt)::BR where {BR<:Unsigned}
    return make_bitmask(bitwidth(hs), BR)
end

function get_bitmask(hs::AbstractHilbertSpace, isite::Integer, ::Type{BR}=UInt) where {BR<:Unsigned}
    off = bitoffset(hs, isite)
    bw = bitwidth(hs, isite)
    return make_bitmask(off+bw, off, BR)
end


"""
    update(hs, binrep, isite, new_state_index)

Update the binary representation of a basis state
by changing the state at site `isite` to a new local state specified by
`new_state_index`.
"""
@inline function update(
    hs::AbstractHilbertSpace,
    binrep::BR,
    isite::Integer,
    new_state_index::Integer
) where {BR<:Unsigned}
    @boundscheck if !(1 <= new_state_index <= dimension(get_site(hs, isite)))
        throw(BoundsError(1:dimension(get_site(hs, isite)), new_state_index))
    end
    mask = get_bitmask(hs, isite, BR)
    offset = bitoffset(hs, isite)
    return (binrep & (~mask)) | (BR(new_state_index-1) << offset)
end



"""
    get_state_index(hs, binrep, isite)

Get the *index of* the local state at site `isite` for the basis state
represented by `binrep`.
"""
function get_state_index(hs::AbstractHilbertSpace, binrep::BR, isite::Integer) where {BR<:Unsigned}
    return Int((binrep >> bitoffset(hs, isite)) & make_bitmask(bitwidth(hs, isite), BR)) + 1
end


"""
    get_state(hs, binrep, isite)

Get the local state at site `isite` for the basis state
represented by `binrep`. Returns an object of type `State`
"""
function get_state(hs::AbstractHilbertSpace, binrep::BR, isite::Integer) where {BR<:Unsigned}
    return get_site(hs, isite).states[get_state_index(hs, binrep, isite)]
end
