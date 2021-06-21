export ProductHilbertSpace

struct ProductHilbertSpace{
    QN<:Tuple{Vararg{AbstractQuantumNumber}},
    S<:Tuple{Vararg{AbstractHilbertSpace{QN}}}
} <: AbstractHilbertSpace{QN}
    subspaces::S

    function ProductHilbertSpace(subspaces::S) where {QN, S<:Tuple{Vararg{AbstractHilbertSpace{QN}}}}
        return new{QN, S}(subspaces)
    end
end

scalartype(::Type{<:ProductHilbertSpace}) = Bool
Base.valtype(::Type{<:ProductHilbertSpace}) = Bool
qntype(::Type{<:ProductHilbertSpace{QN, S}}) where {QN, S} = QN
basespace(hs::ProductHilbertSpace) = hs

numsites(hs::ProductHilbertSpace) = mapreduce(numsites, +, hs.subspaces)
bitwidth(hs::ProductHilbertSpace) = mapreduce(bitwidth, +, hs.subspaces)

function get_site(hs::ProductHilbertSpace, isite::Integer)
    jsite = isite
    for sub in hs.subspaces
        if isite <= numsites(sub)
            return get_site(sub, isite)
        else
            isite -= numsites(sub)
        end
    end
    throw(BoundsError(hs, jsite))
end

function bitwidth(hs::ProductHilbertSpace, isite::Integer)
    for sub in hs.subspaces
        if isite <= numsites(sub)
            return bitwidth(sub, isite)
        else
            isite -= numsites(sub)
        end
    end
    return 0
end

function bitoffset(hs::ProductHilbertSpace, isite::Integer)
    offset = 0
    for sub in hs.subspaces
        if isite <= numsites(sub)
            return offset + bitoffset(sub, isite)
        else
            isite -= numsites(sub)
        end
        offset += bitwidth(sub)
    end
    return 0
end

function get_bitmask(hs::ProductHilbertSpace, isite::Integer, ::Type{BR}=UInt) where {BR<:Unsigned}
    offset = 0
    for sub in hs.subspaces
        if isite <= numsites(sub)
            return get_bitmask(sub, isite, BR) << offset
        else
            isite -= numsites(sub)
        end
        offset += bitwidth(sub)
    end
    return 0
end

function get_bitmask(hs::ProductHilbertSpace, ::Type{BR}=UInt)::BR where {BR<:Unsigned}
    return make_bitmask(bitwidth(hs), BR)
end

function get_quantum_numbers(hs::ProductHilbertSpace)
    combine(x, y) = unique(sort([a .+ b for b in y for a in x]))
    mapreduce(get_quantum_numbers, combine, hs.subspaces)
end

function get_quantum_number(hs::ProductHilbertSpace, binrep::Unsigned)
    bvecs = bitsplit(bitwidth.(hs.subspaces), binrep)
    qns = get_quantum_number.(hs.subspaces, bvecs)
    return mapreduce(identity, (x,y) -> x .+ y, qns)
end

function get_tags(hs::ProductHilbertSpace)
    qns = vcat(collect(Iterators.product(get_tags.(hs.subspaces)...))...)
    return sort(qns, by=reverse)
end

function get_tag(hs::ProductHilbertSpace, binrep::Unsigned)
    bvecs = bitsplit(bitwidth.(hs.subspaces), binrep)
    return get_tag.(hs.subspaces, bvecs)
end

function extract(hs::ProductHilbertSpace, binrep::Unsigned)
    bvecs = bitsplit(bitwidth.(hs.subspaces), binrep)
    indices = ((x,y) -> collect(extract(x,y).I)).(hs.subspaces, bvecs)
    return CartesianIndex(vcat(indices...)...)
end

function arraysplit(sizes::AbstractVector{<:Integer}, array::AbstractVector{T}) where {T}
    if sum(sizes) != length(array)
        throw(ArgumentError("size mismatch"))
    end
    out = Vector{T}[]
    i = 1
    for s in sizes
        push!(out, array[i:i+s-1])
        i += s
    end
    return out
end


function compress(
    hs::ProductHilbertSpace{QN},
    indexarray::CartesianIndex,
    ::Type{BR}=UInt
) where {QN, BR<:Unsigned}
    ns = collect(map(numsites, hs.subspaces))
    splitindices = map(x -> CartesianIndex(x...), arraysplit(ns, collect(indexarray.I)))
    bvecs = tuple([compress(x, y, BR) for (x,y) in zip(hs.subspaces, splitindices)]...)
    return bitjoin(bitwidth.(hs.subspaces), bvecs)
end


@inline function update(
    hs::ProductHilbertSpace,
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


@inline function get_state_index(hs::ProductHilbertSpace, binrep::BR, isite::Integer) where {BR<:Unsigned}
    return Int((binrep >> bitoffset(hs, isite)) & make_bitmask(bitwidth(hs, isite), BR)) + 1
end


@inline function get_state(hs::ProductHilbertSpace, binrep::BR, isite::Integer) where {BR<:Unsigned}
    return get_site(hs, isite).states[get_state_index(hs, binrep, isite)]
end


#=
function Base.keys end
=#