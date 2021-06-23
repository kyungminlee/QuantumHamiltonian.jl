export ProductHilbertSpace

struct ProductHilbertSpace{
    QN<:Tuple{Vararg{AbstractQuantumNumber}},
    N,
    S<:Tuple{Vararg{AbstractHilbertSpace{QN}, N}}
}<:AbstractHilbertSpace{QN}
    subspaces::S
    numsites::NTuple{N, Int}
    bitwidths::NTuple{N, Int}

    function ProductHilbertSpace(subspaces::S) where {QN, N, S<:Tuple{Vararg{AbstractHilbertSpace{QN}, N}}}
        ns = numsites.(subspaces)
        bitwidths = bitwidth.(subspaces)
        return new{QN, N, S}(subspaces, ns, bitwidths)
    end
    function ProductHilbertSpace(subspaces::Vararg{AbstractHilbertSpace})
        return ProductHilbertSpace(subspaces)
    end
end

qntype(::Type{<:ProductHilbertSpace{QN, N, S}}) where {QN, N, S} = QN
tagtype(::Type{<:ProductHilbertSpace{QN, N, S}}) where {QN, N, S} = NTuple{N, QN}

basespace(hs::ProductHilbertSpace) = hs

numsites(hs::ProductHilbertSpace) = sum(hs.numsites)
bitwidth(hs::ProductHilbertSpace) = sum(hs.bitwidths)

function get_site(hs::ProductHilbertSpace, isite::Integer)
    jsite = isite
    for (isub, sub) in enumerate(hs.subspaces)
        if isite <= hs.numsites[isub]
            return get_site(sub, isite)
        else
            isite -= hs.numsites[isub]
        end
    end
    throw(BoundsError(hs, jsite))
end

function bitwidth(hs::ProductHilbertSpace, isite::Integer)
    for (isub, sub) in enumerate(hs.subspaces)
        if isite <= hs.numsites[isub]
            return bitwidth(sub, isite)
        else
            isite -= hs.numsites[isub]
        end
    end
    return 0
end

function bitoffset(hs::ProductHilbertSpace, isite::Integer)
    offset = 0
    for (isub, sub) in enumerate(hs.subspaces)
        if isite <= hs.numsites[isub]
            return offset + bitoffset(sub, isite)
        else
            isite -= hs.numsites[isub]
        end
        offset += hs.bitwidths[isub]
    end
    return 0
end

function get_bitmask(hs::ProductHilbertSpace, isite::Integer, ::Type{BR}=UInt) where {BR<:Unsigned}
    offset = 0
    for (isub, sub) in enumerate(hs.subspaces)
        if isite <= hs.numsites[isub]
            return get_bitmask(sub, isite, BR) << offset
        else
            isite -= hs.numsites[isub]
        end
        offset += hs.bitwidths[isub]
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
    hs::ProductHilbertSpace,
    indexarray::CartesianIndex,
    ::Type{BR}=UInt
) where {BR<:Unsigned}
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


function get_state_index(hs::ProductHilbertSpace, binrep::BR, isite::Integer) where {BR<:Unsigned}
    return Int((binrep >> bitoffset(hs, isite)) & make_bitmask(bitwidth(hs, isite), BR)) + 1
end


function get_state(hs::ProductHilbertSpace, binrep::BR, isite::Integer) where {BR<:Unsigned}
    return get_site(hs, isite).states[get_state_index(hs, binrep, isite)]
end


#=
function Base.keys end
=#