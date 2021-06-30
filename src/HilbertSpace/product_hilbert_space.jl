export ProductHilbertSpace

struct ProductHilbertSpace{
    QN<:Tuple{Vararg{AbstractQuantumNumber}},
    N,
    S<:Tuple{Vararg{AbstractHilbertSpace{QN}, N}},
}<:AbstractHilbertSpace{QN}
    subspaces::S
    numsites::NTuple{N, Int}
    bitwidths::NTuple{N, Int}

    function ProductHilbertSpace(subspaces::S) where {QN, N, S<:Tuple{Vararg{AbstractHilbertSpace{QN}, N}}}
        ns = numsites.(subspaces)
        bitwidths = bitwidth.(subspaces)
        return new{QN, N, S}(subspaces, ns, bitwidths)
    end
    function ProductHilbertSpace(subspaces::AbstractHilbertSpace...)
        return ProductHilbertSpace(subspaces)
    end
end

# qntype(::Type{<:ProductHilbertSpace{QN, N, TS, S, TT}}) where {QN, N, TS, S, TT} = QN
# tagtype(::Type{<:ProductHilbertSpace{QN, N, TS, S, TT}}) where {QN, N, S, TT} = NTuple{N, QN}

tagtype(::Type{<:ProductHilbertSpace{QN, N, S}}, x::Val) where {QN, N, S}= Tuple{tagtype.(S.parameters, x)...}


basespace(hs::ProductHilbertSpace) = hs

numsites(hs::ProductHilbertSpace) = sum(hs.numsites)
bitwidth(hs::ProductHilbertSpace) = sum(hs.bitwidths)

function get_site(hs::ProductHilbertSpace, isite::Integer)
    isite > 0 || throw(BoundsError(hs, isite))
    jsite = isite
    for (isub, sub) in enumerate(hs.subspaces)
        if jsite <= hs.numsites[isub]
            return get_site(sub, jsite)
        else
            jsite -= hs.numsites[isub]
        end
    end
    throw(BoundsError(hs, isite))
end

function bitwidth(hs::ProductHilbertSpace, isite::Integer)
    isite > 0 || throw(BoundsError(hs, isite))
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
    isite > 0 || throw(BoundsError(hs, isite))
    offset = 0
    for (isub, sub) in enumerate(hs.subspaces)
        if isite <= hs.numsites[isub]
            return offset + bitoffset(sub, isite)
        else
            isite -= hs.numsites[isub]
        end
        offset += hs.bitwidths[isub]
    end
    return offset
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

function get_tags(hs::ProductHilbertSpace{QN, N, S}, ::Val{:QuantumNumberAsTag}) where {QN, N, S}
    qns = vcat(collect(Iterators.product(get_tags.(hs.subspaces, Val(:QuantumNumberAsTag))...))...)
    return sort(qns, by=reverse)
end

function get_tag(hs::ProductHilbertSpace{QN, N, S}, binrep::Unsigned, ::Val{:QuantumNumberAsTag}) where {QN, N, S}
    bvecs = bitsplit(bitwidth.(hs.subspaces), binrep)
    return get_tag.(hs.subspaces, bvecs, Val(:QuantumNumberAsTag))
end

function extract(hs::ProductHilbertSpace, binrep::Unsigned)
    bvecs = bitsplit(bitwidth.(hs.subspaces), binrep)
    indices = ((x,y) -> collect(extract(x,y).I)).(hs.subspaces, bvecs)
    return CartesianIndex(vcat(indices...)...)
end


# COV_EXCL_START
function _arraysplit(sizes::NTuple{N, <:Integer}, array::Vector{T}) where {N, T}
    if sum(sizes) != length(array)
        throw(ArgumentError("size mismatch"))
    end
    out = Vector{T}[]
    i = 1
    for s in sizes
        push!(out, array[i:i+s-1])
        i += s
    end
    return tuple(out...)
end
# COV_EXCL_STOP


function compress(
    hs::ProductHilbertSpace,
    indexarray::CartesianIndex,
    ::Type{BR}=UInt
) where {BR<:Unsigned}
    ns = numsites.(hs.subspaces)
    splitindices = map(x -> CartesianIndex(x...), _arraysplit(ns, collect(indexarray.I)))
    bvecs = ((x,y) -> compress(x, y, BR)).(hs.subspaces, splitindices)
    return bitjoin(bitwidth.(hs.subspaces), bvecs)
end
