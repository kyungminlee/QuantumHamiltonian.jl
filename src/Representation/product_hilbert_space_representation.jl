export ProductHilbertSpaceRepresentation
export get_basis, get_basis_list, get_basis_iterator, get_basis_state, get_basis_index_amplitude

struct ProductHilbertSpaceRepresentation{
    BR<:Unsigned,
    S<:Number,
    N,
    HS<:ProductHilbertSpace{<:Any, N, <:Any},
    HSR<:Tuple{Vararg{AbstractHilbertSpaceRepresentation, N}},
}<:AbstractHilbertSpaceRepresentation{BR, S}
    hilbertspace::HS
    subrepresentations::HSR
    dimensions::NTuple{N, Int}
    bitwidths::NTuple{N, Int}
    indices::CartesianIndices{N, NTuple{N, Base.OneTo{Int}}}

    function ProductHilbertSpaceRepresentation(
        hilbertspace::HS,
        subrepresentations::HSR,
        ::Type{BR}=UInt
    ) where {
        BR<:Unsigned,
        N,
        HS<:ProductHilbertSpace{<:Any, N, <:Any},
        HSR<:Tuple{Vararg{AbstractHilbertSpaceRepresentation, N}},
    }
        dimensions = dimension.(subrepresentations)
        bitwidths = bitwidth.(hilbertspace.subspaces)
        indices = CartesianIndices(tuple(dimensions...))
        S = promote_type(scalartype.(subrepresentations)...)
        return new{BR, S, N, HS, HSR}(hilbertspace, subrepresentations, dimensions, bitwidths, indices)
    end

    function ProductHilbertSpaceRepresentation(hsr_list::AbstractHilbertSpaceRepresentation...)
        space = ProductHilbertSpace(basespace.(hsr_list))
        return ProductHilbertSpaceRepresentation(space, hsr_list)
    end
end

Base.valtype(::Type{<:ProductHilbertSpaceRepresentation{BR, S, N, HS, HSR}}) where {BR, S, N, HS, HSR} = S
scalartype(::Type{<:ProductHilbertSpaceRepresentation{BR, S, N, HS, HSR}}) where {BR, S, N, HS, HSR} = S
bintype(::Type{<:ProductHilbertSpaceRepresentation{BR, S, N, HS, HSR}}) where {BR, S, N, HS, HSR} = BR

qntype(::Type{<:ProductHilbertSpaceRepresentation{BR, S, N, HS, HSR}}) where {BR, S, N, HS, HSR} = qntype(HS)

basespace(hsr::ProductHilbertSpaceRepresentation) = hsr.hilbertspace
basespace(::Type{<:ProductHilbertSpaceRepresentation{BR, S, N, HS, HSR}}) where {BR, S, N, HS, HSR} = HS

dimension(p::ProductHilbertSpaceRepresentation) = mapreduce(dimension, *, p.subrepresentations)


function get_basis_state(hsr::ProductHilbertSpaceRepresentation{BR, S, N, HS, HSR}, index::Integer) where {BR, S, N, HS, HSR}
    bvecs = ((x,y) -> get_basis_state(x, y)).(hsr.subrepresentations, hsr.indices[index].I)
    return bitjoin(hsr.bitwidths, bvecs, BR)
end


function get_basis_list(hsr::ProductHilbertSpaceRepresentation{BR, S, N, HS, HSR}) where {BR, S, N, HS, HSR}
    out = Vector{BR}(undef, dimension(hsr))
    i = 0
    for x in Iterators.product(get_basis_list.(hsr.subrepresentations)...)
        i += 1
        out[i] = bitjoin(hsr.bitwidths, x)
    end
    return out
end

function get_basis_iterator(hsr::ProductHilbertSpaceRepresentation{BR, S, N, HS, HSR}) where {BR, S, N, HS, HSR}
    return Iterators.flatten(
        bitjoin(hsr.bitwidths, x)
            # for x in Iterators.product(get_basis_list.(hsr.subrepresentations)...)
            for x in Iterators.product( (get_basis_list(z) for z in hsr.subrepresentations)...)
    )
end


function get_basis_index_amplitude(hsr::ProductHilbertSpaceRepresentation{BR, S, N, HS, HSR}, bvec::Unsigned) where {BR, S, N, HS, HSR}
    iszero((~get_bitmask(hsr)) & bvec) || return (index=-1, amplitude=zero(S))
    bvecs = bitsplit(hsr.bitwidths, bvec)
    iat = get_basis_index_amplitude.(hsr.subrepresentations, bvecs)
    indices = CartesianIndex(map(x->x[1], iat))
    amplitude = mapreduce(x -> x[2], *, iat)
    if any(x -> x <= 0, indices.I)
        return (index=-1, amplitude=zero(S))
    end
    return (index=LinearIndices(hsr.indices)[indices], amplitude=S(amplitude))
end

