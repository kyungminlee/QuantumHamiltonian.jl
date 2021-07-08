export DecomposedHilbertSpaceRepresentation
export dimension
export get_basis_list, get_basis_state, get_basis_index_amplitude
export sectorslice


"""
    DecomposedHilbertSpaceRepresentation

Represents the decomposed Hilbert space representation.
i.e. subspace of the total Hilbert space, with basis states grouped by tags.
"""
struct DecomposedHilbertSpaceRepresentation{
    BR<:Unsigned,
    S<:Number,
    HS<:AbstractHilbertSpace,
    HSR<:AbstractHilbertSpaceRepresentation,
    TagStrategy,
    TagType,
}<:AbstractHilbertSpaceRepresentation{BR, S}
    hilbertspace::HS
    tags::DictIndexedVector{TagType}
    components::Vector{HSR}
    offsets::Vector{Int}

    function DecomposedHilbertSpaceRepresentation(
        hs::HS,
        tags::AbstractVector{TT},
        components::AbstractVector{HSR},
        tagstrategy::Val{TS}
    ) where {HS<:AbstractHilbertSpace, TT, S, BR, HSR<:AbstractHilbertSpaceRepresentation{BR, S}, TS}
        if length(tags) != length(components)
            throw(ArgumentError("number of tags should match number of components"))
        end
        if tagtype(HS, tagstrategy) != TT
            throw(ArgumentError("Tag strategy not compatible with tags ($(tagtype(HS, tagstrategy)) != $(TT))"))
        end
        offsets = Vector{Int}(undef, length(components)+1)
        offsets[1] = 0
        for (icompo, compo) in enumerate(components)
            offsets[icompo+1] = offsets[icompo] + dimension(compo)
        end
        return new{BR, S, HS, HSR, TS, TT}(hs, DictIndexedVector(tags), components, offsets)
    end

end




# specializations
function tagtype(
    ::Type{<:DecomposedHilbertSpaceRepresentation{BR, S, HS, HSR, TS, TT}},
    ::Val{TS2}=Val(TS)
) where {BR, S, HS, HSR, TS, TT, TS2}
    TS == TS2 || throw(ArgumentError("tag strategies do not match ($TS != $TS2)"))
    return TT
end


function tagtype(
    ::DecomposedHilbertSpaceRepresentation{BR, S, HS, HSR, TS, TT},
    ::Val{TS2}=Val(TS)
) where {BR, S, HS, HSR, TS, TT, TS2}
    TS == TS2 || throw(ArgumentError("tag strategies do not match ($TS != $TS2)"))
    return TT
end


function get_tag(
    hsr::DecomposedHilbertSpaceRepresentation{BR, S, HS, HSR, TS, TT},
    binrep::Unsigned,
    ::Val{TS2}=Val(TS),
) where {BR, S, HS, HSR, TS, TT, TS2}
    TS == TS2 || throw(ArgumentError("tag strategies do not match ($TS != $TS2)"))
    return get_tag(hsr.hilbertspace, binrep, Val(TS))
end


function checkvalid(dhsr::DecomposedHilbertSpaceRepresentation{BR, S, HS, HSR, TS, TT}) where {BR, S, HS, HSR, TS, TT}
    for (tag, compo) in zip(dhsr.tags, dhsr.components)
        bl = get_basis_list(compo)
        for b in bl
            @assert tag == get_tag(dhsr.hilbertspace, b, Val(TS))
        end
    end
end


dimension(hsr::DecomposedHilbertSpaceRepresentation) = hsr.offsets[end] #mapreduce(length, +, hsr.components)


get_basis_list(hsr::DecomposedHilbertSpaceRepresentation) = vcat(get_basis_list.(hsr.components)...)
get_basis_iterator(hsr::DecomposedHilbertSpaceRepresentation) = Iterators.flatten(get_basis_list(x) for x in hsr.components)


function get_basis_state(hsr::DecomposedHilbertSpaceRepresentation, index::Integer)
    (1 <= index <= dimension(hsr)) || throw(BoundsError(hsr, index))
    # @show index, hsr.offsets
    index -= 1
    c = searchsortedlast(hsr.offsets, index)
    index -= hsr.offsets[c]
    return get_basis_state(hsr.components[c], index+1)
end


function get_basis_index_amplitude(hsr::DecomposedHilbertSpaceRepresentation{BR, S, HS, HSR, TS, TT}, bvec::Unsigned) where {BR, S, HS, HSR, TS, TT}
    try
        tag = get_tag(hsr, bvec, Val(TS))
        icompo = findindex(hsr.tags, tag)
        icompo > 0 || return (index=-1, amplitude=zero(S))
        (i, a) = get_basis_index_amplitude(hsr.components[icompo], bvec)
        i > 0 || return (index=-1, amplitude=zero(S))
        return (index=hsr.offsets[icompo] + i, amplitude=a)
    catch e
        if isa(e, BoundsError)
            return (index=-1, amplitude=zero(S))
        else
            rethrow() # COV_EXCL_LINE
        end
    end
end


basespace(hsr::DecomposedHilbertSpaceRepresentation) = hsr.hilbertspace


"""
    Predicate on the tags
"""
function sectorslice(
    hsr::DecomposedHilbertSpaceRepresentation{BR, S, HS, HSR, TS, TT},
    predicate::Function
) where {BR, S, HS, HSR, TS, TT}
    indices = findall(predicate, hsr.tags)
    tags = hsr.tags[indices]
    components = hsr.components[indices]
    DecomposedHilbertSpaceRepresentation(hsr.hilbertspace, tags, components, Val(TS))
end


"""
    Predicate on the tags
"""
function sectorslice(
    hsr::DecomposedHilbertSpaceRepresentation{BR, S, HS, HSR, TS, TT},
    tags::AbstractVector{TT}
) where {BR, S, HS, HSR, TS, TT}
    indices = findall(x -> x in tags, hsr.tags)
    tags_sel = tags[indices]
    components = hsr.components[indices]
    DecomposedHilbertSpaceRepresentation(hsr.hilbertspace, tags_sel, components, Val(TS))
end