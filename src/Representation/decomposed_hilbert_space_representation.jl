export DecomposedHilbertSpaceRepresentation
export represent_decompose_dict
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

    function DecomposedHilbertSpaceRepresentation(
        hsr::HilbertSpaceRepresentation{BR, HS, BT},
        tagstrategy::Val{TS}=Val(:QuantumNumberAsTag)
    ) where {BR, HS, BT, TS}
        S = Bool
        HSR = HilbertSpaceRepresentation{BR, HS, BT}
        TT = tagtype(HS, tagstrategy)

        hs = basespace(hsr)
        tag_basis_list = Dict{TT, Vector{BR}}()
        for b in get_basis_list(hsr)
            tag = get_tag(hs, b, tagstrategy)
            if haskey(tag_basis_list, tag)
                push!(tag_basis_list[tag], b)
            else
                tag_basis_list[tag] = [b]
            end
        end
        tags = DictIndexedVector(unique(sort(collect(keys(tag_basis_list)))))
        components = Vector{HSR}(undef, length(tags))
        offsets = Vector{Int}(undef, length(tags)+1)
        offsets[1] = 0
        for (itag, tag) in enumerate(tags)
            blist = HilbertSpaceRepresentation(hs, BT(tag_basis_list[tag]))
            components[itag] = blist
            offsets[itag+1] = offsets[itag] + dimension(blist)
        end
        return new{BR, S, HS, HSR, TS, TT}(hs, tags, components, offsets)
    end
end


# function represent_decompose_dict(hs::HilbertSpace{QN, TT}, ::Type{BR}=UInt) where {QN, TT, BR<:Unsigned}
#     S = Bool
#     HS = HilbertSpace{QN, TT}
#     HSR = HilbertSpaceRepresentation{BR, HS, DictIndexedVector{BR}}
#     TagType = TT
    
#     tag_basis_list = Dict{TagType, Vector{BR}}()
#     for b in QuantumHamiltonian.hs_get_basis_list(hs)
#         tag = get_tag(hs, b)
#         if haskey(tag_basis_list, tag)
#             push!(tag_basis_list[tag], b)
#         else
#             tag_basis_list[tag] = [b]
#         end
#     end
#     tags = unique(sort(collect(keys(tag_basis_list))))
#     components = Vector{HSR}(undef, length(tags))
#     offsets = Vector{Int}(undef, length(tags)+1)
#     offsets[1] = 0
#     for (itag, tag) in enumerate(tags)
#         blist = represent_dict(hs,tag_basis_list[tag])
#         components[itag] = blist
#         offsets[itag+1] = offsets[itag] + dimension(blist)
#     end
#     return DecomposedHilbertSpaceRepresentation{BR, S, HS, HSR, TagType}(hs, tags, components, offsets)
# end


dimension(hsr::DecomposedHilbertSpaceRepresentation) = hsr.offsets[end] #mapreduce(length, +, hsr.components)

# get_tags(hsr::DecomposedHilbertSpaceRepresentation) = sort(collect(hsr.tags))

get_basis_list(hsr::DecomposedHilbertSpaceRepresentation) = vcat(get_basis_list.(hsr.components)...)


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
        # hsr.tags[icompo] == tag || return (index=-1, amplitude=zero(S))
        (i, a) = get_basis_index_amplitude(hsr.components[icompo], bvec)
        i > 0 || return (index=-1, amplitude=zero(S))
        return (index=hsr.offsets[icompo] + i, amplitude=a)
    catch e
        if isa(e, BoundsError)
            return (index=-1, amplitude=zero(S))
        else
            rethrow()
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
    DecomposedHilbertSpaceRepresentation(hsr.hilbertspace, tags_sel, components)
end