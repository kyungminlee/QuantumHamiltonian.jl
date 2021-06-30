export decompose

function decompose(
    hsr::HilbertSpaceRepresentation{BR, HS, BT},
    tagstrategy::Val{TS}=Val(:QuantumNumberAsTag)
) where {BR, HS, BT, TS}
    HSR = HilbertSpaceRepresentation{BR, HS, BT}
    TT = tagtype(HSR, tagstrategy)
    
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
    for (itag, tag) in enumerate(tags)
        blist = HilbertSpaceRepresentation(hs, BT(tag_basis_list[tag]))
        components[itag] = blist
    end
    return DecomposedHilbertSpaceRepresentation(hs, tags, components, tagstrategy)
end

function decompose(
    hsr::ProductHilbertSpaceRepresentation{BR, S, N, HS, HSR},
    tagstrategy::Val{TS}=Val(:QuantumNumberAsTag)
) where {BR, S, N, HS, HSR, TS}
    hs = basespace(hsr)
    subs = decompose.(hsr.subrepresentations, tagstrategy)
    size_list = (x -> length(x.components)).(subs)
    components = []
    tags = []
    for ci in CartesianIndices(size_list)
        tag = tuple([sub.tags[i] for (i, sub) in zip(ci.I, subs)]...)
        compo = tuple([sub.components[i] for (i, sub) in zip(ci.I, subs)]...)
        push!(tags, tag)
        push!(components, ProductHilbertSpaceRepresentation(hs, compo))
    end
    tags = [tags...]
    components = [components...]
    return DecomposedHilbertSpaceRepresentation(hs, tags, components, tagstrategy)
end