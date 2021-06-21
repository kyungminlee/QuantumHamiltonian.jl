export product_state

function product_state(
    hs::HilbertSpace,
    local_states::AbstractVector,
    ::Type{BR}=UInt
) where BR
    return product_state(hs, (local_states...,), BR)
end


function product_state(
    hs::HilbertSpace,
    local_states::Tuple{Vararg{<:AbstractVector{<:Number}}},
    ::Type{BR}=UInt
) where BR
    n_sites = length(hs.sites)
    if length(local_states) != n_sites
        throw(ArgumentError("length of local states does not match the number of sites in the Hilbert space"))
    end
    for (isite, site) in enumerate(hs.sites)
        if length(site.states) != length(local_states[isite])
            throw(ArgumentError("length of the local state at site $isite does not match the dimension of the local Hilbert site"))
        end
    end

    S = promote_type(eltype.(typeof.(local_states))...)
    out = SparseState{S, BR}()

    index_iterator = Iterators.product((
        collect(
            istate
                for istate in 1:length(site.states) if !iszero(local_states[isite][istate])
        )
            for (isite, site) in enumerate(hs.sites)
    )...,)
    for indexarray in index_iterator
        bvec = compress(hs, CartesianIndex(indexarray), BR)
        coeff = prod(local_states[i_site][indexarray[i_site]] for i_site in Base.OneTo(n_sites))
        out[bvec] = coeff
    end
    return out
end


function product_state(hsr::HilbertSpaceRepresentation, local_states::AbstractVector)
    return product_state(hsr, (local_states...,))
end


function product_state(
    hsr::HilbertSpaceRepresentation{HS, BR, DT},
    local_states::Tuple{Vararg{<:AbstractVector{<:Number}}}
) where {HS, BR, DT}
    hs = basespace(hsr)
    n_sites = length(hs.sites)
    if length(local_states) != n_sites
        throw(ArgumentError("length of local states does not match the number of sites in the Hilbert space"))
    end
    for (isite, site) in enumerate(hs.sites)
        if length(site.states) != length(local_states[isite])
            throw(ArgumentError("length of the local state at site $isite does not match the dimension of the local Hilbert site"))
        end
    end

    S = promote_type(eltype.(typeof.(local_states))...)
    out = zeros(S, dimension(hsr))

    index_iterator = Iterators.product((
        collect(
            istate
                for istate in 1:length(site.states) if !iszero(local_states[isite][istate])
        )
            for (isite, site) in enumerate(hs.sites)
    )...,)
    for indexarray in index_iterator
        bvec = compress(hs, CartesianIndex(indexarray), BR)
        # idx = get(hsr.basis_lookup, bvec, -1)
        idx = findindex(hsr.basis, bvec)
        idx <= 0 && continue
        coeff = prod(local_states[i_site][indexarray[i_site]] for i_site in Base.OneTo(n_sites))
        out[idx] = coeff
    end
    return out
end
