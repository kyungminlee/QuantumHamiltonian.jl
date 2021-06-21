export HilbertSpaceRepresentation
export dimension
export represent, represent_array, represent_dict
export get_basis, get_basis_list, get_basis_state, get_basis_index_amplitude

@inline represent(::Site, istate::Integer, ::Type{BR}=UInt) where {BR} = BR(istate-1)

"""
    HilbertSpaceRepresentation{HS, BR, DictType}

# Fields
```
hilbert_space :: HS
basis_list    :: Vector{BR}
basis_lookup  :: DictType
```
"""
struct HilbertSpaceRepresentation{
    HS<:AbstractHilbertSpace,
    BR<:Unsigned,
    BasisType<:AbstractIndexedVector{BR}
}<:AbstractHilbertSpaceRepresentation{Bool}

    hilbert_space::HS
    basis::BasisType

    function HilbertSpaceRepresentation(
        hilbert_space::AbstractHilbertSpace,
        basis::BT,
    ) where {BR<:Unsigned, BT<:AbstractIndexedVector{BR}}
        if sizeof(BR)*8 <= bitwidth(hilbert_space)
            # equality added such that the MSB checks overflow
            throw(ArgumentError(
                "type $(BR) not enough to represent the hilbert space"*
                " (need $(bitwidth(hilbert_space)) bits)"
            ))
        end
        hs = basespace(hilbert_space)
        HS = typeof(hs)
        return new{HS, BR, BT}(hs, basis)
    end
end


Base.valtype(::Type{HilbertSpaceRepresentation{HS, BR, BT}}) where {HS, BR, BT} = Bool
scalartype(::Type{HilbertSpaceRepresentation{HS, BR, BT}}) where {HS, BR, BT} = Bool
bintype(::Type{HilbertSpaceRepresentation{HS, BR, BT}}) where {HS, BR, BT} = BR


basespace(lhs::HilbertSpaceRepresentation{HS, BR, BT}) where {HS, BR, BT} = lhs.hilbert_space::HS


"""
    dimension

Dimension of the Concrete Hilbert space, i.e. number of basis vectors.
"""
dimension(hsr::HilbertSpaceRepresentation) = length(hsr.basis)


function Base.:(==)(
    lhs::HilbertSpaceRepresentation{HS1, BR1, BT1},
    rhs::HilbertSpaceRepresentation{HS2, BR2, BT2},
) where {HS1, BR1, BT1, HS2, BR2, BT2}
    return basespace(lhs) == basespace(rhs) && (lhs.basis == rhs.basis)
end

@inline get_basis(hsr::HilbertSpaceRepresentation) = hsr.basis

@inline get_basis_list(hsr::HilbertSpaceRepresentation) = collect(hsr.basis)

"""
    get_basis_state(hsr, index)

Get the basis state representation at index.
"""
@inline function get_basis_state(hsr::HilbertSpaceRepresentation, index::Integer)
    return hsr.basis[index]
end


"""
    get_basis_index_amplitude(hsr, bvec)

Get the index of the basis state that overlaps with bvec, and the value of the overlap.
Currentiy it is guaranteed to be at most one.
Returns (i, ⟨b|ϕᵢ⟩). For the unsymmetrized `HilbertSpaceRepresentation`, the amplitude is
1 of Int type.
If no such basis vector exists, return (-1, 0).
"""
@inline function get_basis_index_amplitude(hsr::HilbertSpaceRepresentation, bvec::Unsigned)
    index = findindex(hsr.basis, bvec)
    # index = get(hsr.basis_lookup, bvec, -1)
    return (index=index, amplitude=(index <= 0) ? 0 : 1)
end




function checkvalidbasis(hsr::HilbertSpaceRepresentation{HS, BR, BT}) where {HS, BR, BT}
    for (ivec, bvec) in enumerate(hsr.basis)
        # ivec2 = hsr.basis_lookup[bvec]
        ivec2 = findindex(hsr.basis, bvec)
        @assert ivec == ivec2
    end
end




"""
    represent(hs, binary_type=UInt)

Make a HilbertSpaceRepresentation with all the basis vectors of the specified HilbertSpace.
This function defaults to [`represent_array`](@ref).
"""
function represent(hs::AbstractHilbertSpace, ::Type{BR}=UInt) where {BR <:Unsigned}
    return represent_array(hs, BR)
end


"""
    represent_array(hs, binary_type=UInt)

Make a HilbertSpaceRepresentation with all the basis vectors of the specified HilbertSpace
using [`FrozenSortedArrayIndex`](@ref).
"""
function represent_array(hs::AbstractHilbertSpace, ::Type{BR}=UInt) where {BR <:Unsigned}
    basis_list = hs_get_basis_list(hs, BR)
    # basis_lookup = FrozenSortedArrayIndex{BR}(basis_list)
    basis = SortedIndexedVector(basis_list)
    return HilbertSpaceRepresentation(basespace(hs), basis)
end


"""
    represent(hs, basis_list)

Make a HilbertSpaceRepresentation with the provided list of basis vectors.
This defaults to [`represent_array`](@ref).
"""
function represent(hs::AbstractHilbertSpace, basis_list::AbstractVector{BR}) where {BR<:Unsigned}
    return represent_array(hs, basis_list)
end


"""
    represent_array(hs, basis_list)

Make a HilbertSpaceRepresentation with the provided list of basis vectors
using [`FrozenSortedArrayIndex`](@ref).
"""
function represent_array(hs::AbstractHilbertSpace, basis_list::AbstractVector{BR}) where {BR<:Unsigned}
    # if !issorted(basis_list)
    #     basis_list = sort(basis_list)
    # end
    # basis_lookup = FrozenSortedArrayIndex{BR}(basis_list)
    basis = SortedIndexedVector(basis_list)
    return HilbertSpaceRepresentation(basespace(hs), basis)
end


"""
    represent_dict(hs, binary_type=UInt)

Make a HilbertSpaceRepresentation with the provided list of basis vectors
using `Dict{binary_type, Int}`.
"""
function represent_dict(hs::AbstractHilbertSpace, binary_type::Type{BR}=UInt) where {BR<:Unsigned}
    basis_list = hs_get_basis_list(hs, BR)
    basis = DictIndexedVector(basis_list)
    return HilbertSpaceRepresentation(basespace(hs), basis)
end


"""
    represent_dict(hs, basis_list)

Make a HilbertSpaceRepresentation with the provided list of basis vectors using `Dict`.
"""
function represent_dict(hs::AbstractHilbertSpace, basis_list::AbstractVector{BR}) where {BR<:Unsigned}
    basis = DictIndexedVector(basis_list)
    return HilbertSpaceRepresentation(basespace(hs), basis)
end
