export HilbertSpaceRepresentation
export basespace
export dimension
export get_basis, get_basis_list, get_basis_state, get_basis_index_amplitude
export represent, represent_array, represent_dict

# @inline represent(::Site, istate::Integer, ::Type{BR}=UInt) where {BR} = BR(istate-1)

"""
    HilbertSpaceRepresentation{BR, HS, BasisType}

"""
struct HilbertSpaceRepresentation{
    BR<:Unsigned,
    HS<:AbstractHilbertSpace,
    BasisType<:AbstractIndexedVector{BR},
}<:AbstractHilbertSpaceRepresentation{BR, Bool}

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
        return new{BR, HS, BT}(hs, basis)
    end
end


basespace(lhs::HilbertSpaceRepresentation) = lhs.hilbert_space
basespace(::Type{<:HilbertSpaceRepresentation{BR, HS, BT}}) where {BR, HS, BT} = HS


"""
    dimension

Dimension of the Concrete Hilbert space, i.e. number of basis vectors.
"""
dimension(hsr::HilbertSpaceRepresentation) = length(hsr.basis)


function Base.:(==)(lhs::HilbertSpaceRepresentation, rhs::HilbertSpaceRepresentation)
    return basespace(lhs) == basespace(rhs) && (lhs.basis == rhs.basis)
end


get_basis(hsr::HilbertSpaceRepresentation) = hsr.basis


get_basis_list(hsr::HilbertSpaceRepresentation) = collect(hsr.basis)


"""
    get_basis_state(hsr, index)

Get the basis state representation at index.
"""
get_basis_state(hsr::HilbertSpaceRepresentation, index::Integer) = hsr.basis[index]


"""
    get_basis_index_amplitude(hsr, bvec)

Get the index of the basis state that overlaps with bvec, and the value of the overlap.
Currentiy it is guaranteed to be at most one.
Returns (i, ⟨b|ϕᵢ⟩). For the unsymmetrized `HilbertSpaceRepresentation`, the amplitude is
1 of Int type.
If no such basis vector exists, return (-1, 0).
"""
function get_basis_index_amplitude(hsr::HilbertSpaceRepresentation, bvec::Unsigned)
    index = findindex(hsr.basis, bvec)
    return (index=index, amplitude=(index <= 0) ? 0 : 1)
end


"""
    represent(hs, binary_type=UInt, basis_type=SortedIndexedVector)

Make a HilbertSpaceRepresentation with all the basis vectors of the specified HilbertSpace.
This function defaults to [`represent_array`](@ref).
"""
function represent(
    hs::AbstractHilbertSpace,
    ::Type{BR}=UInt,
    ::Type{BT}=SortedIndexedVector
) where {BR<:Unsigned, BT<:AbstractIndexedVector}
    basis_list = hs_get_basis_list(hs, BR)
    basis = BT(basis_list)
    return HilbertSpaceRepresentation(basespace(hs), basis)
end


function represent(
    hs::AbstractHilbertSpace{QN},
    allowed_quantum_numbers::AbstractVector{QN},
    ::Type{BR}=UInt,
    ::Type{BT}=SortedIndexedVector
) where {QN, BR, BT}
    basis_list = hs_get_basis_list(hs, allowed_quantum_numbers, BR)
    basis = BT(basis_list)
    return HilbertSpaceRepresentation(basespace(hs), basis)
end


function represent(
    hs::AbstractHilbertSpace{Tuple{I}},
    allowed_quantum_numbers::AbstractVector{I},
    ::Type{BR}=UInt,
    ::Type{BT}=SortedIndexedVector
) where {I<:AbstractQuantumNumber, BR, BT}
    return represent(hs, [(x,) for x in allowed_quantum_numbers], BR, BT)
end


"""
    represent(hs, basis_list, basis_type=SortedIndexedVector)

Make a HilbertSpaceRepresentation with the provided list of basis vectors.
This defaults to [`represent_array`](@ref).
"""
function represent(
    hs::AbstractHilbertSpace,
    basis_list::AbstractVector{BR},
    ::Type{BT}=SortedIndexedVector
) where {BR<:Unsigned, BT<:AbstractIndexedVector}
    basis = BT(basis_list)
    return HilbertSpaceRepresentation(basespace(hs), basis)
end


# COV_EXCL_START

"""
    represent_array(hs, binary_type=UInt)

Make a HilbertSpaceRepresentation with all the basis vectors of the specified HilbertSpace
using [`FrozenSortedArrayIndex`](@ref).
"""
function represent_array(hs::AbstractHilbertSpace, ::Type{BR}=UInt) where {BR<:Unsigned}
    @warn "represent_array deprecated. Use reprent(hs, BR, BT) instead." maxlog=1
    basis_list = hs_get_basis_list(hs, BR)
    basis = SortedIndexedVector(basis_list)
    return HilbertSpaceRepresentation(basespace(hs), basis)
end


"""
    represent_array(hs, basis_list)

Make a HilbertSpaceRepresentation with the provided list of basis vectors
using [`FrozenSortedArrayIndex`](@ref).
"""
function represent_array(hs::AbstractHilbertSpace, basis_list::AbstractVector{BR}) where {BR<:Unsigned}
    @warn "represent_array deprecated. Use reprent(hs, BR, BT) instead." maxlog=1
    basis = SortedIndexedVector(basis_list)
    return HilbertSpaceRepresentation(basespace(hs), basis)
end


"""
    represent_dict(hs, binary_type=UInt)

Make a HilbertSpaceRepresentation with the provided list of basis vectors
using `Dict{binary_type, Int}`.
"""
function represent_dict(hs::AbstractHilbertSpace, ::Type{BR}=UInt) where {BR<:Unsigned}
    @warn "represent_dict deprecated. Use reprent(hs, BR, BT) instead." maxlog=1
    basis_list = hs_get_basis_list(hs, BR)
    basis = DictIndexedVector(basis_list)
    return HilbertSpaceRepresentation(basespace(hs), basis)
end


"""
    represent_dict(hs, basis_list)

Make a HilbertSpaceRepresentation with the provided list of basis vectors using `Dict`.
"""
function represent_dict(hs::AbstractHilbertSpace, basis_list::AbstractVector{BR}) where {BR<:Unsigned}
    @warn "represent_dict deprecated. Use reprent(hs, BR, BT) instead." maxlog=1
    basis = DictIndexedVector(basis_list)
    return HilbertSpaceRepresentation(basespace(hs), basis)
end

# COV_EXCL_STOP
