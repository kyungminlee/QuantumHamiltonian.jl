export ReducedHilbertSpaceRepresentation
export bintype
export symmetry_reduce, symmetry_unreduce
export get_basis_state, get_basis_index_amplitude

"""
    ReducedHilbertSpaceRepresentation{HSR, BR, C}

Representation of the symmetry-reduced hilbert space.
Currently only supports Translation group (i.e. Abelian group).
```
"""
struct ReducedHilbertSpaceRepresentation{
    HSR<:HilbertSpaceRepresentation,
    BR<:Unsigned,
    C<:Number,
    IV<:AbstractIndexedVector{BR}
}<:AbstractHilbertSpaceRepresentation{BR, C}
    parent::HSR
    basis::IV
    basis_mapping_index::Vector{Int} # has size of parent dimension. each index item contains index at reduced basis, or -1 if not included
    basis_mapping_amplitude::Vector{C}

    function ReducedHilbertSpaceRepresentation(
        hsr::HSR,
        basis::IV,
        basis_mapping_index::Vector{Int},
        basis_mapping_amplitude::Vector{C}
    ) where {HSR<:HilbertSpaceRepresentation, BR<:Unsigned, C<:Number, IV<:AbstractIndexedVector{BR}}
        return new{HSR, BR, C, IV}(hsr, basis, basis_mapping_index, basis_mapping_amplitude)
    end
end


Base.valtype(::Type{ReducedHilbertSpaceRepresentation{HSR, BR, C, IV}}) where {HSR, BR, C, IV} = C
scalartype(::Type{ReducedHilbertSpaceRepresentation{HSR, BR, C, IV}}) where {HSR, BR, C, IV} = C
bintype(::Type{ReducedHilbertSpaceRepresentation{HSR, BR, C, IV}}) where {HSR, BR, C, IV} = BR


basespace(lhs::ReducedHilbertSpaceRepresentation) = basespace(lhs.parent)


@inline function get_basis_state(hsr::ReducedHilbertSpaceRepresentation, index::Integer)
    return hsr.basis[index]
end

@inline function get_basis_index_amplitude(hsr::ReducedHilbertSpaceRepresentation{HSR, BR, C, IV}, bvec::Unsigned) where {HSR, BR, C, IV}
    i, a = get_basis_index_amplitude(hsr.parent, bvec)
    if i <= 0
        return (index=i, amplitude=zero(C))
    else
        return (index=hsr.basis_mapping_index[i], amplitude=hsr.basis_mapping_amplitude[i] * a)
    end
end


"""
    dimension(arg::ReducedHilbertSpaceRepresentation{HSR, BR, C}) -> Int

Dimension of the given reduced hilbert space representation, i.e. number of basis elements.
"""
dimension(arg::ReducedHilbertSpaceRepresentation) = length(arg.basis)
