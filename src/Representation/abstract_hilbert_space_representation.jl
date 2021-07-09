export AbstractHilbertSpaceRepresentation
export scalartype, bintype, basespace

export get_basis_list, get_basis_iterator, get_basis_state, get_basis_index_amplitude

"""
    AbstractHilbertSpaceRepresentation{BR, S}
"""
abstract type AbstractHilbertSpaceRepresentation{BR<:Unsigned, S<:Number} end

Base.valtype(::T) where {T<:AbstractHilbertSpaceRepresentation} = valtype(T)
scalartype(::T) where {T<:AbstractHilbertSpaceRepresentation} = scalartype(T)
bintype(::T) where {T<:AbstractHilbertSpaceRepresentation} = bintype(T)

Base.valtype(::Type{<:AbstractHilbertSpaceRepresentation{BR, S}}) where {BR, S} = S
scalartype(::Type{<:AbstractHilbertSpaceRepresentation{BR, S}}) where {BR, S} = S
bintype(::Type{<:AbstractHilbertSpaceRepresentation{BR, S}}) where {BR, S} = BR


"""
    basespace(x::AbstractHilbertSpaceRepresentation)
    basespace(x::Type{T}) where {T<:AbstractHilbertSpaceRepresentation}

If the argument is an instance of `AbstractHilbertSpaceRepresentation`, return the underlying Hilbert space of the Hilbert space representation.
If the argument is a subtype of `AbstractHilbertSpaceRepresentation`, return the type of the underlying Hilbert space.
Subtypes of AbstractHilbertSpace must implement this method.
"""
function basespace end


"""
    get_basis_list(hsr::AbstractHilbertSpaceRepresentation)

Return a Vector of the list of basis binaries.
Subtypes of AbstractHilbertSpaceRepresentation must implement this method.
"""
function get_basis_list end


"""
    get_basis_iterator(hsr::AbstractHilbertSpaceRepresentation)

Return an iterator of the list of basis binaries.
Subtypes of AbstractHilbertSpaceRepresentation must implement this method.
"""
function get_basis_iterator end


"""
    get_basis_state(hsr::AbstractHilbertSpaceRepresentation, index::Integer)

Return the state at index `index`.
Subtypes of AbstractHilbertSpaceRepresentation must implement this method.
"""
function get_basis_state end


"""
    get_basis_index_amplitude(hsr::AbstractHilbertSpaceRepresentation, bin::Unsigned)

Return a tuple `(index=index, amplitude=amplitude)` that corresponds to the binary `bin`,
i.e., the index of the basis state that overlaps with `bin` and the value of the overlap ⟨b|ϕᵢ⟩.
Return `(index=-1, amplitude=0)` if `bin` is not in the representation.
Subtypes of AbstractHilbertSpaceRepresentation must implement this method.
"""
function get_basis_index_amplitude end


for fname in [
    :numsites,
    :get_site,
    :bitwidth,
    :bitoffset,
    :get_quantum_number,
    :get_tag,
    :extract,
    :compress,
    :uncompress,
    :update,
    :get_state,
    :get_state_index,
    :get_bitmask,
]
    @eval begin
        """
            $($fname)(hsr::AbstractHilbertSpaceRepresentation, args...;kwargs...)

        Call `$($fname)` with basespace of `hsr`.
        """
        @inline $fname(hsr::AbstractHilbertSpaceRepresentation, args...; kwargs...) = $fname(basespace(hsr), args...; kwargs...)
        
    end
end


for fname in [
    :qntype,
    :tagtype,
]
    @eval begin
        """
            $($fname)(hsr::AbstractHilbertSpaceRepresentation, args...;kwargs...)

        Call `$($fname)` with basespace of `hsr`.
        """
        @inline $fname(hsr::AbstractHilbertSpaceRepresentation, args...; kwargs...) = $fname(basespace(hsr), args...; kwargs...)
        
        @inline $fname(::Type{HSR}, args...; kwargs...) where {HSR<:AbstractHilbertSpaceRepresentation} = $fname(basespace(HSR), args...; kwargs...)
    end
end