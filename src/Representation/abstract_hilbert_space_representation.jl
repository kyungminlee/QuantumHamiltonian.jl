export AbstractHilbertSpaceRepresentation
export scalartype, bintype


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