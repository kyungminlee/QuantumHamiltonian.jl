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

bitwidth(lhs::AbstractHilbertSpaceRepresentation) = bitwidth(basespace(lhs))
