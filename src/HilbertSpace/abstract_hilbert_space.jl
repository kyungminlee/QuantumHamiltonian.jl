export AbstractQuantumNumber
export AbstractHilbertSpace

## TODO: Think about this
AbstractQuantumNumber = Integer

abstract type AbstractHilbertSpace{QN<:Tuple{Vararg{<:AbstractQuantumNumber}}} end

qntype(::T) where {T<:AbstractHilbertSpace} = qntype(T)
tagtype(::T) where {T<:AbstractHilbertSpace} = tagtype(T)