export AbstractQuantumNumber
export AbstractHilbertSpace

## TODO: Think about this
AbstractQuantumNumber = Integer

abstract type AbstractHilbertSpace{QN<:Tuple{Vararg{<:AbstractQuantumNumber}}} end

scalartype(::T) where {T<:AbstractHilbertSpace} = scalartype(T)
Base.valtype(::T) where {T<:AbstractHilbertSpace} = valtype(T)
qntype(::T) where {T<:AbstractHilbertSpace} = qntype(T)
