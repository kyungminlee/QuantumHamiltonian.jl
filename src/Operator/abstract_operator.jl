export AbstractOperator
export get_row_iterator, get_column_iterator, get_iterator
export scalartype
export bintype


import LinearAlgebra


"""
    AbstractOperator{S<:Number}

Represent an abstract operator in Hilbert space.
"""
abstract type AbstractOperator{S<:Number} end


"""
    scalartype(lhs::Type{<:AbstractOperator{S}})

Returns the scalar type of the given AbstractOperator.
"""
scalartype(::Type{<:AbstractOperator{S}}) where S = S
scalartype(::AbstractOperator{S}) where S = S


"""
    valtype(lhs::Type{<:AbstractOperator{S}})

Returns the `valtype` (scalar type) of the given AbstractOperator.
"""
Base.valtype(::Type{<:AbstractOperator{S}}) where S = S
Base.valtype(::AbstractOperator{S}) where S = S


bintype(lhs::AbstractOperator{S}) where S = bintype(typeof(lhs))::DataType


Base.:(-)(lhs::AbstractOperator, rhs::AbstractOperator) = (lhs) + (-rhs)
Base.:(-)(lhs, rhs::AbstractOperator) = (lhs) + (-rhs)
Base.:(-)(lhs::AbstractOperator, rhs) = (lhs) + (-rhs)

Base.:(+)(op::AbstractOperator) = op


function LinearAlgebra.issymmetric(arg::AbstractOperator{S}) where S
    return isa(simplify(arg - transpose(arg)), NullOperator)
end


function LinearAlgebra.ishermitian(arg::AbstractOperator{S}) where S
    return isa(simplify(arg - adjoint(arg)), NullOperator)
end


function Base.:(^)(lhs::AbstractOperator{S}, p::Integer) where S
    p <= 0 && error("Non-positive power for AbstractOperator not supported")

    # smallest nonzero power
    pow = simplify(lhs)
    while (p & 0x1) == 0
        pow = simplify(pow * pow)
        if isa(pow, NullOperator)
            return pow
        end
        p = p >> 1
    end

    out = pow
    p = p >> 1
    while p > 0
        pow = simplify(pow * pow)
        if (p & 0b1) != 0
            out = simplify(out * pow)
        end
        if isa(out, NullOperator)
            return out
        end
        p = p >> 1
    end
    return out
end
