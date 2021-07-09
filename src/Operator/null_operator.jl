export NullOperator


"""
    NullOperator

A null operator, i.e. 0.
"""
struct NullOperator<:AbstractOperator{Bool} end


bintype(::Type{<:NullOperator}) = UInt8


Base.:(-)(op::NullOperator) = op

Base.:(*)(::NullOperator, ::NullOperator) = NullOperator()

Base.:(*)(::AbstractOperator, ::NullOperator) = NullOperator()
Base.:(*)(::NullOperator, ::AbstractOperator) = NullOperator()

Base.:(*)(::Number, ::NullOperator)::NullOperator = NullOperator()
Base.:(*)(::NullOperator, ::Number)::NullOperator = NullOperator()
Base.:(\)(::Number, ::NullOperator)::NullOperator = NullOperator()
Base.:(/)(::NullOperator, ::Number)::NullOperator = NullOperator()
Base.:(//)(::NullOperator, ::Number)::NullOperator = NullOperator()

Base.:(+)(::NullOperator, ::NullOperator) = NullOperator()
Base.:(+)(lhs::AbstractOperator, ::NullOperator) = lhs
Base.:(+)(::NullOperator, rhs::AbstractOperator) = rhs

Base.:(==)(lhs::NullOperator, rhs::NullOperator) = true

Base.zero(::Type{NullOperator}) = NullOperator()
Base.zero(::NullOperator) = NullOperator()
Base.iszero(::NullOperator) = true

Base.real(arg::NullOperator) = arg
Base.imag(arg::NullOperator) = arg
Base.conj(arg::NullOperator) = arg
Base.transpose(arg::NullOperator) = arg
Base.adjoint(arg::NullOperator) = arg

# null operator is less than any other operators
Base.isless(lhs::NullOperator, rhs::NullOperator) = false
Base.isless(lhs::NullOperator, rhs::AbstractOperator) = true
Base.isless(lhs::AbstractOperator, rhs::NullOperator) = false
