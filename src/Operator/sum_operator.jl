export SumOperator
export bintype


"""
    SumOperator{Scalar, BR}

Represents a sum of pure operators.

# Fields
- `terms::Vector{PureOperator{Scalar,BR}}`
"""
struct SumOperator{Scalar<:Number, BR<:Unsigned} <:AbstractOperator{Scalar}
    terms::Vector{PureOperator{Scalar, BR}}

    function SumOperator{S, BR}(terms) where {S, BR}
        return new{S, BR}(terms)
    end

    function SumOperator(terms::AbstractVector{PureOperator{S, BR}}) where {S<:Number, BR<:Unsigned}
        return new{S, BR}(terms)
    end
end


bintype(lhs::Type{SumOperator{S, BR}}) where {S, BR} = BR


# === 1/6 Equality ===

function Base.:(==)(lhs::SumOperator{S1, BR}, rhs::SumOperator{S2, BR}) where {S1, S2, BR}
    return (lhs.terms == rhs.terms)
end


# === 2/6 Unary functions ===

Base.real(arg::SumOperator{S, BR}) where {S<:Real, BR} = arg
Base.imag(arg::SumOperator{S, BR}) where {S<:Real, BR} = SumOperator{S, BR}([])
Base.conj(arg::SumOperator{S, BR}) where {S<:Real, BR} = arg

Base.real(arg::SumOperator{Complex{S}, BR}) where {S<:Real, BR} = SumOperator{S, BR}(real.(arg.terms))
Base.imag(arg::SumOperator{Complex{S}, BR}) where {S<:Real, BR} = SumOperator{S, BR}(imag.(arg.terms))
Base.conj(arg::SumOperator{Complex{S}, BR}) where {S<:Real, BR} = SumOperator{Complex{S}, BR}(conj.(arg.terms))

Base.transpose(arg::SumOperator{S, BR}) where {S, BR} = SumOperator{S, BR}(transpose.(arg.terms))
Base.adjoint(arg::SumOperator{S, BR}) where {S, BR} = SumOperator{S, BR}(adjoint.(arg.terms))


# === 3/6 Scalar Operators ===

function Base.:(-)(arg::SumOperator{S, BR}) where {S, BR}
    return SumOperator{S, BR}(-arg.terms)
end


function Base.:(*)(lhs::S1, rhs::SumOperator{S2, BR}) where {S1<:Number, S2<:Number, BR}
    return SumOperator(lhs .* rhs.terms)
end

function Base.:(*)(lhs::SumOperator{S1, BR}, rhs::S2) where {S1<:Number, S2<:Number, BR}
    return SumOperator(lhs.terms .* rhs)
end

function Base.:(\)(lhs::S1, rhs::SumOperator{S2, BR}) where {S1<:Number, S2<:Number, BR}
    return SumOperator(lhs .\ rhs.terms)
end

function Base.:(/)(lhs::SumOperator{S1, BR}, rhs::S2) where {S1<:Number, S2<:Number, BR}
    return SumOperator(lhs.terms ./ rhs)
end

function Base.:(//)(lhs::SumOperator{S1, BR}, rhs::S2) where {S1<:Number, S2<:Number, BR}
    return SumOperator(lhs.terms .// rhs)
end


# === 4/6 Operator Products ===

function Base.:(*)(lhs::SumOperator{S1, BR}, rhs::PureOperator{S2, BR}) where {S1, S2, BR}
    S3 = promote_type(S1, S2)
    terms::Vector{PureOperator{S3, BR}} = filter((x) -> !isa(x, NullOperator), [t * rhs for t in lhs.terms])
    return SumOperator{S3, BR}(terms)
end

function Base.:(*)(lhs::PureOperator{S1, BR}, rhs::SumOperator{S2, BR}) where {S1, S2, BR}
    S3 = promote_type(S1, S2)
    terms::Vector{PureOperator{S3, BR}} = filter((x) -> !isa(x, NullOperator), [lhs * t for t in rhs.terms])
    return SumOperator{S3, BR}(terms)
end

function Base.:(*)(lhs::SumOperator{S1, BR}, rhs::SumOperator{S2, BR}) where {S1, S2, BR}
    S3 = promote_type(S1, S2)
    terms::Vector{PureOperator{S3, BR}} = filter((x) -> !isa(x, NullOperator), [tl * tr for tl in lhs.terms, tr in rhs.terms])
    return SumOperator{S3, BR}(vec(terms))
end


# === 5/6 Operator Sums ===

function Base.:(+)(lhs::PureOperator{S1, BR}, rhs::PureOperator{S2, BR}) where {S1, S2, BR}
    S = promote_type(S1, S2)
    return SumOperator{S, BR}(PureOperator{S, BR}[lhs, rhs])
end

function Base.:(+)(lhs::SumOperator{S1, BR}, rhs::PureOperator{S2, BR}) where {S1, S2, BR}
    S3 = promote_type(S1, S2)
    return SumOperator{S3, BR}(PureOperator{S3, BR}[lhs.terms..., rhs])
end

function Base.:(+)(lhs::PureOperator{S1, BR}, rhs::SumOperator{S2, BR}) where {S1, S2, BR}
    S3 = promote_type(S1, S2)
    return SumOperator{S3, BR}(PureOperator{S3, BR}[lhs, rhs.terms...])
end

function Base.:(+)(lhs::SumOperator{S1, BR}, rhs::SumOperator{S2, BR}) where {S1, S2, BR}
    S3 = promote_type(S1, S2)
    return SumOperator{S3, BR}(PureOperator{S3, BR}[lhs.terms..., rhs.terms...])
end


# === 6/6 Conversion ===

function Base.promote_rule(lhs::Type{SumOperator{S1, B1}}, rhs::Type{SumOperator{S2, B2}}) where {S1, S2, B1, B2}
    S3 = promote_type(S1, S2)
    B3 = promote_type(B1, B2)
    return SumOperator{S3, B3}
end

function Base.convert(type::Type{SumOperator{S1, B1}}, obj::SumOperator{S2, B2}) where {S1, S2, B1, B2}
    return SumOperator{S1, B1}([convert(PureOperator{S1, B1}, t) for t in obj.terms])
end
