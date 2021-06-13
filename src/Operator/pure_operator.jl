export PureOperator
export pure_operator
export bintype

using LinearAlgebra


"""
    PureOperator{Scalar, BR}

Represents an operator ``α (P₁ ⊗ P₂ ⊗ … ⊗ Pₙ)`` where
``Pᵢ`` is either identity (when bitmask is set to zero),
or projection ``|rᵢ⟩⟨cᵢ|`` (when bitmask is set to one).

See also: [`pure_operator`](@ref)

# Fields
```
bitmask   :: BR
bitrow    :: BR
bitcol    :: BR
amplitude :: Scalar
```
"""
struct PureOperator{Scalar<:Number, BR<:Unsigned} <:AbstractOperator{Scalar}
    bitmask::BR
    bitrow::BR
    bitcol::BR
    amplitude::Scalar

    function PureOperator{S, BR}(bitmask::Unsigned, bitrow::Unsigned, bitcol::Unsigned, amplitude::Number) where {S, BR}
        if (~bitmask) & bitrow != zero(BR)
            throw(ArgumentError("every bit of bitrow not in bitmask should be set to zero"))
        elseif (~bitmask) & bitcol != zero(BR)
            throw(ArgumentError("every bit of bitcol not in bitmask should be set to zero"))
        end
        return new{S, BR}(bitmask, bitrow, bitcol, amplitude)
    end

    function PureOperator(bitmask::BR, bitrow::BR, bitcol::BR, amplitude::S) where {S<:Number, BR<:Unsigned}
        if (~bitmask) & bitrow != zero(BR)
            throw(ArgumentError("every bit of bitrow not in bitmask should be set to zero"))
        elseif (~bitmask) & bitcol != zero(BR)
            throw(ArgumentError("every bit of bitcol not in bitmask should be set to zero"))
        end
        return new{S, BR}(bitmask, bitrow, bitcol, amplitude)
    end

    function PureOperator{S, BR}(s::UniformScaling{S}) where {S, BR}
        return new{S, BR}(zero(BR), zero(BR), zero(BR), s.λ)
    end
end


bintype(::Type{PureOperator{S, BR}}) where {S, BR} = BR


Base.zero(::Type{PureOperator{S, BR}}) where {S, BR} = PureOperator{S, BR}(zero(BR), zero(BR), zero(BR), zero(S))
Base.zero(::PureOperator{S, BR}) where {S, BR} = PureOperator{S, BR}(zero(BR), zero(BR), zero(BR), zero(S))
Base.one(::Type{PureOperator{S, BR}}) where {S, BR} = PureOperator{S, BR}(zero(BR), zero(BR), zero(BR), one(S))
Base.one(::PureOperator{S, BR}) where {S, BR} = PureOperator{S, BR}(zero(BR), zero(BR), zero(BR), one(S))

Base.iszero(x::PureOperator) = iszero(x.amplitude)


# === 1/6 (In)equality ===

function Base.:(==)(lhs::PureOperator{S1, BR}, rhs::PureOperator{S2, BR}) where {S1, S2, BR}
    return (
        (lhs.bitmask == rhs.bitmask) &&
        (lhs.bitrow == rhs.bitrow) &&
        (lhs.bitcol == rhs.bitcol) &&
        (lhs.amplitude == rhs.amplitude)
    )
end


function Base.isless(lhs::PureOperator{S1, BR}, rhs::PureOperator{S2, BR}) where {S1, S2, BR}
    lhs.bitmask < rhs.bitmask && return true
    lhs.bitmask > rhs.bitmask && return false
    lhs.bitrow < rhs.bitrow && return true
    lhs.bitrow > rhs.bitrow && return false
    lhs.bitcol < rhs.bitcol && return true
    lhs.bitcol > rhs.bitcol && return false
    real(lhs.amplitude) < real(rhs.amplitude) && return true
    real(lhs.amplitude) > real(rhs.amplitude) && return false
    imag(lhs.amplitude) < imag(rhs.amplitude) && return true
    imag(lhs.amplitude) > imag(rhs.amplitude) && return false
    return false
end


# === 2/6 Unary functions ===

Base.real(arg::PureOperator{R, BR}) where {R<:Real, BR} = arg
Base.imag(arg::PureOperator{R, BR}) where {R<:Real, BR} = PureOperator{R, BR}(arg.bitmask, arg.bitrow, arg.bitcol, zero(R))
Base.conj(arg::PureOperator{R, BR}) where {R<:Real, BR} = arg

Base.real(arg::PureOperator{Complex{R}, BR}) where {R<:Real, BR} = PureOperator{R, BR}(arg.bitmask, arg.bitrow, arg.bitcol, real(arg.amplitude))
Base.imag(arg::PureOperator{Complex{R}, BR}) where {R<:Real, BR} = PureOperator{R, BR}(arg.bitmask, arg.bitrow, arg.bitcol, imag(arg.amplitude))
Base.conj(arg::PureOperator{Complex{R}, BR}) where {R<:Real, BR} = PureOperator{Complex{R}, BR}(arg.bitmask, arg.bitrow, arg.bitcol, conj(arg.amplitude))

Base.transpose(arg::PureOperator{S, BR}) where {S, BR} = PureOperator{S, BR}(arg.bitmask, arg.bitcol, arg.bitrow, arg.amplitude)
Base.adjoint(arg::PureOperator{S, BR}) where {S, BR} = PureOperator{S, BR}(arg.bitmask, arg.bitcol, arg.bitrow, conj(arg.amplitude))


# === 3/6 Scalar Operators ===

function Base.:(-)(op::PureOperator{S, BR}) where {S, BR}
    return PureOperator{S, BR}(op.bitmask, op.bitrow, op.bitcol, -op.amplitude)
end

function Base.:(*)(lhs::S1, rhs::PureOperator{S2, BR}) where {S1<:Number, S2<:Number, BR}
    return PureOperator(rhs.bitmask, rhs.bitrow, rhs.bitcol, lhs * rhs.amplitude)
end

function Base.:(*)(lhs::PureOperator{S1, BR}, rhs::S2) where {S1<:Number, S2<:Number, BR}
    return PureOperator(lhs.bitmask, lhs.bitrow, lhs.bitcol, lhs.amplitude * rhs)
end

function Base.:(\)(lhs::S1, rhs::PureOperator{S2, BR}) where {S1<:Number, S2<:Number, BR}
    return PureOperator(rhs.bitmask, rhs.bitrow, rhs.bitcol, lhs \ rhs.amplitude)
end

function Base.:(/)(lhs::PureOperator{S1, BR}, rhs::S2) where {S1<:Number, S2<:Number, BR}
    return PureOperator(lhs.bitmask, lhs.bitrow, lhs.bitcol, lhs.amplitude / rhs)
end

function Base.:(//)(lhs::PureOperator{S1, BR}, rhs::S2) where {S1<:Number, S2<:Number, BR}
    return PureOperator(lhs.bitmask, lhs.bitrow, lhs.bitcol, lhs.amplitude // rhs)
end


# === 4/6 Operator Products ===

function Base.:(*)(
    lhs::PureOperator{S1, BR},
    rhs::PureOperator{S2, BR}
)::AbstractOperator where {S1<:Number, S2<:Number, BR}
    onlylhs_bitmask   =   lhs.bitmask  & (~rhs.bitmask)
    onlyrhs_bitmask   = (~lhs.bitmask) &   rhs.bitmask
    intersect_bitmask =   lhs.bitmask  &   rhs.bitmask
    union_bitmask     =   lhs.bitmask  |   rhs.bitmask

    if (lhs.bitcol & intersect_bitmask) != (rhs.bitrow & intersect_bitmask)
        return NullOperator()
    else
        new_bitmask = union_bitmask
        new_bitrow = lhs.bitrow | (rhs.bitrow & onlyrhs_bitmask)
        new_bitcol = rhs.bitcol | (lhs.bitcol & onlylhs_bitmask)
        new_amplitude = lhs.amplitude * rhs.amplitude
        return PureOperator(new_bitmask, new_bitrow, new_bitcol, new_amplitude)
    end
end

# === 5/6 Operator Sums ===
# implemented in sum_operator.jl

# === 6/6 Conversion ===

function Base.promote_rule(::Type{PureOperator{S1, B1}}, ::Type{PureOperator{S2, B2}}) where {S1, S2, B1, B2}
    S3 = promote_type(S1, S2)
    B3 = promote_type(B1, B2)
    return PureOperator{S3, B3}
end

function Base.convert(::Type{PureOperator{S1, B1}}, obj::PureOperator{S2, B2}) where {S1, S2, B1, B2}
    return PureOperator{S1, B1}(
        convert(B1, obj.bitmask),
        convert(B1, obj.bitrow),
        convert(B1, obj.bitcol),
        convert(S1, obj.amplitude)
    )
end

function Base.promote_rule(::Type{PureOperator{S1, B1}}, ::UniformScaling{S2}) where {S1, S2, B1}
    S3 = promote_type(S1, S2)
    return PureOperator{S3, B1}
end

function Base.convert(::Type{PureOperator{S1, B1}}, obj::UniformScaling{S2}) where {S1, S2, B1}
    return PureOperator{S1, B1}(zero(B1), zero(B1), zero(B1), obj.λ)
end

function Base.convert(::Type{AbstractOperator{S1}}, obj::PureOperator{S2, B}) where {S1, S2, B}
    return convert(PureOperator{S1, B}, obj)
end

"""
    pure_operator(hilbert_space, isite, istate_row, istate_col, amplitude=1, binary_type=UInt)

Creates a pure operator where projection is at one of the sites.
"""
function pure_operator(
    hilbert_space::HilbertSpace,
    isite::Integer,
    istate_row::Integer,
    istate_col::Integer,
    amplitude::S=1,
    binary_type::Type{BR}=UInt,
) where {S<:Number, BR<:Unsigned}
    @boundscheck let
        site = hilbert_space.sites[isite]
        state_row = site.states[istate_row]
        state_col = site.states[istate_col]
    end
    bm = get_bitmask(hilbert_space, isite, BR)
    br = BR(istate_row - 1) << hilbert_space.bitoffsets[isite]
    bc = BR(istate_col - 1) << hilbert_space.bitoffsets[isite]
    return PureOperator(bm, br, bc, amplitude)
end
