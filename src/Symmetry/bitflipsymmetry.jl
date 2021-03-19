# Can be used as spin flips for spin-half systems

export GlobalBitFlip
export symmetry_apply

import LatticeTools

"""
    GlobalBitFlip

Global bit flip operation.
For spin half systems, this amounts to ℤ₂ spin flip.
"""
struct GlobalBitFlip <: AbstractSymmetryOperation
    value::Bool
    GlobalBitFlip() = new(false)
    GlobalBitFlip(value::Bool) = new(value)
end

Base.:(*)(lhs::GlobalBitFlip, rhs::GlobalBitFlip) = GlobalBitFlip(lhs.value ⊻ rhs.value)
LatticeTools.isidentity(arg::GlobalBitFlip) = !arg.value

# How GlobalBitFlip combines with SitePermutation
Base.:(*)(lhs::GlobalBitFlip, rhs::SitePermutation) = DirectProductOperation(lhs, rhs)
Base.:(*)(lhs::SitePermutation, rhs::GlobalBitFlip) = DirectProductOperation(lhs, rhs)

Base.:(==)(lhs::GlobalBitFlip, rhs::GlobalBitFlip) = lhs.value == rhs.value

function symmetry_apply(
    hs::HilbertSpace,
    op::GlobalBitFlip,
    bitrep::BR
) where {BR<:Unsigned}
    if op.value
        @boundscheck (sizeof(BR)*8 >= bitwidth(hs)) || throw(ArgumentError("binary representation too small"))
        mask = get_bitmask(hs, BR)
        return (mask & (~bitrep), 1)
    else
        return (bitrep, 1)
    end
end
