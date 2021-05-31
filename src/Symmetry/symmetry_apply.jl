export symmetry_apply
export isinvariant

using GroupTools
# import LatticeTools.AbstractSpaceSymmetryOperationEmbedding
import LatticeTools.SymmetryEmbedding

## AbstractSpaceSymmetryOperationEmbedding

### HilbertSpaceSector
function symmetry_apply(
    hss::HilbertSpaceSector{QN},
    symop::AbstractSymmetryOperation,
    args...;
    kwargs...
) where {QN}
    return symmetry_apply(hss.parent, symop, args...; kwargs...)
end

function isinvariant(
    hss::HilbertSpaceSector{QN},
    symop::AbstractSymmetryOperation,
    args...;
    kwargs...
) where {QN}
    return isinvariant(hss.parent, symop, args...; kwargs...)
end

function isinvariant(
    hss::HilbertSpaceSector{QN},
    symbed::SymmetryEmbedding,
    args...;
    kwargs...
) where {QN}
    return isinvariant(hss.parent, symbed, args...; kwargs...)
end

function isinvariant(
    hss::HilbertSpaceSector{QN},
    symbed::SymmorphicSymmetryEmbedding,
    args...;
    kwargs...
) where {QN}
    return isinvariant(hss.parent, symbed, args...; kwargs...)
end


### generic symmetry operations for NullOperator and SumOperator
function symmetry_apply(
    hs::HilbertSpace{QN},
    symop::AbstractSymmetryOperation,
    op::NullOperator
) where {QN}
    return op
end

function symmetry_apply(
    hs::HilbertSpace,
    symop::AbstractSymmetryOperation,
    op::SumOperator
)
    terms = collect(symmetry_apply(hs, symop, t) for t in op.terms)
    return SumOperator(terms)
end

function symmetry_apply(
    hs::AbstractHilbertSpace,
    dop::DirectProductOperation,
    bitrep::BR,
    amplitude::S=one(Int)
) where {BR<:Unsigned, S<:Number}
    # assumption is that the operations commute.
    for op in reverse(dop.operations)  # (ABC)(ψ) = A(B(C(ψ)))
        bitrep, amplitude = symmetry_apply(hs, op, bitrep, amplitude)
    end
    return (bitrep, amplitude)
end


## Permutation
### Binary Representation
function symmetry_apply(
    hs::HilbertSpace{QN},
    permutation::SitePermutation,
    bitrep::BR,
    amplitude::S=one(Int)
    # TODO: bitmask::BR,
) where {QN, BR<:Unsigned, S<:Number}
    out = zero(BR)
    for (i, j) in enumerate(permutation.permutation.map)
        out |= ( (bitrep >> hs.bitoffsets[i]) & make_bitmask(hs.bitwidths[i]) ) << hs.bitoffsets[j]
    end
    return (out, amplitude)
end


### Operator
function symmetry_apply(
    hs::HilbertSpace{QN},
    permutation::SitePermutation,
    op::PureOperator{S, BR}
) where {QN, S<:Number, BR<:Unsigned}
    bm, _ = symmetry_apply(hs, permutation, op.bitmask)
    br, _ = symmetry_apply(hs, permutation, op.bitrow)
    bc, _ = symmetry_apply(hs, permutation, op.bitcol)
    am = op.amplitude
    return PureOperator{S, BR}(bm, br, bc, am)
end


## isinvariant
function isinvariant(
    hs::HilbertSpace{QN},
    symop::AbstractSymmetryOperation,
    op::AbstractOperator
) where {QN}
    return simplify(op - symmetry_apply(hs, symop, op)) == NullOperator()
end

function isinvariant(
    hs::HilbertSpace{QN},
    symbed::SymmetryEmbedding,
    op::AbstractOperator
) where {QN}
    return all(isinvariant(hs, g, op) for g in generator_elements(symbed))
end

function isinvariant(
    hs::HilbertSpace{QN},
    symbed::SymmorphicSymmetryEmbedding,
    op::AbstractOperator
) where {QN}
    return (
        all(isinvariant(hs, g, op) for g in generator_elements(symbed.normal)) &&
        all(isinvariant(hs, g, op) for g in generator_elements(symbed.rest))
    )
end
