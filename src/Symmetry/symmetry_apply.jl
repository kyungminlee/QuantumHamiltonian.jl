export symmetry_apply
export isinvariant

using GroupTools
# import LatticeTools.AbstractSpaceSymmetryOperationEmbedding
import LatticeTools.SymmetryEmbedding

## AbstractSpaceSymmetryOperationEmbedding

### HilbertSpaceSector
function symmetry_apply(
    hss::HilbertSpaceSector,
    symop::AbstractSymmetryOperation,
    args...;
    kwargs...
)
    return symmetry_apply(hss.parent, symop, args...; kwargs...)
end

function isinvariant(
    hss::HilbertSpaceSector,
    symop::AbstractSymmetryOperation,
    args...;
    kwargs...
)
    return isinvariant(hss.parent, symop, args...; kwargs...)
end

function isinvariant(
    hss::HilbertSpaceSector,
    symbed::SymmetryEmbedding,
    args...;
    kwargs...
)
    return isinvariant(hss.parent, symbed, args...; kwargs...)
end

function isinvariant(
    hss::HilbertSpaceSector,
    symbed::SymmorphicSymmetryEmbedding,
    args...;
    kwargs...
)
    return isinvariant(hss.parent, symbed, args...; kwargs...)
end


### generic symmetry operations for NullOperator and SumOperator
function symmetry_apply(
    hs::AbstractHilbertSpace,
    symop::AbstractSymmetryOperation,
    op::NullOperator
)
    return op
end

function symmetry_apply(
    hs::AbstractHilbertSpace,
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

function symmetry_apply(
    hs::AbstractHilbertSpace,
    dop::ProductOperation,
    bitrep::Unsigned,
    amplitude::Number=one(Int)
)
    for op in reverse(dop.operations)  # (ABC)(ψ) = A(B(C(ψ)))
        bitrep, amplitude = symmetry_apply(hs, op, bitrep, amplitude)
    end
    return (bitrep, amplitude)
end

function symmetry_apply(
    hs::AbstractHilbertSpace,
    p::SemidirectProductOperation,
    bitrep::Unsigned,
    amplitude::Number=one(Int)
)
    (bitrep, amplitude) = symmetry_apply(hs, p.normal, bitrep, amplitude)
    (bitrep, amplitude) = symmetry_apply(hs, p.rest, bitrep, amplitude)
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


function symmetry_apply(
    hs::AbstractHilbertSpace,
    perm::LocalGeneralizedPermutation{A},
    bitrep::BR,
    amplitude::S=one(Int)
) where {A, BR<:Unsigned, S<:Number}
    @boundscheck let
        if length(hs.sites) != length(perm.operations)
            throw(ArgumentError("number of sites should match"))
        end
        for (isite, (site, op)) in enumerate(zip(hs.sites, perm.operations))
            if dimension(site) != length(op.map)
                throw(ArgumentError("dimension mismatch at site $isite: $(dimension(site)) != $(length(op.map))"))
            end
        end
    end
    indexarray = extract(hs, bitrep)
    new_indexarray = CartesianIndex([
        let 
            istate_new, amplitude = op(istate, amplitude)
            istate_new
        end
        for (istate, op) in zip(indexarray.I, perm.operations)
    ]...)
    return (compress(hs, new_indexarray), amplitude)
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
    hs::HilbertSpace,
    symop::AbstractSymmetryOperation,
    op::AbstractOperator
)
    return simplify(op - symmetry_apply(hs, symop, op)) == NullOperator()
end

function isinvariant(
    hs::HilbertSpace,
    symbed::SymmetryEmbedding,
    op::AbstractOperator
)
    return all(isinvariant(hs, g, op) for g in generator_elements(symbed))
end

function isinvariant(
    hs::HilbertSpace,
    symbed::SymmorphicSymmetryEmbedding,
    op::AbstractOperator
)
    return (
        all(isinvariant(hs, g, op) for g in generator_elements(symbed.normal)) &&
        all(isinvariant(hs, g, op) for g in generator_elements(symbed.rest))
    )
end
