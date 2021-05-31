export LocalGeneralizedPermutation

struct LocalGeneralizedPermutation{A} <: AbstractSymmetryOperation
    operations::Vector{GeneralizedPermutation{A}}
    LocalGeneralizedPermutation(operations::AbstractVector{GeneralizedPermutation{A}}) where {A} = new{A}(operations)
end

function Base.:(*)(lhs::LocalGeneralizedPermutation, rhs::LocalGeneralizedPermutation)
    return LocalGeneralizedPermutation(lhs.operations .* rhs.operations)
end

function LatticeTools.isidentity(arg::LocalGeneralizedPermutation)
    return all(LatticeTools.isidentity, arg.operations)
end

Base.inv(arg::LocalGeneralizedPermutation) = LocalGeneralizedPermutation(inv.(arg.operations))

function Base.:(*)(lhs::LocalGeneralizedPermutation, rhs::SitePermutation)
    return DirectProductOperation(lhs, rhs)
end

# Simpler version
function Base.:(*)(lhs::SitePermutation, rhs::LocalGeneralizedPermutation)
    return DirectProductOperation(lhs, rhs)
end

# # SitePermutation P: ψ[i] ↦ ψ[P(i)] = ψ'[i]
# # LocalUnitary    U: ψ[i] ↦ U[i] * ψ[i] = ψ'[i]
# # U * P : ψ''[i] = U[i] * ψ'[i] = U[i] * ψ[P(i)]
# # P * U : ψ''[i] = ψ'[P(i)] = U[P(i)] * ψ[P(i)]
# function Base.:(*)(lhs::SitePermutation, rhs::LocalGeneralizedPermutation)
#     lhs_new = lhs
#     rhs.operations[lhs_new]
#     return DirectProductOperation(lhs_new, rhs_new)
# end

function Base.:(==)(lhs::LocalGeneralizedPermutation, rhs::LocalGeneralizedPermutation)
    return all(lhs.operations .== rhs.operations)
end

