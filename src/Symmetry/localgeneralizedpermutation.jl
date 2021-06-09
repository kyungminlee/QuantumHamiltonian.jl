export LocalGeneralizedPermutation

using GroupTools

struct LocalGeneralizedPermutation{A} <: GroupTools.AbstractSymmetryOperation
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

Base.hash(arg::L, h::UInt) where {L<:LocalGeneralizedPermutation} = hash(L, hash(arg.operations, h))



# NOTE! NOT DIRECT PRODUCT!!!
# function Base.:(*)(lhs::LocalGeneralizedPermutation, rhs::SitePermutation)
#     return DirectProductOperation(lhs, rhs)
# end

# # Simpler version
# function Base.:(*)(lhs::SitePermutation, rhs::LocalGeneralizedPermutation)
#     return DirectProductOperation(lhs, rhs)
# end

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


# function Base.:(*)(
#     x::ProductOperation{<:Tuple{<:SitePermutation, <:LocalGeneralizedPermutation}},
#     y::ProductOperation{<:Tuple{<:SitePermutation, <:LocalGeneralizedPermutation}}
# )
#     px, ux = x.operations
#     py, uy = y.operations

#     n = length(px.permutation.map)
#     if n != length(py.permutation.map)
#         throw(ArgumentError("the two operations should have the same number of sites"))
#     end

#     pz = px * py
#     # inv_py = inv(py)
#     uz = LocalGeneralizedPermutation([
#         ux.operations[ py(i) ] * uy.operations[i]
#         for i in 1:n
#     ])
#     return ProductOperation(pz, uz)
# end



function (permutation::SitePermutation)(op::LocalGeneralizedPermutation)
    out = similar(op.operations)
    for (i, x) in enumerate(op.operations)
        out[permutation(i)] = x
    end
    return LocalGeneralizedPermutation(out)
end




# function Base.:(*)(
#     x::ProductOperation{<:Tuple{<:SitePermutation, <:LocalGeneralizedPermutation}},
#     y::ProductOperation{<:Tuple{<:SitePermutation, <:LocalGeneralizedPermutation}}
# )
#     px, ux = x.operations
#     py, uy = y.operations

#     n = length(px.permutation.map)
#     if n != length(py.permutation.map)
#         throw(ArgumentError("the two operations should have the same number of sites"))
#     end
#     pz = px * py
#     uz = py(ux) * uy
#     return ProductOperation(pz, uz)
# end

