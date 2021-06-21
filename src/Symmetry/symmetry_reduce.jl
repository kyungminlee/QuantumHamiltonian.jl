export symmetry_reduce, symmetry_reduce_serial, symmetry_reduce_parallel
export symmetry_reduce!, symmetry_reduce_serial!, symmetry_reduce_parallel!
export symmetry_unreduce


"""
    symmetry_reduce(rhsr, large_vector)

Reduce a large vector into the reduced hilbert space representation.
Simply throw away components that don't fit.
"""
function symmetry_reduce(
    rhsr::ReducedHilbertSpaceRepresentation{HSR, BR, C},
    large_vector::AbstractVector{Si}
) where {HSR, BR, C, Si<:Number}
    symred = Threads.nthreads() == 1 ? symmetry_reduce_serial : symmetry_reduce_parallel
    return symred(rhsr, large_vector)
end


function symmetry_reduce_serial(
    rhsr::ReducedHilbertSpaceRepresentation{HSR, BR, C},
    large_vector::AbstractVector{Si}
) where {HSR, BR, C, Si<:Number}
    if length(large_vector) != dimension(rhsr.parent)
        throw(DimensionMismatch("Dimension of the input vector should match the larger representation"))
    end
    So = promote_type(C, Si)
    small_vector = zeros(So, dimension(rhsr))
    # (i_p | i_r | ampl) indicates : U_(p, r) = ampl
    for (i_p, i_r) in enumerate(rhsr.basis_mapping_index)
        if i_r > 0
            @inbounds ampl = rhsr.basis_mapping_amplitude[i_p]
            # H_r = U† H U
            @inbounds small_vector[i_r] += conj(ampl) * large_vector[i_p]
        end
    end
    return small_vector
end


function symmetry_reduce_parallel(
    rhsr::ReducedHilbertSpaceRepresentation{HSR, BR, C},
    large_vector::AbstractVector{Si}
) where {HSR, BR, C, Si<:Number}
    if length(large_vector) != dimension(rhsr.parent)
        throw(DimensionMismatch("Dimension of the input vector should match the larger representation"))
    end
    So = promote_type(C, Si)
    # (i_p | i_r | ampl) indicates : U_(p, r) = ampl
    n_p = length(rhsr.basis_mapping_index)
    local_small_vectors = [zeros(So, dimension(rhsr)) for tid in 1:Threads.nthreads()]
    Threads.@threads for i_p in 1:n_p
        @inbounds i_r = rhsr.basis_mapping_index[i_p]
        if i_r > 0
            tid = Threads.threadid()
            @inbounds ampl = rhsr.basis_mapping_amplitude[i_p]
            # H_r = U† H U
            @inbounds val = conj(ampl) * large_vector[i_p]
            @inbounds local_small_vectors[tid][i_r] += val
        end
    end
    for tid in 2:Threads.nthreads()
        @inbounds local_small_vectors[1] += local_small_vectors[tid]
    end
    return local_small_vectors[1]
end

"""
    symmetry_reduce!(out, rhsr, largevector)

Adds and not overwrites.
"""
function symmetry_reduce!(
    out::AbstractVector{So},
    rhsr::ReducedHilbertSpaceRepresentation{HSR, BR, C},
    large_vector::AbstractVector{Si}
) where {HSR, BR, C, Si<:Number, So<:Number}
    symred = Threads.nthreads() == 1 ? symmetry_reduce_serial! : symmetry_reduce_parallel!
    return symred(out, rhsr, large_vector)
end


function symmetry_reduce_serial!(
    out::AbstractVector{So},
    rhsr::ReducedHilbertSpaceRepresentation{HSR, BR, C},
    large_vector::AbstractVector{Si}
) where {HSR, BR, C, Si<:Number, So<:Number}
    if length(large_vector) != dimension(rhsr.parent)
        throw(DimensionMismatch("Dimension of the input vector should match the larger representation"))
    end
    if length(out) != dimension(rhsr)
        throw(DimensionMismatch("Dimension of the output vector should match the smaller representation"))
    end
    # (i_p | i_r | ampl) indicates : U_(p, r) = ampl
    for (i_p, i_r) in enumerate(rhsr.basis_mapping_index)
        if i_r > 0
            @inbounds ampl = rhsr.basis_mapping_amplitude[i_p]
            # H_r = U† H U
            @inbounds val = conj(ampl) * large_vector[i_p]
            @inbounds out[i_r] += val
        end
    end
    return out
end


function symmetry_reduce_parallel!(
    out::AbstractVector{So},
    rhsr::ReducedHilbertSpaceRepresentation{HSR, BR, C},
    large_vector::AbstractVector{Si},
) where {HSR, BR, C, Si<:Number, So<:Number}
    if length(large_vector) != dimension(rhsr.parent)
        throw(DimensionMismatch("Dimension of the input vector should match the larger representation"))
    end
    if length(out) != dimension(rhsr)
        throw(DimensionMismatch("Dimension of the output vector should match the smaller representation"))
    end
    # (i_p | i_r | ampl) indicates : U_(p, r) = ampl
    n_p = length(rhsr.basis_mapping_index)
    local_small_vectors = zeros(So, (dimension(rhsr), Threads.nthreads()))
    Threads.@threads for i_p in 1:n_p
        @inbounds i_r = rhsr.basis_mapping_index[i_p]
        if i_r > 0
            tid = Threads.threadid()
            @inbounds ampl = rhsr.basis_mapping_amplitude[i_p]
            # H_r = U† H U
            @inbounds val = conj(ampl) * large_vector[i_p]
            local_small_vectors[i_r, tid] += val
        end
    end
    for tid in 1:Threads.nthreads()
        @inbounds out[:] += view(local_small_vectors, :, tid)
    end
    return out
end


# function symmetry_reduce_parallel!(
#     out::AbstractVector{Complex{Ro}},
#     rhsr::ReducedHilbertSpaceRepresentation{HSR, SIC, BR, C},
#     large_vector::AbstractVector{ComplexF64}
# ) where {HSR, SIC, BR, C, Si<:Number, Ro<:Real}
#     if length(large_vector) != dimension(rhsr.parent)
#         throw(DimensionMismatch("Dimension of the input vector should match the larger representation"))
#     end
#     # So = promote_type(C, Si)

#     # (i_p | i_r | ampl) indicates : U_(p, r) = ampl
#     n_p = length(rhsr.basis_mapping_index)
#     small_vector = Matrix{Threads.Atomic{Ro}}(undef, (2, dimension(rhsr)))
#     for i in 1:dimension(rhsr)
#         small_vector[1, i] = Threads.Atomic{Ro}(0)
#         small_vector[2, i] = Threads.Atomic{Ro}(0)
#     end

#     Threads.@threads for i_p in 1:n_p
#         i_r = rhsr.basis_mapping_index[i_p]
#         if i_r > 0
#             tid = Threads.threadid()
#             ampl = rhsr.basis_mapping_amplitude[i_p]
#             # H_r = U† H U
#             v = conj(ampl) * large_vector[i_p]
#             Threads.atomic_add!(small_vector[1, i_r], real(v))
#             Threads.atomic_add!(small_vector[2, i_r], imag(v))
#         end
#     end
#     for i in 1:dimension(rhsr)
#         out[i] = small_vector[1,i][] + im * small_vector[2, i][]
#     end
#     return out
# end


raw"""
    symmetry_unreduce

```math
\begin{pmatrix} l_1 \\ l_2 \\ l_3 \\ \vdots \\ l_n \end{pmatrix}
=
\begin{pmatrix}
. & \cdots & . \\
. & \cdots & . \\
. & \cdots & . \\
  & \dots & \\
. & \cdots &
\end{pmatrix}
\begin{pmatrix} s_1 \\ \vdots \\ s_m \end{pmatrix}
```
"""
function symmetry_unreduce(
    rhsr::ReducedHilbertSpaceRepresentation{HSR, BR, C},
    small_vector::AbstractVector{Si}
) where {HSR, BR, C, Si<:Number}
    if length(small_vector) != dimension(rhsr)
        throw(DimensionMismatch("Dimension of the input vector should match the reduced representation"))
    end
    So = promote_type(C, Si)
    large_vector = zeros(So, dimension(rhsr.parent))
    for (i_p, i_r) in enumerate(rhsr.basis_mapping_index)
        if i_r > 0
            ampl = rhsr.basis_mapping_amplitude[i_p]
            # H_r = U† H U
            large_vector[i_p] += ampl * small_vector[i_r]
        end
    end
    return large_vector
end
