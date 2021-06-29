
export AbstractOperatorRepresentation
export spacetype, operatortype
export bintype
export get_space
export dimension, bitwidth
export get_row, get_column
export sparse_serial, sparse_parallel

export apply!, apply_serial!, apply_parallel!


"""
    AbstractOperatorRepresentation{S}
"""
abstract type AbstractOperatorRepresentation{S} <: AbstractMatrix{S} end


## typetraits

Base.valtype(::Type{<:AbstractOperatorRepresentation{T}}) where T = T
Base.valtype(::AbstractOperatorRepresentation{T}) where T = T

scalartype(::Type{<:AbstractOperatorRepresentation{T}}) where T = T
scalartype(::AbstractOperatorRepresentation{T}) where T = T

bintype(lhs::Type{<:AbstractOperatorRepresentation}) = bintype(spacetype(lhs))
bintype(lhs::AbstractOperatorRepresentation) = bintype(typeof(lhs))


# a subclass of AbstractOperatorRepresentation should implement
# spacetype, operatortype, and get_space.
spacetype(lhs::AbstractOperatorRepresentation{T}) where T = spacetype(typeof(lhs))::Type{<:AbstractHilbertSpaceRepresentation}
operatortype(lhs::AbstractOperatorRepresentation{T}) where T = operatortype(typeof(lhs))::Type{<:AbstractOperator}


dimension(lhs::AbstractOperatorRepresentation{S}) where S = dimension(get_space(lhs))
bitwidth(lhs::AbstractOperatorRepresentation{S}) where S = bitwidth(get_space(lhs))


function Base.size(arg::AbstractOperatorRepresentation{T})::Tuple{Int, Int} where T
    dim = dimension(get_space(arg))
    return (dim, dim)
end

Base.size(arg::AbstractOperatorRepresentation{T}, i::Integer) where T = Base.size(arg)[i]


function Base.:(==)(
    lhs::AbstractOperatorRepresentation{T1},
    rhs::AbstractOperatorRepresentation{T2}
) where {T1, T2}
    return ((get_space(lhs) == get_space(rhs)) && (lhs.operator == rhs.operator))
end


for uniop in [:+, :-]
    @eval begin
        function Base.$uniop(lhs::AbstractOperatorRepresentation{T}) where T
            return represent(get_space(lhs), ($uniop)(lhs.operator))
        end
    end
end

for binop in [:+, :-, :*]
    @eval begin
        function Base.$binop(
            lhs::AbstractOperatorRepresentation{T1},
            rhs::AbstractOperatorRepresentation{T2}
        ) where {T1, T2}
            @boundscheck if (get_space(lhs) != get_space(rhs))
                throw(ArgumentError("The two OperatorRepresentation's do not have the same HilbertSpaceRepresentation"))
            end
            return represent(get_space(lhs), simplify(($binop)(lhs.operator, rhs.operator)))
        end
    end
end

function Base.:(*)(lhs::AbstractOperatorRepresentation, rhs::Number)
    return represent(get_space(lhs), simplify(lhs.operator * rhs))
end

function Base.:(*)(lhs::Number, rhs::AbstractOperatorRepresentation)
    return represent(get_space(rhs), simplify(lhs * rhs.operator))
end

function Base.:(/)(lhs::AbstractOperatorRepresentation, rhs::Number)
    return represent(get_space(lhs), simplify(lhs.operator / rhs))
end

function Base.:(\)(lhs::Number, rhs::AbstractOperatorRepresentation)
    return represent(get_space(rhs), simplify(lhs \ rhs.operator))
end

function Base.:(//)(lhs::AbstractOperatorRepresentation, rhs::Number)
    return represent(get_space(lhs), simplify(lhs.operator // rhs))
end


function simplify(arg::AbstractOperatorRepresentation)
    return represent(get_space(arg), simplify(arg.operator))
end


function LinearAlgebra.ishermitian(arg::AbstractOperatorRepresentation)
    return LinearAlgebra.ishermitian(arg.operator)
end


function LinearAlgebra.mul!(
    out::AbstractVector{S1},
    opr::AbstractOperatorRepresentation,
    state::AbstractVector
) where {S1<:Number}
    fill!(out, zero(S1))
    apply!(out, opr, state)
    return out
end


"""
    C(i, k) = A(i, j) * B(j, k) 
    C(:, k) = A(:, :) * B(:, k)
"""
function LinearAlgebra.mul!(
    C::AbstractMatrix{S1},
    A::AbstractOperatorRepresentation,
    B::AbstractMatrix
) where {S1<:Number}
    fill!(C, zero(S1))
    Threads.@threads for k in 1:size(B, 2)
        apply!(view(C, :, k), A, view(B, :, k))
    end
    return C
end


function Base.:(*)(A::AbstractOperatorRepresentation{S}, B::AbstractMatrix{T}) where {S, T}
    size(A, 2) == size(B, 1) || throw(DimensionMismatch("A has size $(size(A)) and B has size $(size(B))"))
    U = promote_type(S, T)
    C = Matrix{U}(undef, (size(A, 1), size(B, 2)))
    LinearAlgebra.mul!(C, A, B)
    return C
end


function Base.:(*)(A::AbstractOperatorRepresentation{S}, B::AbstractVector{T}) where {S, T}
    size(A, 2) == length(B) || throw(DimensionMismatch("A has size $(size(A)) and B has size $(size(B))"))
    U = promote_type(S, T)
    C = Vector{U}(undef, size(A, 1))
    LinearAlgebra.mul!(C, A, B)
    return C
end

# C(i, k) = A(i, j) * B(j, k) 
# C(:, k) = A(:, :) * B(:, k)
# C(i, :) = A(i, :) * B(:, :)
function Base.:(*)(A::AbstractMatrix{T}, B::AbstractOperatorRepresentation{S}) where {S, T}
    size(A, 2) == size(B, 1) || throw(DimensionMismatch("A has size $(size(A)) and B has size $(size(B))"))
    U = promote_type(S, T)
    C = zeros(U, (size(A, 1), size(B, 2)))
    Threads.@threads for i in 1:size(A, 1)
        apply!(view(C, i, :), view(A, i, :), B)
    end
    return C
end


function Base.:(*)(A::AbstractVector{T}, B::AbstractOperatorRepresentation{S}) where {S, T}
    length(A) == size(B, 1) || throw(DimensionMismatch("A has size $(size(A)) and B has size $(size(B))"))
    U = promote_type(S, T)
    C = zeros(U, size(B, 2))
    apply!(C, A, B)
    return C
end


# function Base.:(*)(A::AbstractOperatorRepresentation{S}, B::AbstractArray{T}) where {S, T}
#     size(A, 2) == size(B, 1) || throw(DimensionMismatch("A has size $(size(A)) and B has size $(size(B))"))
#     U = promote_type(S, T)
#     iter_size = size(B)[2:end]
#     C = zeros(U, (size(A, 1), iter_size...))
#     for i3 in CartesianIndices(iter_size)
#         apply!(view(C, :, i3.I...), A, view(B, :, i3.I...))
#     end
#     return C
# end


function Base.Matrix(opr::AbstractOperatorRepresentation{S}) where S
    m, n = size(opr)
    out = zeros(S, (m, n))
    Threads.@threads for icol in 1:n
        for (irow, ampl) in get_column_iterator(opr, icol)
            if 1 <= irow <= m
                out[irow, icol] += ampl
            end
        end
    end
    return out
end


import SparseArrays
function SparseArrays.sparse(
    opr::AbstractOperatorRepresentation{S};
    tol::Real=Base.rtoldefault(real(S))
) where S
    sp = Threads.nthreads() == 1 ? sparse_serial : sparse_parallel
    return sp(opr; tol=tol)
end


function sparse_serial(
    opr::AbstractOperatorRepresentation{S};
    tol::Real=Base.rtoldefault(real(S))
) where S
    m, n = size(opr)
    colptr = zeros(Int, n+1)
    rowval = Int[]
    nzval = S[]

    colptr[1] = 1
    for icol in 1:n
        colvec = Dict{Int, S}()
        for (irow, ampl) in get_column_iterator(opr, icol)
            if 1 <= irow <= m
                colvec[irow] = get(colvec, irow, zero(S)) + ampl
            end
        end
        choptol!(colvec, tol)
        colptr[icol+1] = colptr[icol] + length(colvec)
        sorted_items = sort(collect(colvec), by = item -> item[1])
        append!(rowval, irow for (irow, ampl) in sorted_items)
        append!(nzval, ampl for (irow, ampl) in sorted_items)
    end
    return SparseMatrixCSC{S, Int}(m, n, colptr, rowval, nzval)
end


function sparse_parallel(
    opr::AbstractOperatorRepresentation{S};
    tol::Real=Base.rtoldefault(real(S))
) where S
    m, n = size(opr)

    colsize = zeros(Int, n)

    nblocks = Threads.nthreads()
    block_ranges = splitrange(1:n, nblocks) # guaranteed to be ordered

    local_rowval = Vector{Int}[Int[] for i in 1:nblocks]
    local_nzval = Vector{S}[S[] for i in 1:nblocks]

    Threads.@threads for iblock in 1:nblocks
        subrange = block_ranges[iblock]
        for icol in subrange
            colvec = Dict{Int, S}()
            for (irow, ampl) in get_column_iterator(opr, icol)
                if 1 <= irow <= m
                    colvec[irow] = get(colvec, irow, zero(S)) + ampl
                end
            end
            choptol!(colvec, tol)
            sorted_items = sort(collect(colvec), by=item->item[1])
            colsize[icol] = length(sorted_items)
            append!(local_rowval[iblock], irow for (irow, ampl) in sorted_items)
            append!(local_nzval[iblock], ampl for (irow, ampl) in sorted_items)
        end
    end

    colptr = cumsum(Int[1, colsize...])
    rowval = vcat(local_rowval...)
    nzval = vcat(local_nzval...)
    return SparseMatrixCSC{S, Int}(m, n, colptr, rowval, nzval)
end


function get_row(opr::AbstractOperatorRepresentation{S}, irow::Integer) where S
    Z = zero(S)
    dim = size(opr, 2)
    items = Dict{Int, S}()
    for (icol, val) in get_row_iterator(opr, irow)
        if 1 <= icol <= dim
            items[icol] = get(items, icol, Z) + val
        end
    end
    choptol!(items, Base.rtoldefault(Float64))
    return sparsevec(items, dim)
end


function get_column(opr::AbstractOperatorRepresentation{S}, icol::Integer) where S
    Z = zero(S)
    dim = size(opr, 1)
    items = Dict{Int, S}()
    for (irow, val) in get_column_iterator(opr, icol)
        if 1 <= irow <= dim
            items[irow] = get(items, irow, Z) + val
        end
    end
    choptol!(items, Base.rtoldefault(Float64))
    return sparsevec(items, dim)
end


@inline function Base.getindex(oprep::AbstractOperatorRepresentation, ::Colon, icol::Integer)
    return get_column(oprep, icol)
end

@inline function Base.getindex(oprep::AbstractOperatorRepresentation, irow::Integer, ::Colon)
    return get_row(oprep, irow)
end

@inline function Base.getindex(oprep::AbstractOperatorRepresentation, ::Colon, ::Colon)
    return sparse(oprep)
end

@inline function Base.getindex(oprep::AbstractOperatorRepresentation, irow::Integer, icol::Integer)
    return get_element(oprep, irow, icol)
end


"""
    apply!(out, opr, state)

Perform `out += opr * state`. Apply the operator representation `opr` to the
column vector `state` and *add* it to the column vector `out`.
Return sum of errors and sum of error-squared.
Call [`apply_serial!`](@ref) if `Threads.nthreads() == 1`, and [`apply_parallel!`](@ref) otherwise.
"""
function apply!(
    out::AbstractVector,
    opr::AbstractOperatorRepresentation,
    state::AbstractVector
)
    a! = Threads.nthreads() == 1 ? apply_serial! : apply_parallel!
    return a!(out, opr, state)
end


"""
    apply!(out, opr, state)

Perform `out += opr * state`. Apply the operator representation `opr` to the
row vector `state` and *add* it to the row vector `out`.
Return sum of errors and sum of error-squared.
Call [`apply_serial!`](@ref) if `Threads.nthreads() == 1`, and [`apply_parallel!`](@ref) otherwise.
"""
function apply!(
    out::AbstractVector,
    state::AbstractVector,
    opr::AbstractOperatorRepresentation
)
    a! = Threads.nthreads() == 1 ? apply_serial! : apply_parallel!
    return a!(out, state, opr)
end


"""
    apply_serial!(out, opr, state)

Perform `out += opr * state`. Apply the operator representation `opr` to the
column vector `state` and *add* it to the column vector `out`.
Return sum of errors and sum of error-squared.
Single-threaded version.
"""
function apply_serial!(
    out::AbstractVector,
    opr::AbstractOperatorRepresentation{S},
    state::AbstractVector
) where {S}
    # w_r += sum_r  A_rc v_c
    nrows, ncols = size(opr)
    if length(out) != nrows
        throw(DimensionMismatch("out has length $(length(out)) != dimension $(nrows)"))
    elseif length(state) != ncols
        throw(DimensionMismatch("state has length $(length(state)) != dimension $(ncols)"))
    end
    for irow in 1:nrows
        for (icol::Int, amplitude::S) in get_row_iterator(opr, irow)
            if 1 <= icol <= ncols
                @inbounds out[irow] += amplitude * state[icol]
            end
        end
    end
    return out
end


"""
    apply_serial!(out, state, opr)

Perform `out += state * opr`. Apply the operator representation `opr` to the
row vector `state` and *add* it to the row vector `out`.
Return sum of errors and sum of error-squared.
Single-threaded version.
"""
function apply_serial!(
    out::AbstractVector,
    state::AbstractVector,
    opr::AbstractOperatorRepresentation{S}
) where {S}
    # w_c += sum_r v_r A_rc
    nrows, ncols = size(opr)
    if length(out) != ncols
        throw(DimensionMismatch("out has length $(length(out)) != dimension $(ncols)"))
    elseif length(state) != nrows
        throw(DimensionMismatch("state has length $(length(state)) != dimension $(nrows)"))
    end
    for icol in 1:ncols
        for (irow::Int, amplitude::S) in get_column_iterator(opr, icol)
            if 1 <= irow <= nrows
                @inbounds out[icol] += state[irow] * amplitude
            end
        end
    end
    return out
end


"""
    apply_parallel!(out, opr, state)

Perform `out += opr * state`. Apply the operator representation `opr` to the
column vector `state` and *add* it to the column vector `out`.
Return sum of errors and sum of error-squared.
Multi-threaded version.
"""
function apply_parallel!(
    out::AbstractVector,
    opr::AbstractOperatorRepresentation{S},
    state::AbstractVector
) where {S}
    # w_r += sum_r  A_rc v_c
    nrows, ncols = size(opr)
    if length(out) != nrows
        throw(DimensionMismatch("out has length $(length(out)) != dimension $(nrows)"))
    elseif length(state) != ncols
        throw(DimensionMismatch("state has length $(length(state)) != dimension $(ncols)"))
    end
    Threads.@threads for irow in 1:nrows
        for (icol::Int, amplitude::S) in get_row_iterator(opr, irow)
            if 1 <= icol <= ncols
                @inbounds out[irow] += amplitude * state[icol]
            end
        end
    end
    return out
end


"""
    apply_parallel!(out, state, opr)

Perform `out += state * opr`. Apply the operator representation `opr` to the
row vector `state` and *add* it to the row vector `out`.
Return sum of errors and sum of error-squared.
Multi-threaded version.
"""
function apply_parallel!(
    out::AbstractVector,
    state::AbstractVector,
    opr::AbstractOperatorRepresentation{S}
) where {S}
    # w(c) += sum_(r) v(r) A(r,c)
    nrows, ncols = size(opr)
    if length(out) != ncols
        throw(DimensionMismatch("out has length $(length(out)) != dimension $(ncols)"))
    elseif length(state) != nrows
        throw(DimensionMismatch("state has length $(length(state)) != dimension $(nrows)"))
    end
    Threads.@threads for icol in 1:ncols
        for (irow::Int, amplitude::S) in get_column_iterator(opr, icol)
            if 1 <= irow <= nrows
                @inbounds out[icol] += state[irow] * amplitude
            end
        end
    end
    return out
end
