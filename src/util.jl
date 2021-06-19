export IntegerModulo

export make_bitmask
export choptol!
export merge_vec

"""
    IntegerModulo{N}

Implement Zₙ.
"""
struct IntegerModulo{N}
    value::Int
    IntegerModulo{N}(value::Integer) where N = new{N}(mod(value, N))
end


Base.:(+)(lhs::IntegerModulo{N}) where N = lhs
Base.:(-)(lhs::IntegerModulo{N}) where N = IntegerModulo{N}(-lhs.value)

Base.:(+)(lhs::IntegerModulo{N}, rhs::IntegerModulo{N}) where N = IntegerModulo{N}(lhs.value + rhs.value)
Base.:(+)(lhs::IntegerModulo{N}, rhs::Integer) where N = IntegerModulo{N}(lhs.value + rhs)
Base.:(+)(lhs::Integer, rhs::IntegerModulo{N}) where N = IntegerModulo{N}(lhs + rhs.value)

Base.:(-)(lhs::IntegerModulo{N}, rhs::IntegerModulo{N}) where N = IntegerModulo{N}(lhs.value - rhs.value)
Base.:(-)(lhs::IntegerModulo{N}, rhs::Integer) where N = IntegerModulo{N}(lhs.value - rhs)
Base.:(-)(lhs::Integer, rhs::IntegerModulo{N}) where N = IntegerModulo{N}(lhs - rhs.value)

Base.:(*)(lhs::IntegerModulo{N}, rhs::IntegerModulo{N}) where N = IntegerModulo{N}(lhs.value * rhs.value)
Base.:(*)(lhs::IntegerModulo{N}, rhs::Integer) where N = IntegerModulo{N}(lhs.value * rhs)
Base.:(*)(lhs::Integer, rhs::IntegerModulo{N}) where N = IntegerModulo{N}(lhs * rhs.value)

Base.:(==)(lhs::IntegerModulo{N}, rhs::IntegerModulo{N}) where N = lhs.value == rhs.value
Base.:(==)(lhs::IntegerModulo{N}, rhs::Integer) where N = lhs.value == rhs
Base.:(==)(lhs::Integer, rhs::IntegerModulo{N}) where N = lhs == rhs.value

tuplelength(::Type{<:NTuple{N, Any}}) where N = N
tuplelength(::NTuple{N, Any}) where N = N

tupleadd(l::T, r::T) where {T<:Tuple} = l .+ r
tuplesubtract(l::T, r::T) where {T<:Tuple} = l .- r

tuplezero(::Type{T}) where {T<:Tuple} = ((zero(S) for S in T.parameters)...,)
tupleone(::Type{T}) where {T<:Tuple} = ((one(S) for S in T.parameters)...,)

tuplezero(::T) where {T<:Tuple} = ((zero(S) for S in T.parameters)...,)
tupleone(::T) where {T<:Tuple} = ((one(S) for S in T.parameters)...,)

function make_bitmask(msb::Integer, ::Type{BR}=UInt) where {BR<:Unsigned}
    return one(BR) << msb - one(BR)
end

function make_bitmask(msb::Integer, lsb::Integer, ::Type{BR}=UInt) where {BR<:Unsigned}
    mask = one(BR) << msb - one(BR)
    submask = one(BR) << lsb - one(BR)
    return mask ⊻ submask
end

function merge_vec(x::AbstractVector{T}, y::AbstractVector{T})::Vector{T} where {T}
    (nx, ny) = (1, 1)
    z = T[]
    sizehint!(z, length(x) + length(y))
    while nx <= length(x) && ny <= length(y)
        if x[nx] < y[ny]
            push!(z, x[nx])
            nx += 1
        elseif y[ny] < x[nx]
            push!(z, y[ny])
            ny += 1
        else
            push!(z, x[nx])
            push!(z, y[ny])
            nx += 1
            ny += 1
        end
    end
    if nx <= length(x)
        append!(z, x[n] for n in nx:length(x))
    elseif ny <= length(y)
        append!(z, y[n] for n in ny:length(y))
    end
    return z
end

function choptol!(d::Dict{K, V}, tol::Real) where {K, V<:Number}
    to_delete = K[k for (k, v) in d if abs(v) < tol]
    for k in to_delete
        delete!(d, k)
    end
    return d
end


"""
    splitblock

Split n into b blocks.
"""
function splitblock(n::Integer, b::Integer)::Vector{Int}
    (n < 0) && throw(ArgumentError("n cannot be negative"))
    (b <= 0) && throw(ArgumentError("b must be positive"))
    blocksize = n ÷ b
    blocks = blocksize * ones(Int, b)
    r = n - (blocksize * b)
    blocks[1:r] .+= 1
    return blocks
end


"""
    splitrange
"""
function splitrange(range::AbstractVector{<:Integer}, b::Integer)
    nblocks = b
    block_sizes = splitblock(length(range), nblocks)
    block_offsets = cumsum([1, block_sizes...])
    return [range[block_offsets[iblock]:block_offsets[iblock+1]-1] for iblock in 1:nblocks]
end


"""
    compress(bitwidths, data, BR)

Compress data array into a binary integer of type BR.
"""
function compress(bitwidths::AbstractVector{<:Integer}, data::AbstractVector{<:Integer}, ::Type{BR}) where {BR<:Unsigned}
    n = length(bitwidths)
    @boundscheck if length(data) != n
        throw(ArgumentError("lengths of bitwidths and data have to match ($(length(bitwidths)) != $(length(data)))"))
    end
    bitoffsets = [0, cumsum(bitwidths)...]
    @boundscheck if bitoffsets[end] > sizeof(BR) * 8
        throw(ArgumentError("type $BR too small to represent data. Need $(bitoffsets[end]) bits."))
    end

    out = zero(BR)
    for i in 1:n
        @boundscheck if data[i] < 0 || data[i] >= (1<<bitwidths[i])
            throw(ArgumentError("value $(data[i]) too large to be represented with $(bitwidths[i]) bits."))
        end
        out |= BR(data[i]) << bitoffsets[i]
    end
    return out
end


# # Not used (as of 2020.04.21)
# export elmax, elmin
# export elmaximum, elminimum
#
# elmax(x::S, y::S) where {S<:Number} = max(x,y)
# elmax(x::T, y::T) where {T<:Tuple{<:Number}} = (max(first(x), first(y)),)
#
# function elmax(x::T, y::T) where {T<:Tuple{<:Number, <:Number, Vararg{<:Number}}}
#   return (max(first(x), first(y)), elmax(x[2:end], y[2:end])...)
# end
#
# elmin(x::S, y::S) where {S<:Number} = max(x,y)
# elmin(x::T, y::T) where {T<:Tuple{<:Number}} = (max(first(x), first(y)),)
#
# function elmin(x::T, y::T) where {T<:Tuple{<:Number, <:Number, Vararg{<:Number}}}
#   return (max(first(x), first(y)), elmax(x[2:end], y[2:end])...)
# end
#
# elmaximum(arr::AbstractArray) = mapreduce(identity, elmax, arr)
# elminimum(arr::AbstractArray) = mapreduce(identity, elmin, arr)
