export AbstractIndexedVector
export DictIndexedVector
export SortedIndexedVector
export findindex


abstract type AbstractIndexedVector{E} <: AbstractVector{E} end

"""
    findindex(a::AbstractIndexedVector{E}, key::E)

Return the index of the element `key` if it exists in `a`. Otherwise return -1.
"""
function findindex end

"""
"""
struct DictIndexedVector{E} <: AbstractIndexedVector{E}
    elements::Vector{E}
    lookup::Dict{E, Int}
    function DictIndexedVector(elements::AbstractVector{E}) where {E}
        lookup = Dict{E, Int}()
        sizehint!(lookup, length(elements) + length(elements) ÷ 4)
        for (i, v) in enumerate(elements)
            lookup[v] = i
        end
        return new{E}(elements, lookup)
    end
    function DictIndexedVector{E}(elements) where {E}
        lookup = Dict{E, Int}()
        sizehint!(lookup, length(elements) + length(elements) ÷ 4)
        for (i, v) in enumerate(elements)
            lookup[v] = i
        end
        return new{E}(elements, lookup)
    end
    function DictIndexedVector{E}(elements::AbstractVector{E}, lookup::AbstractDict{E, Int}) where {E}
        return new{E}(elements, lookup)
    end
end

@inline findindex(obj::DictIndexedVector, val) = get(obj.lookup, val, -1)
Base.getindex(obj::DictIndexedVector, idx::Integer) = obj.elements[idx]
Base.in(val, obj::DictIndexedVector) = haskey(obj.lookup, val)

Base.iterate(obj::DictIndexedVector) = iterate(obj.elements)
Base.iterate(obj::DictIndexedVector, state) = iterate(obj.elements, state)

# Base.IteratorEltype(::Type{<:DictIndexedVector}) = Base.HasEltype()
# Base.IteratorSize(::Type{<:DictIndexedVector}) = Base.HasLength()
Base.eltype(::Type{DictIndexedVector{E}}) where E = E
Base.length(obj::DictIndexedVector) = length(obj.elements)
Base.size(obj::DictIndexedVector, args...) = size(obj.elements, args...)

Base.hash(obj::V, h::UInt) where {V<:DictIndexedVector} = hash(V, hash(obj.elements, h))
Base.:(==)(lhs::DictIndexedVector, rhs::DictIndexedVector) = lhs.elements == rhs.elements
# Base.:(==)(lhs::DictIndexedVector, rhs::AbstractVector) = lhs.elements == rhs
# Base.:(==)(lhs::AbstractVector, rhs::DictIndexedVector) = lhs == rhs.elements

function checkvalid(obj::DictIndexedVector)
    for (i, x) in enumerate(obj.elements)
        @assert obj.lookup[x] == i
    end
    for (x, i) in obj.lookup
        @assert obj.elements[i] == x
    end
end


"""
"""
struct SortedIndexedVector{E} <: AbstractIndexedVector{E}
    elements::Vector{E}

    function SortedIndexedVector{E}(elements) where E
        issorted(elements) || throw(ArgumentError("elements must be sorted"))
        for i in 2:length(elements)
            @inbounds elements[i-1] == elements[i] && throw(ArgumentError("elements contains duplicates $(keys[i])"))
        end
        return new{E}(elements)
    end

    function SortedIndexedVector(elements::AbstractVector{E}) where E
        issorted(elements) || throw(ArgumentError("elements must be sorted"))
        for i in 2:length(elements)
            @inbounds elements[i-1] == elements[i] && throw(ArgumentError("elements contains duplicates $(keys[i])"))
        end
        return new{E}(elements)
    end
end


@noinline function findindex(arr::SortedIndexedVector, val)::Int
    idx = searchsortedfirst(arr.elements, val)
    return (idx <= length(arr.elements) && @inbounds arr.elements[idx] == val) ? idx : -1
end
Base.getindex(obj::SortedIndexedVector, idx::Integer) = obj.elements[idx]
Base.in(val, obj::SortedIndexedVector) = (findindex(obj, val) > 0)

Base.iterate(obj::SortedIndexedVector) = iterate(obj.elements)
Base.iterate(obj::SortedIndexedVector, state) = iterate(obj.elements, state)

# Base.IteratorEltype(::Type{<:SortedIndexedVector}) = Base.HasEltype()
# Base.IteratorSize(::Type{<:SortedIndexedVector}) = Base.HasLength()
Base.eltype(::Type{SortedIndexedVector{E}}) where E = E
Base.length(obj::SortedIndexedVector) = length(obj.elements)
Base.size(obj::SortedIndexedVector, args...) = size(obj.elements, args...)

Base.hash(obj::V, h::UInt) where {V<:SortedIndexedVector} = hash(V, hash(obj.elements, h))
Base.:(==)(lhs::SortedIndexedVector, rhs::SortedIndexedVector) = lhs.elements == rhs.elements
# Base.:(==)(lhs::SortedIndexedVector, rhs::AbstractVector) = lhs.elements == rhs
# Base.:(==)(lhs::AbstractVector, rhs::SortedIndexedVector) = lhs == rhs.elements

function checkvalid(obj::SortedIndexedVector)
    @assert issorted(obj.elements)
    for i in 2:length(obj.elements)
        @assert obj.elements[i-1] != obj.elements[i]
    end
end
