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
        lookup = Dict(v => k for (k, v) in enumerate(elements))
        return new{E}(elements, lookup)
    end
    function DictIndexedVector{E}(elements) where {E}
        lookup = Dict(v => k for (k, v) in enumerate(elements))
        return new{E}(elements, lookup)
    end
    function DictIndexedVector{E}(elements::AbstractVector{E}, lookup::AbstractDict{E, Int}) where {E}
        return new{E}(elements, lookup)
    end
end

findindex(obj::DictIndexedVector{E}, key::E) where E = get(obj.lookup, key, -1)
Base.getindex(obj::DictIndexedVector, idx::Integer) = obj.elements[idx]
Base.in(obj::DictIndexedVector, key) = haskey(obj.lookup, key)

Base.iterate(obj::DictIndexedVector) = iterate(obj.elements)
Base.iterate(obj::DictIndexedVector, state) = iterate(obj.elements, state)

Base.IteratorSize(::Type{<:DictIndexedVector}) = Base.HasLength()
Base.IteratorEltype(::Type{<:DictIndexedVector}) = Base.HasEltype()
Base.eltype(::Type{DictIndexedVector{E}}) where E = E
Base.length(obj::DictIndexedVector) = length(obj.elements)
Base.size(obj::DictIndexedVector, args...) = size(obj.elements, args...)


"""
"""
struct SortedIndexedVector{E}
    elements::Vector{E}

    function SortedIndexedVector{E}(elements) where E
        issorted(elements) || throw(ArgumentError("elements must be sorted"))
        if length(elements) > 1
            for i in 2:length(elements)
                @inbounds elements[i-1] == elements[i] && throw(ArgumentError("elements contains duplicates $(keys[i])"))
            end
        end
        return new{E}(elements)
    end

    function SortedIndexedVector(elements::AbstractVector{E}) where E
        issorted(elements) || throw(ArgumentError("elements must be sorted"))
        if length(elements) > 1
            for i in 2:length(elements)
                @inbounds elements[i-1] == elements[i] && throw(ArgumentError("elements contains duplicates $(keys[i])"))
            end
        end
        return new{E}(elements)
    end
end


@noinline function findindex(arr::SortedIndexedVector, val)::Int
    idx = searchsortedfirst(arr.elements, val)
    return (idx <= length(arr.elements) && @inbounds arr.elements[idx] == val) ? idx : -1
end
Base.getindex(obj::SortedIndexedVector, idx::Integer) = obj.elements[idx]
Base.in(obj::SortedIndexedVector, val) = (findindex(obj, val) > 0)

Base.iterate(obj::SortedIndexedVector) = iterate(obj.elements)
Base.iterate(obj::SortedIndexedVector, state) = iterate(obj.elements, state)

Base.IteratorSize(::Type{<:SortedIndexedVector}) = Base.HasLength()
Base.IteratorEltype(::Type{<:SortedIndexedVector}) = Base.HasEltype()
Base.eltype(::Type{SortedIndexedVector{E}}) where E = E
Base.length(obj::SortedIndexedVector) = length(obj.elements)
Base.size(obj::SortedIndexedVector, args...) = size(obj.elements, args...)
