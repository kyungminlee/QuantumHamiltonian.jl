export HilbertSpaceSector
export scalartype
export qntype
export basespace
export bitwidth


"""
    HilbertSpaceSector{QN}

Hilbert space sector.
"""
struct HilbertSpaceSector{
    QN<:Tuple{Vararg{AbstractQuantumNumber}},
    HS<:AbstractHilbertSpace{QN}
}<:AbstractHilbertSpace{QN}

    parent::HS
    allowed_quantum_numbers::Set{QN}

    """
        HilbertSpaceSector(parent::HilbertSpace)
    """
    function HilbertSpaceSector(parent::HS) where {HS<:AbstractHilbertSpace}
        @warn "Use of HilbertSpaceSector deprecated" maxlog=1
        QN = qntype(HS)
        sectors = quantum_number_sectors(parent)
        return new{QN, HS}(parent, Set(sectors))
    end

    """
        HilbertSpaceSector(parent::HilbertSpace{QN, TT}, allowed::Integer) where {QN<:Tuple{<:Integer}}
    """
    function HilbertSpaceSector(parent::HS, allowed::Integer) where {HS<:AbstractHilbertSpace{<:Tuple{<:Integer}}}
        @warn "Use of HilbertSpaceSector deprecated" maxlog=1
        QN = qntype(HS)
        sectors = Set{QN}(quantum_number_sectors(parent))
        return new{QN, HS}(parent, intersect(sectors, Set([(allowed,)])))
    end

    """
        HilbertSpaceSector(parent::HilbertSpace{QN, TT}, allowed::QN)
    """
    function HilbertSpaceSector(parent::AbstractHilbertSpace{QN}, allowed::QN) where {QN}
        @warn "Use of HilbertSpaceSector deprecated" maxlog=1
        HS = typeof(parent)
        sectors = Set{QN}(quantum_number_sectors(parent))
        return new{QN, HS}(parent, intersect(sectors, Set([allowed])))
    end

    """
        HilbertSpaceSector(parent::HilbertSpace{QN, TT}, allowed::Union{AbstractSet{<:Integer}, AbstractVector{<:Integer}})
    """
    function HilbertSpaceSector(
        parent::AbstractHilbertSpace{QN},
        allowed::Union{AbstractSet{<:Integer}, AbstractVector{<:Integer}}
    ) where {QN}
        @warn "Use of HilbertSpaceSector deprecated" maxlog=1
        HS = typeof(parent)
        sectors = Set{QN}(quantum_number_sectors(parent))
        return new{QN, HS}(parent, intersect(sectors, Set((x,) for x in allowed)))
    end

    function HilbertSpaceSector(
        parent::AbstractHilbertSpace{QN},
        allowed::Union{AbstractSet{QN}, AbstractVector{QN}}
    ) where QN
        @warn "Use of HilbertSpaceSector deprecated" maxlog=1
        HS = typeof(parent)
        sectors = Set{QN}(quantum_number_sectors(parent))
        return new{QN, HS}(parent, intersect(sectors, Set(allowed)))
    end
end


"""
    qntype(arg::Type{HilbertSpaceSector{QN}})

Returns the quantum number type of the given hilbert space sector type.
"""
qntype(::Type{<:HilbertSpaceSector{QN, HS}}) where {QN, HS} = QN
tagtype(::Type{<:HilbertSpaceSector{QN, HS}}, strategy::Val) where {QN, HS} = tagtype(HS, strategy)

function Base.:(==)(lhs::HilbertSpaceSector{QN, HS}, rhs::HilbertSpaceSector{QN, HS}) where {QN, HS}
    return (
        basespace(lhs) == basespace(rhs) &&
        lhs.allowed_quantum_numbers == rhs.allowed_quantum_numbers
    )
end


get_quantum_numbers(hss::HilbertSpaceSector) = sort(collect(hss.allowed_quantum_numbers))
get_tags(hss::HilbertSpaceSector, ::Val{:QuantumNumberAsTag}) = get_quantum_numbers(hss)
quantum_number_sectors(hss::HilbertSpaceSector) = get_quantum_numbers(hss)


for fname in [
    :basespace,
    :numsites,
    :get_site,
    :bitwidth,
    :bitoffset,
    :get_quantum_number,
    :get_tag,
    :extract,
    :compress,
]
    @eval begin
        """
            $($fname)(hss::HilbertSpaceSector, args...;kwargs...)

        Call `$($fname)` with basespace of `hss`.
        """
        @inline $fname(hss::HilbertSpaceSector, args...;kwargs...) = $fname(hss.parent, args...; kwargs...)
    end
end

Base.keys(hs::HilbertSpaceSector) = error("Base.keys unsupported for HilbertSpaceSector") # COV_EXCL_LINE
