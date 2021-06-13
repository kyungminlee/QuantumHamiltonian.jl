export OperatorRepresentation
export represent
export apply!, apply_serial!, apply_parallel!

import LinearAlgebra


"""
    OperatorRepresentation{HSR, S, O}

Operator representation of given operator of type `O`.
"""
struct OperatorRepresentation{
    HSR<:HilbertSpaceRepresentation,
    S<:Number,
    O<:AbstractOperator
}<:AbstractOperatorRepresentation{S}

    hilbert_space_representation::HSR
    operator::O

    function OperatorRepresentation(hsr::HSR, op::O) where {HSR<:HilbertSpaceRepresentation, O<:AbstractOperator}
        S = valtype(op)
        return new{HSR, S, O}(hsr, op)
    end
end


"""
    represent(hilbert_space_representation, operator)

Create an `OperatorRepresentation` of the `operator` in the `hilbert_space_representation`.
"""
function represent(hsr::HSR, op::O) where {HSR<:HilbertSpaceRepresentation, O<:AbstractOperator}
    return OperatorRepresentation(hsr, op)
end


spacetype(::Type{OperatorRepresentation{HSR, S, O}}) where {HSR, S, O} = HSR
operatortype(::Type{OperatorRepresentation{HSR, S, O}}) where {HSR, S, O} = O
function get_space(lhs::OperatorRepresentation{HSR, S, O})::HSR where {HSR, S, O}
    return lhs.hilbert_space_representation
end


function LinearAlgebra.issymmetric(arg::OperatorRepresentation)
    return LinearAlgebra.issymmetric(arg.operator)
end


function Base.show(io::IO, ::MIME"text/plain", arg::OperatorRepresentation{HSR, S, O}) where {HSR, S, O}
    print(io, string(typeof(arg)), "(", arg.hilbert_space_representation, ", ", arg.operator, ")")
end


## iterators

"""
    get_row_iterator(opr, irow)

Returns an iterator which generates a list of elements of the row `irow`.
Each element is represented as (icol, amplitude).
**May contain duplicates and invalid elements.**
Invalid elements are represented as (-1, amplitude).
"""
function get_row_iterator(
    opr::OperatorRepresentation{HSR, S, O},
    irow::Integer
) where {HSR, S, O}
    hsr = opr.hilbert_space_representation
    brow = hsr.basis_list[irow]
    basis_lookup = hsr.basis_lookup
    operator = opr.operator
    iter = (
        get(basis_lookup, bcol, -1) => amplitude
            for (bcol, amplitude) in get_row_iterator(operator, brow)
    )
    return iter
end


"""
    get_column_iterator(opr, icol)

Returns an iterator which generates a list of elements of the column `icol`.
Each element is represented as (irow, amplitude).
**May contain duplicates and invalid elements.**
Invalid elements are represented as (-1, amplitude).
"""
function get_column_iterator(
    opr::OperatorRepresentation{HSR, S, O},
    icol::Integer
) where {HSR, S, O}
    hsr = opr.hilbert_space_representation
    bcol = hsr.basis_list[icol]
    basis_lookup = hsr.basis_lookup
    operator = opr.operator
    iter = (
        get(basis_lookup, brow, -1) => amplitude
            for (brow, amplitude) in get_column_iterator(operator, bcol)
    )
    return iter
end


"""
    get_element(opr, irow, icol)
"""
function get_element(opr::OperatorRepresentation{HSR, S, O}, irow::Integer, icol::Integer) where {HSR, S, O}
    hsr = opr.hilbert_space_representation
    @boundscheck let dim = length(hsr.basis_list)
        if irow <= 0 || irow > dim || icol <= 0 || icol > dim
            throw(BoundsError(opr, [irow, icol]))
        end
    end
    @inbounds brow = hsr.basis_list[irow]
    @inbounds bcol = hsr.basis_list[icol]
    return get_element(opr.operator, brow, bcol)
end


function Base.:(*)(opr::OperatorRepresentation{HSR, SO, O}, state::AbstractVector{SV}) where {HSR, O, SO, SV<:Number}
    hsr = opr.hilbert_space_representation
    n = dimension(hsr)
    T = promote_type(SO, SV)
    out = zeros(T, n)
    err = apply!(out, opr, state)
    return out
end


function Base.:(*)(state::AbstractVector{SV}, opr::OperatorRepresentation{HSR, SO, O}) where {HSR, SO, O, SV<:Number}
    hsr = opr.hilbert_space_representation
    n = dimension(hsr)
    T = promote_type(SO, SV)
    out = zeros(T, n)
    err = apply!(out, state, opr)
    return out
end


function Base.convert(::Type{AbstractMatrix{T}}, arg::OperatorRepresentation) where {T}
    return OperatorRepresentation(arg.hilbert_space_representation, convert(AbstractOperator{T}, arg.operator))
end
