export ReducedOperatorRepresentation
export bintype
export represent
export get_row_iterator
export get_column_iterator
export get_row
export get_column


"""
    ReducedOperatorRepresentation{RHSR, O, S, BR}

Representation of an operator of type `O` in the symmetry-reduced hilbert space representation of type `RHSR`.
"""
struct ReducedOperatorRepresentation{
    RHSR<:ReducedHilbertSpaceRepresentation,
    O<:AbstractOperator,
    S<:Number,
    BR<:Unsigned
}<:AbstractOperatorRepresentation{S}

    reduced_hilbert_space_representation::RHSR
    operator::O

    function ReducedOperatorRepresentation(rhsr::RHSR, op::O) where {RHSR<:ReducedHilbertSpaceRepresentation, O<:AbstractOperator}
        S = promote_type(valtype(RHSR), valtype(O))
        BR = bintype(RHSR)
        return new{RHSR, O, S, BR}(rhsr, op)
    end
end


spacetype(::Type{ReducedOperatorRepresentation{RHSR, O, S, BR}}) where {RHSR, O, S, BR} = RHSR
operatortype(::Type{ReducedOperatorRepresentation{RHSR, O, S, BR}}) where {RHSR, O, S, BR} = O
function get_space(lhs::ReducedOperatorRepresentation{RHSR, O, S, BR})::RHSR where {RHSR, O, S, BR}
    return lhs.reduced_hilbert_space_representation
end


function represent(rhsr::RHSR, op::O) where {RHSR<:ReducedHilbertSpaceRepresentation, O<:AbstractOperator}
    return ReducedOperatorRepresentation(rhsr, op)
end


function Base.show(
    io::IO, ::MIME"text/plain", arg::ReducedOperatorRepresentation{RHSR, O, S, BR}
) where {RHSR, O, S, BR}
    print(io, string(typeof(arg)))
    print(io, "(", arg.reduced_hilbert_space_representation, ", ", arg.operator, ")")
end

"""
    get_row_iterator(ropr::ROR, irow_r::Integer)

Get the row iterator for the reduced operator representation
"""
function get_row_iterator(
    opr::ReducedOperatorRepresentation{RHSR, O, S, BR},
    irow_r::Integer
) where {RHSR, O, S, BR}
    rhsr = opr.reduced_hilbert_space_representation
    hsr = rhsr.parent

    brow::BR = rhsr.basis[irow_r]
    irow_p::Int = findindex(hsr.basis, brow)
    ampl_row::S = rhsr.basis_mapping_amplitude[irow_p]

    # H_r = U H U†
    full_iter = let inv_ampl_row = one(S) / ampl_row,
                    # basis_lookup = rhsr.parent.basis_lookup,
                    basis = rhsr.parent.basis,
                    basis_mapping_index::Vector{Int} = rhsr.basis_mapping_index,
                    basis_mapping_amplitude::Vector{S} = rhsr.basis_mapping_amplitude,
                    operator = opr.operator
        function element(bcol::BR, ampl::S)::Pair{Int, S}
            # icol_p = get(basis_lookup, bcol, -1)
            icol_p = findindex(basis, bcol)
            (icol_p > 0) || return (-1 => ampl)
            icol_r = basis_mapping_index[icol_p]
            (icol_r > 0) || return (-1 => ampl)
            ampl_col = basis_mapping_amplitude[icol_p]
            return icol_r => ampl * ampl_col * inv_ampl_row
        end
        (element(bcol, ampl) for (bcol::BR, ampl::S) in get_row_iterator(operator, brow))
    end
    return full_iter
end


function get_column_iterator(
    opr::ReducedOperatorRepresentation{RHSR, O, S, BR},
    icol_r::Integer
) where {RHSR, O, S, BR}
    rhsr = opr.reduced_hilbert_space_representation
    hsr = rhsr.parent

    bcol::BR = rhsr.basis[icol_r]
    # icol_p::Int = hsr.basis_lookup[bcol]
    icol_p::Int = findindex(hsr.basis, bcol)
    ampl_col::S = rhsr.basis_mapping_amplitude[icol_p]

    # H_r = U† H U
    full_iter = let inv_ampl_col = one(S) / conj(ampl_col),
                    # basis_lookup = rhsr.parent.basis_lookup,
                    basis = rhsr.parent.basis,
                    basis_mapping_index::Vector{Int} = rhsr.basis_mapping_index,
                    basis_mapping_amplitude::Vector{S} = rhsr.basis_mapping_amplitude,
                    operator = opr.operator
        function element(brow::BR, ampl::S)::Pair{Int, S}
            # irow_p = get(basis_lookup, brow, -1)
            irow_p = findindex(basis, brow)
            (irow_p > 0) || return (-1 => ampl)
            irow_r = basis_mapping_index[irow_p]
            (irow_r > 0) || return (-1 => ampl)
            ampl_row = conj(basis_mapping_amplitude[irow_p])
            return irow_r => ampl * ampl_row * inv_ampl_col
        end
        (element(brow, ampl) for (brow::BR, ampl::S) in get_column_iterator(operator, bcol))
    end
    return full_iter
end


# TODO: better implementation?
function get_element(
    opr::ReducedOperatorRepresentation{RHSR, O, S, BR},
    irow_r::Integer, icol_r::Integer
) where {RHSR, O, S, BR}
    rhsr = opr.reduced_hilbert_space_representation
    @boundscheck let
        dim = length(rhsr.basis)
        if irow_r <= 0 || irow_r > dim || icol_r <= 0 || icol_r > dim
            throw(BoundsError(opr, [irow_r, icol_r]))
        end
    end
    #return @inbounds sum(ampl for (irow_r2, ampl::S) in get_column_iterator(opr, icol_r) if irow_r2 == irow_r)
    return @inbounds mapreduce(
        identity,
        +,
        ampl for (irow_r2, ampl::S) in get_column_iterator(opr, icol_r) if irow_r2 == irow_r;
        init=zero(S)
    )
end
