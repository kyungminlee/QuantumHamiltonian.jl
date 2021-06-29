export simplify


"""
    simplify

Simplify the given operator.
"""
simplify(op::NullOperator; tol::Real=0) = op


function simplify(op::PureOperator{S, BR}; tol::Real=Base.rtoldefault(S)) where {S<:Real, BR}
    if isapprox(op.amplitude, zero(S); atol=tol)
        return NullOperator()
    else
        return op
    end
end


function simplify(op::PureOperator{Complex{S}, BR}; tol::Real=Base.rtoldefault(S)) where {S<:Real, BR}
    if isapprox(abs(op.amplitude), zero(S); atol=tol)
        return NullOperator()
    elseif isapprox(imag(op.amplitude), zero(S); atol=tol)
        return real(op)
    else
        return op
    end
end


function simplify(so::SumOperator{S, BR}; tol::Real=Base.rtoldefault(real(S))) where {S, BR}
    terms::Vector{PureOperator{S, BR}} = filter((x) -> !isa(x, NullOperator), simplify.(so.terms))

    isempty(terms) && return NullOperator()

    sort!(terms; lt=isless)
    new_terms = PureOperator{S, BR}[]

    bm::BR = terms[1].bitmask
    br::BR = terms[1].bitrow
    bc::BR = terms[1].bitcol
    am::S  = terms[1].amplitude

    for term in terms[2:end]
        if (bm == term.bitmask) && (br == term.bitrow) && (bc == term.bitcol)
            am += term.amplitude
        else
            if !isapprox(am, zero(S); atol=tol)
                push!(new_terms, PureOperator{S, BR}(bm, br, bc, am))
            end
            bm = term.bitmask
            br = term.bitrow
            bc = term.bitcol
            am = term.amplitude
        end
    end

    if !isapprox(am, zero(S); atol=tol)
        push!(new_terms, PureOperator{S, BR}(bm, br, bc, am))
    end

    if isempty(new_terms)
        return NullOperator()
    end

    R = real(S)
    if S <: Complex && isapprox(maximum(abs(imag(t.amplitude)) for t in new_terms), zero(R); atol=tol)
        if length(new_terms) == 1
            return real(new_terms[1])
        else
            return SumOperator{R, BR}(real.(new_terms))
        end
    else
        if length(new_terms) == 1
            return new_terms[1]
        else
            return SumOperator{S, BR}(new_terms)
        end
    end
end


# function simplify_combine(so ::SumOperator{S, BR}; tol::Real=Base.rtoldefault(Float64)) where {S, BR}
#   bw = sizeof(BR)*8
#   diagonal_terms = [t for t in so.terms if t.bitrow == t.bitcol]
#   offdiagonal_terms = [t for t in so.terms if t.bitrow != t.bitcol]
#
#   @label simplify_combine_loop_start
#
#   for ib in 0:(bw-1)
#     mask = make_bitmask(ib+1, ib)
#
#     for (i1, t1) in enumerate(diagonal_terms)
#       (t1.bitmask & mask) == 0x0 && continue
#       t1.bitrow != t1.bitcol && continue
#
#       for i2 in (i1+1):length(diagonal_terms)
#         t2 = diagonal_terms[i2]
#         (t2.bitmask & mask) == 0x0 && continue
#         t2.bitrow != t2.bitcol && continue
#
#         if ( ((t1.bitmask & ~mask) == (t2.bitmask & ~mask)) &&
#              ((t1.bitrow  & ~mask) == (t2.bitrow  & ~mask)) &&
#              ((t1.bitrow  &  mask) != (t2.bitrow  &  mask)) )
#
#           if isapprox(t1.amplitude, t2.amplitude; atol=tol)
#             i_small, i_large = i1, i2
#             new_bitmask = t1.bitmask & (~mask)
#             new_bitrow  = t1.bitrow  & (~mask)
#             new_bitcol  = t1.bitcol  & (~mask)
#             new_amplitude = t1.amplitude
#
#             diagonal_terms[i_small] = PureOperator{S, BR}(new_bitmask, new_bitrow, new_bitcol, new_amplitude)
#             deleteat!(diagonal_terms, i_large)
#
#             @goto simplify_combine_loop_start
#             #dirty = true
#             #break
#           else
#             # different amplitude. do nothing
#           end # if amplitude
#         end # if bits
#
#       end # for i2
#     end # for i1
#   end # for ib
#
#   #end # while dirty
#   SumOperator{S, BR}([diagonal_terms..., offdiagonal_terms])
# end
