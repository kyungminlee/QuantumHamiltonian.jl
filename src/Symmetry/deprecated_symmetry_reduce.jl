for symred in [
    :symmetry_reduce,
    :symmetry_reduce_serial,
    :symmetry_reduce_parallel
]
    @eval begin
        function $symred(
            hsr::AbstractHilbertSpaceRepresentation,
            psic::AbstractSymmetryIrrepComponent,
            ::Type{BT}=DictIndexedVector;
            tol::Real=Base.rtoldefault(Float64)
        ) where {BT<:AbstractIndexedVector}
            @warn "Use of IrrepComponent with symmetry_reduce deprecated. Use simpler symmetry_reduce function that takes symops and amplitudes as input" maxlog=1
            symops_and_amplitudes = [(x, y) for (x, y) in get_irrep_iterator(psic) if !isapprox(y, zero(y); atol=tol)]
            symops = [x for (x, y) in symops_and_amplitudes]
            amplitudes = [y for (x, y) in symops_and_amplitudes]
            return $symred(hsr, symops, amplitudes, BT; tol=tol)
        end

        function $symred(
            hsr::AbstractHilbertSpaceRepresentation,
            symops_and_amplitudes::AbstractArray{<:Tuple{<:AbstractSymmetryOperation, <:Number}},
            ::Type{BT}=DictIndexedVector;
            tol::Real=Base.rtoldefault(Float64)
        ) where {BT<:AbstractIndexedVector}
            @warn "Use of array of (operator, phase) tuple with symmetry_reduce deprecated. Use simpler symmetry_reduce function that takes symops and amplitudes as input" maxlog=1
            symops = [x for (x, y) in symops_and_amplitudes]
            amplitudes = [y for (x, y) in symops_and_amplitudes]
            return $symred(hsr, symops, amplitudes, BT; tol=tol)
        end
    end
end