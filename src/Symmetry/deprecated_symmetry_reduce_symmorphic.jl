export symmetry_reduce_serial, symmetry_reduce_parallel

"""
    symmetry_reduce_serial(hsr, trans_group, frac_momentum, complex_type=ComplexF64, tol=√ϵ)

Symmetry-reduce the HilbertSpaceRepresentation using translation group (single threaded).

"""
function symmetry_reduce_serial(
    hsr::HilbertSpaceRepresentation,
    ssic::SymmorphicIrrepComponent{
        <:SymmetryEmbedding{<:FiniteTranslationSymmetry},
        <:SymmetryEmbedding{<:PointSymmetry}
    };
    tol::Real=Base.rtoldefault(Float64)
)
    @warn "Method deprecated. use simpler symmetry_reduce functions that takes symops_and_amplitudes as input" maxlog=1
    tsym_symops_and_amplitudes = [(x, y) for (x, y) in get_irrep_iterator(ssic.normal)]
    tsym_group_size = group_order(ssic.normal.symmetry)
    @assert length(tsym_symops_and_amplitudes) == tsym_group_size

    psym_symops_and_amplitudes = [
        (x, y)
            for (x, y) in get_irrep_iterator(ssic.rest)
                if !isapprox(y, zero(y); atol=tol)
    ]
    psym_subgroup_size = length(psym_symops_and_amplitudes)
    @assert psym_subgroup_size <= group_order(ssic.rest.symmetry)

    symops_and_amplitudes = [
        (p*t, ϕp*ϕt)
            for (t, ϕt) in tsym_symops_and_amplitudes
                #if !isapprox(ϕt, zero(ϕt); atol=tol)
            for (p, ϕp) in psym_symops_and_amplitudes
                if !isapprox(ϕp, zero(ϕp); atol=tol)
    ]
    return symmetry_reduce_serial(hsr, symops_and_amplitudes; tol=tol)
end




"""
    symmetry_reduce_parallel(hsr, trans_group, frac_momentum, complex_type=ComplexF64, tol=√ϵ)

Symmetry-reduce the HilbertSpaceRepresentation using translation group (multi-threaded).

"""
function symmetry_reduce_parallel(
    hsr::HilbertSpaceRepresentation,
    ssic::SymmorphicIrrepComponent{
        <:SymmetryEmbedding{<:FiniteTranslationSymmetry},
        <:SymmetryEmbedding{<:PointSymmetry}
    };
    tol::Real=Base.rtoldefault(Float64)
)
    tsym_symops_and_amplitudes = [(x, y) for (x, y) in get_irrep_iterator(ssic.normal)]
    tsym_group_size = group_order(ssic.normal.symmetry)
    @assert length(tsym_symops_and_amplitudes) == tsym_group_size

    psym_symops_and_amplitudes = [(x, y) for (x, y) in get_irrep_iterator(ssic.rest) if !isapprox(y, zero(y); atol=tol)]
    psym_subgroup_size = length(psym_symops_and_amplitudes)
    @assert psym_subgroup_size <= group_order(ssic.rest.symmetry)

    symops_and_amplitudes = [
        (p*t, ϕp*ϕt)
            for (t, ϕt) in tsym_symops_and_amplitudes #if !isapprox(ϕt, zero(ϕt); atol=tol)
            for (p, ϕp) in psym_symops_and_amplitudes if !isapprox(ϕp, zero(ϕp); atol=tol)
    ]
    return symmetry_reduce_parallel(hsr, symops_and_amplitudes; tol=tol)
end
