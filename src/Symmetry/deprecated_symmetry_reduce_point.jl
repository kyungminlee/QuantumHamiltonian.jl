export symmetry_reduce_serial, symmetry_reduce_parallel


"""
    symmetry_reduce_serial(hsr, trans_group, frac_momentum, complex_type=ComplexF64, tol=√ϵ)

Symmetry-reduce the HilbertSpaceRepresentation using translation group (single threaded).

"""
function symmetry_reduce_serial(
    hsr::HilbertSpaceRepresentation,
    psic::IrrepComponent{SymmetryEmbedding{PointSymmetry}};
    tol::Real=Base.rtoldefault(Float64)
)
    @warn "Method deprecated. use simpler symmetry_reduce functions that takes symops_and_amplitudes as input" maxlog=1
    symops_and_amplitudes = [(x, y) for (x, y) in get_irrep_iterator(psic) if !isapprox(y, zero(y); atol=tol)]
    return symmetry_reduce_serial(hsr, symops_and_amplitudes; tol=tol)
end


"""
    symmetry_reduce_parallel(hsr, trans_group, frac_momentum, complex_type=ComplexF64, tol=√ϵ)

Symmetry-reduce the HilbertSpaceRepresentation using translation group (multi-threaded).

"""
function symmetry_reduce_parallel(
    hsr::HilbertSpaceRepresentation,
    psic::IrrepComponent{SymmetryEmbedding{PointSymmetry}};
    tol::Real=Base.rtoldefault(Float64)
)
    @warn "Method deprecated. use simpler symmetry_reduce functions that takes symops_and_amplitudes as input" maxlog=1
    symops_and_amplitudes = [(x, y) for (x, y) in get_irrep_iterator(psic) if !isapprox(y, zero(y); atol=tol)]
    return symmetry_reduce_parallel(hsr, symops_and_amplitudes; tol=tol)
end
