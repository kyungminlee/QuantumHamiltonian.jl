export prettyprintln

prettyprintln(x; kwargs...) = prettyprintln(stdout, x; kwargs...)

function prettyprintln(io::IO, op::NullOperator; prefix::AbstractString="")
  println(io, prefix, "NullOperator")
end

function prettyprintln(io::IO, op::PureOperator{S, BR}; prefix::AbstractString="") where {S, BR}
  println(io, prefix, "PureOperator")
  println(io, prefix, "| M: ", string(op.bitmask, base=2, pad=8*sizeof(BR)))
  println(io, prefix, "| R: ", string(op.bitrow, base=2, pad=8*sizeof(BR)))
  println(io, prefix, "| C: ", string(op.bitcol, base=2, pad=8*sizeof(BR)))
  println(io, prefix, "| A: ", op.amplitude)
end

function prettyprintln(io::IO, op::SumOperator{S, BR}; prefix::AbstractString="") where {S, BR}
  println(io, prefix, "SumOperator")
  for t in op.terms
    prettyprintln(io, t; prefix=string(prefix, "| "))
  end
end

function prettyprintln(io ::IO, hsr::HilbertSpaceRepresentation; prefix::AbstractString="")
  n = hsr.hilbert_space.bitoffsets[end]
  println(io, prefix, "HilbertSpaceRepresentation")
  for bvec in hsr.basis_list
    println(io, prefix, "| ", string(bvec, base=2, pad=n))
  end
end

function prettyprintln(io ::IO, rhsr::ReducedHilbertSpaceRepresentation; prefix::AbstractString="")
  n = rhsr.hilbert_space.bitoffsets[end]
  println(io, prefix, "ReducedHilbertSpaceRepresentation")
  println(io, prefix, "| basis_list")
  for bvec in rhsr.basis_list
    println(io, prefix, "| | ", string(bvec, base=2, pad=n))
  end
  println(io, prefix, "| basis_mapping")
  for (i_r, (i_p, ampl)) in enumerate(rhsr.basis_mapping)
    println(io, prefix, "| | ", i_r, ": ", i_p, ", ", ampl)
  end
end
