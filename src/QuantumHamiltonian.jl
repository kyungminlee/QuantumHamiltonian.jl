module QuantumHamiltonian

using StaticArrays
using SparseArrays

# workaround?
import LatticeTools.dimension

include("util.jl")
include("indexedvector.jl")

include("HilbertSpace.jl")
include("Operator.jl")
include("Representation.jl")
include("Symmetry.jl")
include("Toolkit.jl")

include("prettyprintln.jl")

end
