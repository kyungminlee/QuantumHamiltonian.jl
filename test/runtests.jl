using Test
using QuantumHamiltonian

include("test_util.jl")
include("test_frozensortedarray.jl")
include("test_indexedvector.jl")

include("test_hilbert_space.jl")
include("test_hilbert_space_sector.jl")
include("test_product_hilbert_space.jl")

include("test_operator.jl")

include("test_sparse_state.jl")
include("test_operator_simplify.jl")
include("test_operator_application.jl")

include("test_hilbert_space_representation.jl")
include("test_product_hilbert_space_representation.jl")
include("test_decomposed_hilbert_space_representation.jl")
include("test_decompose.jl")

include("test_operator_representation.jl")
include("test_state_representation.jl")

include("test_symmetry_apply.jl")
include("test_symmetry_reduce.jl")
include("test_reduced_representation.jl")

include("test_bitflip.jl")
include("test_lgp.jl")
include("test_semidirectproductoperation.jl")

include("test_prettyprintln.jl")

include("test_toolkit.jl")
