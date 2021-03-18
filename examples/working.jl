using QuantumHamiltonian

QN = Int
Scalar = ComplexF64

up = State{QN}("Up", 1)
dn = State{QN}("Dn",-1)
spinsite = Site{QN}([up, dn])

n_sites = 4


n1 = 4
n2 = 4
n_sites = n1 * n2
hs = HilbertSpace{QN}([spin_site for i in 1:n_sites])

PAULI_MATRICES = [ Float64[0 1.0; 1.0 0.0], ComplexF64[0.0 -1.0*im; 1.0*im 0.0], Float64[1.0 0.0; 0.0 -1.0]]

sigma(i::Integer, j::Integer) = KroneckerProductOperator{ComplexF64}(hs, 1.0, Dict(i=>PAULI_MATRICES[j]))
sigma_plus(i::Integer) = KroneckerProductOperator{ComplexF64}(hs, 1.0, Dict(i=>[0.0 1.0; 0.0 0.0]))
sigma_minus(i::Integer) = KroneckerProductOperator{ComplexF64}(hs, 1.0, Dict(i=>[0.0 0.0; 1.0 0.0]))


kpo = sigma(3,1) * sigma(1,3)

s = convert(SumOperator{ComplexF64, UInt}, kpo)

ψ1 = SparseState{Float64, UInt}(hs, UInt(0x0))
ψ2 = SparseState{Float64, UInt}(hs, UInt(0x1))



@show ψ
@show sigma(1, 1)
@show ψ * sigma(1, 1)

H = sigma(2, 1)
H2 = GenericOperator{Scalar}(hs)

push!(H2, H)
push!(H2, H)
