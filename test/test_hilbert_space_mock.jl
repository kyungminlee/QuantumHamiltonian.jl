using Test
using QuantumHamiltonian
using LinearAlgebra

#=

Mock test.

N_SITES spin-half sites. The binary representations have opposite order to the state indices

| state index | state name | quantum number | binary representation |
|:-----------:|:----------:|:--------------:|:---------------------:|
|       1     |     up     |    ( 1,)       |           0b1         |
|       2     |    down    |    (-1,)       |           0b0         |

=#
@testset "Mock test" begin
    N_SITES = 2

    struct MockSite<:AbstractHilbertSpace{Tuple{Int}}
    end

    QuantumHamiltonian.bitwidth(s::MockSite) = 1
    QuantumHamiltonian.dimension(s::MockSite) = 2

    function QuantumHamiltonian.get_state(s::MockSite, binrep::Unsigned)
        if binrep == 0b1
            return State("up", (1,))
        elseif binrep == 0b0
            return State("dn", (-1,))
        else
            error("unsupported binrep $binrep")
        end
    end

    function QuantumHamiltonian.get_state_index(s::MockSite, binrep::Unsigned)
        if binrep == 0b1
            return 1
        elseif binrep == 0b0
            return 2
        else
            error("unsupported binrep $binrep")
        end
    end

    function QuantumHamiltonian.compress(s::MockSite, i::Integer, ::Type{BR}=UInt) where {BR<:Unsigned}
        (1 <= i <= 2) || error("i out of bounds")
        return i == 1 ? one(BR) : zero(BR)  # opposite order
    end

    # function QuantumHamiltonian.represent(s::MockSite, i::Integer, ::Type{BR}=UInt) where {BR<:Unsigned}
    #     (1 <= i <= 2) || error("i out of bounds")
    #     return i == 1 ? one(BR) : zero(BR)  # opposite order
    # end

    QuantumHamiltonian.quantum_number_sectors(s::MockSite) = [(-1,), (1,)]

    QuantumHamiltonian.get_quantum_number(s::MockSite, state_index::Integer) = state_index == 1 ? (1,) : (-1,)

    Base.keys(::MockSite) = Base.OneTo(2)

    struct MockHilbertSpace<:AbstractHilbertSpace{Tuple{Int}}
    end

    QuantumHamiltonian.basespace(h::MockHilbertSpace) = h
    QuantumHamiltonian.numsites(h::MockHilbertSpace) = N_SITES
    QuantumHamiltonian.get_site(h::MockHilbertSpace, i::Integer) = MockSite()

    QuantumHamiltonian.bitwidth(h::MockHilbertSpace) = N_SITES
    function QuantumHamiltonian.bitwidth(h::MockHilbertSpace, isite::Integer)
        @boundscheck (1 <= isite <= N_SITES) || error("isite out of bounds $isite")
        return 1
    end

    function QuantumHamiltonian.bitoffset(h::MockHilbertSpace, isite::Integer)
        @boundscheck (1 <= isite <= N_SITES) || error("isite out of bounds $isite")
        return (isite-1)
    end

    QuantumHamiltonian.get_quantum_numbers(h::MockHilbertSpace) = [(i,) for i in -N_SITES:2:N_SITES]

    QuantumHamiltonian.quantum_number_sectors(hs::MockHilbertSpace) = get_quantum_numbers(hs)

    function QuantumHamiltonian.get_quantum_number(h::MockHilbertSpace, binrep::Unsigned)
        return (N_SITES - 2 * count_ones(binrep),)
    end

    function QuantumHamiltonian.get_quantum_number(h::MockHilbertSpace, indexarray::AbstractVector{<:Integer})
        s = mapreduce(x -> x == 1 ? 1 : -1, +, indexarray)
        return (s,)
    end

    function QuantumHamiltonian.extract(h::MockHilbertSpace, binrep::Unsigned)
        out = Int[]
        for i in 1:N_SITES
            push!(out, 2 - Int(binrep & 0b1))
            binrep >>= 1
        end
        return CartesianIndex(out...)
    end

    QuantumHamiltonian.uncompress(h::MockHilbertSpace, binrep::Unsigned) = extract(h, binrep)

    function QuantumHamiltonian.compress(h::MockHilbertSpace, indexarray::CartesianIndex, ::Type{BR}=UInt) where {BR<:Unsigned}
        if length(indexarray) != N_SITES
            error("indexarray should be of length $N_SITES")
        end
        binrep = zero(BR)
        for i in reverse(indexarray.I)
            binrep = (binrep << 1) | (i == 1 ? 0b1 : 0b0)
        end
        return binrep
    end

    Base.keys(::MockHilbertSpace) = CartesianIndices(((2 for i in 1:N_SITES)...,))


    hs = MockHilbertSpace()
    
    sx(isite) = pure_operator(hs, isite, 1, 2, 1) +  pure_operator(hs, isite, 2, 1, 1)
    sy(isite) = pure_operator(hs, isite, 1, 2, -im) + pure_operator(hs, isite, 2, 1, im)
    sz(isite) = pure_operator(hs, isite, 1, 1, 1) + pure_operator(hs, isite, 2, 2, -1)

    s(mu, isite) = if mu == 1
        sx(isite)
    elseif mu == 2
        sy(isite)
    elseif mu == 3
        sz(isite)
    else
        error("mu $mu")
    end

    px = [0 1; 1 0]
    py = [0 -im; im 0]
    pz = [1 0; 0 -1]
    p = [px, py, pz]

    # check order of kron. site 1 should come at the end
    @test kron([1 2; 3 4], [1 1; 1 1]) == [
        1 1 2 2;
        1 1 2 2;
        3 3 4 4;
        3 3 4 4
    ]

    @testset "no quantum number" begin
        hsr = represent(hs)
        bl = get_basis_list(hsr)

        @test bl == [0x0, 0x1, 0x2, 0x3]
        @test [extract(hs, x) for x in bl] == [
            CartesianIndex(2,2), CartesianIndex(1,2),
            CartesianIndex(2,1), CartesianIndex(1,1),
        ]
        @test [compress(hs, extract(hs, x)) for x in bl] == bl

        a = [
            1.0 2.0 3.0;
            4.0 5.0 6.0;
            7.0 8.0 9.0
        ]
        hamiltonian = sum(
            s(i, 1) * a[i, j] * s(j, 2) for i in 1:3 for j in 1:3 
        )    
        hamiltonian_rep = represent(hsr, hamiltonian)
        m1 = Matrix(hamiltonian_rep)

        m2 = sum(
            a[i,j] * kron(p[i], p[j])
            for i in 1:3 for j in 1:3
        )
        m3 = sum(
            a[i,j] * kron(p[j], p[i])
            for i in 1:3 for j in 1:3
        )
        @test m1 != m2
        @test isapprox(eigvals(m1), eigvals(m2))
        @test m1 != m3
        @test isapprox(eigvals(m1), eigvals(m3))
    end

    @testset "with quantum number" begin
        hsr = represent(hs, [(0,)])
        bl = get_basis_list(hsr)
        @test bl == [0b01, 0b10]
        m1 = Matrix(represent(hsr, sx(1) * sx(2) + sy(1) * sy(2) + 0.7 * sz(1)))
        m2large = kron(px, px) + kron(py, py) + 0.7 * kron([1 0; 0 1], pz)
        m2 = m2large[2:3, 2:3]
        @test m1 != m2
        @test m1 == m2large[3:-1:2, 3:-1:2]
        @test eigvals(m1) == eigvals(m2)
    end
end