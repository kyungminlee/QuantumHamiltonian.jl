using QuantumHamiltonian

up = State("Up", +1)
dn = State("Dn", -1)
site = Site([up, dn])

# hs = HilbertSpace([site, site, site, site])
# hsr = represent_dict(hs)
# qns = quantum_number_sectors(hs)

hs1 = HilbertSpace([site, site])
hsr1 = represent_dict(hs1)
hs2 = HilbertSpace([site, site])
hsr2 = represent_dict(hs2)

struct ProductHilbertSpace{QN, T<:Tuple{Vararg{<:AbstractHilbertSpace{QN}}}}
    subspaces::T
    ProductHilbertSpace(subspaces::AbstractHilbertSpace{QN}...) where {QN} = new{QN, typeof(subspaces)}(subspaces)
end


function factorize(h::ProductHilbertSpace{QN, <:Tuple{Vararg{AbstractHilbertSpace{QN}, N}}}) where {QN, N}
    qnsets = quantum_number_sectors.(h.subspaces)  # ([qn1, qn2, ...], [qn3, qn4, ...], [qn5, qn6, ...])
    qn_factorization = Dict{QN, Vector{NTuple{N, QN}}}()
    for qns in Iterators.product(qnsets...)
        qn = reduce( (x,y) -> x .+ y, qns)
        if haskey(qn_factorization, qn)
            push!(qn_factorization[qn], qns)
        else
            qn_factorization[qn] = [qns]
        end
    end
    # return [k => sort(qn_factorization[k]) for k in sort(collect(keys(qn_factorization)))]
    for v in values(qn_factorization)
        sort!(v)
    end
    return qn_factorization
end

phs = ProductHilbertSpace(hs1, hs1, hs1)
phs |> factorize

"""
Represents Hilbert space broken up into sum of sectors of Abelian quantum numbers.
H = H₁ ⊕ H₂ ⊕ H₃ ⊕ H₄
Note that this also affects the "ordering" of the basis states.
"""
struct SectorizedHilbertSpaceRepresentation{QN}
    hilbertspace::HilbertSpace{QN}
    sectors::Dict{QN, Vector{UInt}}
    function SectorizedHilbertSpaceRepresentation(hs::HilbertSpace{QN}) where {QN}
        sectors = Dict{QN, Vector{UInt}}()
        for indexarray in keys(hs)
            bvec = compress(hs, indexarray, UInt)
            qn = reduce((x,y) -> x.+y, get_quantum_number(site, i) for (site, i) in zip(hs.sites, indexarray.I))
            if haskey(sectors, qn)
                push!(sectors[qn], bvec)
            else
                sectors[qn] = [bvec]
            end
        end
        return new{QN}(hs, sectors)
    end
end

s = SectorizedHilbertSpaceRepresentation(hs1)
sss = (s,s,s)

qnf = factorize(phs)
sector_decomposition = qnf[(0,)]
function get_basis_list(sector_decomposition)
    dim = 0
    offsets = Int[0]
    for qns in sector_decomposition
        d = prod(length(s.sectors[qn])
                for (s, qn) in zip(sss, qns))
        dim += d
        push!(offsets, dim)
    end
    return offsets
end
get_basis_list(sector_decomposition)


struct FactorizedHilbertSpaceRepresentation{QN, N}
    subspaces::NTuple{N, SectorizedHilbertSpaceRepresentation{QN}}
    quantum_number_decompositions::Dict{QN, Vector{NTuple{N, QN}}}
    offsets::Dict{QN, Vector{Int}}
    function FactorizedHilbertSpaceRepresentation(
        subspaces::NTuple{N, SectorizedHilbertSpaceRepresentation{QN}},
        quantum_number_decompositions::Dict{QN, Vector{NTuple{N, QN}}},
        offsets::Dict{QN, Vector{Int}},
    ) where {QN, N}
        return new{QN, N}(subspaces, quantum_number_decompositions, offsets)
    end
end


function testfunc2(s::Vararg{SectorizedHilbertSpaceRepresentation{QN}, N}) where {QN, N}
    qnd = Dict{QN, Vector{NTuple{N, QN}}}() # quantum number decompositions
    for qs in Iterators.product([keys(x.sectors) for x in s]...)
        q = reduce((x,y) -> x.+y, qs)
        if haskey(qnd, q)
            push!(qnd[q], qs)
        else
            qnd[q] = [qs]
        end
    end
    for l in values(qnd)
        sort!(l)
    end
    offsets = Dict{QN, Vector{Int}}()
    for (qtotal, qslist) in qnd
        offsets_current = Int[]
        off = 0
        for qs in qslist
            push!(offsets_current, off)
            # x is a tuple of QNs
            @assert N == length(qs) == length(s)
            d = prod(length(h.sectors[q]) for (h, q) in zip(s, qs))
            off += d
        end
        # @show qtotal
        offsets[qtotal] = offsets_current
    end
    return FactorizedHilbertSpaceRepresentation(s, qnd, offsets)
end

foo = testfunc2(s, s)

function decompose(hs, bvec::UInt)
    indexarray = uncompress(hs, bvec)
    p = 2 # TODO 
    i = indexarray.I[1:p]
    j = indexarray.I[p+1:end]
    bvec1 = compress(hs1, CartesianIndex(i...))
    bvec2 = compress(hs2, CartesianIndex(j...))
    
    q1 = reduce((x,y) -> x .+ y, get_quantum_number(site, k) for (site, k) in zip(hs1.sites, i))
    q2 = reduce((x,y) -> x .+ y, get_quantum_number(site, k) for (site, k) in zip(hs2.sites, j))
    # @show q1, q2
    q = q1 .+ q2
    idx = findfirst(x -> x == (q1, q2), foo[1][q])
    offset = foo[2][q][idx]
    basislist1 = s1.sectors[q1]
    basislist2 = s2.sectors[q2]

    # @show bvec1, bvec2
    # @show basislist1
    # @show basislist2
    idx1 = findfirst(x -> x == bvec1, basislist1)
    idx2 = findfirst(x -> x == bvec2, basislist2)
    # @show idx1, idx2
    n1, n2 = length(basislist1), length(basislist2)
    true_idx = offset + (idx2-1) * n1 + (idx1-1) + 1
    @show q, true_idx
end
decompose(hs, UInt(0b1000))


# function factorize(h::ProductHilbertSpace{QN, <:Tuple{Vararg{AbstractHilbertSpace{QN}, N}}}) where {QN, N}
#     qnsets = quantum_number_sectors.(h.subspaces)  # ([qn1, qn2, ...], [qn3, qn4, ...], [qn5, qn6, ...])
#     qn_factorization = Dict{QN, Vector{NTuple{N, QN}}}()
#     for qns in Iterators.product(qnsets...)
#         qn = reduce( (x,y) -> x .+ y, qns)
#         if haskey(qn_factorization, qn)
#             push!(qn_factorization[qn], qns)
#         else
#             qn_factorization[qn] = [qns]
#         end
#     end
#     # return [k => sort(qn_factorization[k]) for k in sort(collect(keys(qn_factorization)))]
#     for v in values(qn_factorization)
#         sort!(v)
#     end
#     return qn_factorization
# end


function decompose(hs, bvec::UInt)
    indexarray = uncompress(hs, bvec)
    p = 2 # TODO 
    i = indexarray.I[1:p]
    j = indexarray.I[p+1:end]
    bvec1 = compress(hs1, CartesianIndex(i...))
    bvec2 = compress(hs2, CartesianIndex(j...))
    
    q1 = reduce((x,y) -> x .+ y, get_quantum_number(site, k) for (site, k) in zip(hs1.sites, i))
    q2 = reduce((x,y) -> x .+ y, get_quantum_number(site, k) for (site, k) in zip(hs2.sites, j))
    # @show q1, q2
    q = q1 .+ q2
    idx = findfirst(x -> x == (q1, q2), foo[1][q])
    offset = foo[2][q][idx]
    basislist1 = s1.sectors[q1]
    basislist2 = s2.sectors[q2]

    # @show bvec1, bvec2
    # @show basislist1
    # @show basislist2
    idx1 = findfirst(x -> x == bvec1, basislist1)
    idx2 = findfirst(x -> x == bvec2, basislist2)
    # @show idx1, idx2
    n1, n2 = length(basislist1), length(basislist2)
    true_idx = offset + (idx2-1) * n1 + (idx1-1) + 1
    @show q, true_idx
end
decompose(hs, UInt(0b1000))


bl = hsr.basis_list
for b in bl
    decompose(hs, b)
end


# function hs_get_basis_list(hs::HilbertSpace{QN}, binary_type::Type{BR}=UInt)::Vector{BR} where {QN, BR<:Unsigned}
#     if sizeof(BR) * 8 <= bitwidth(hs)
#         throw(ArgumentError("type $(BR) not enough to represent the hilbert space (need $(bitwidth(hs)) bits)"))
#     end
#     basis_list = BR[]
#     for indexarray in keys(hs)
#         push!(basis_list, compress(hs, indexarray, BR))
#     end
#     sort!(basis_list)
#     return basis_list
# end
