include("lab/factorization.jl")

up = State("Up", +1)
dn = State("Dn", -1)
site = Site([up, dn])

hs1 = HilbertSpace([site, site])
hsr1 = represent_dict(hs1)
ss1 = represent_sectors_dict(hs1)

# hs2 = HilbertSpace([site, site])
# hsr2 = represent_dict(hs2)


function testfunc3(subspaces::Vararg{HilbertSpaceRepresentationSectorization{QN}, N}) where {QN, N}
    target_quantum_numbers = QN[(0,)]
    qnd = combine_quantum_number(get_quantum_numbers.(subspaces)...)
    select_qtotal = sort(collect(qtotal for qtotal in keys(qnd)
                          if isempty(target_quantum_numbers) || (qtotal in target_quantum_numbers)))
    offsets = Int[]
    qni = NTuple{N, QN}[]
    off = 0
    for qtotal in select_qtotal
        qslist = qnd[qtotal]
        for qs in qslist
            push!(qni, qs)
            push!(offsets, off)
            d = prod(dimension(h.sectors[q]) for (h, q) in zip(subspaces, qs))
            off += d
        end
    end
    qnl = Dict(v => k for (k, v) in enumerate(qni))
    return FactorizedHilbertSpaceRepresentation(subspaces, qni, qnl, offsets)
end
foo = testfunc3(ss1, ss1, ss1)
bl = get_basis_list(foo)

function QuantumHamiltonian.get_basis_list(f::FactorizedHilbertSpaceRepresentation)
    basis_list = UInt[]
    for (qs, offset) in zip(f.quantum_number_list, f.offsets)
        for bs in Iterators.product([get_basis_list(subspace.sectors[q]) for (subspace, q) in zip(f.subspaces, qs)]...)
            b = zero(UInt)
            for ibi in length(bs):-1:1
                bi = bs[ibi]
                b <<= bitwidth(f.subspaces[ibi])
                b |= bi
            end
            push!(basis_list, b)
        end
    end
    return basis_list
end


function breakup(bitwidths::NTuple{N, Integer}, bvec::BR) where {N, BR<:Unsigned}
    out = BR[]
    for bw in bitwidths
        bm = make_bitmask(bw, BR)
        push!(out, bm & bvec)
        bvec >>= bw
    end
    return tuple(out...)
end

function breakup(f::FactorizedHilbertSpaceRepresentation, bvec::BR) where {BR<:Unsigned}
    bitwidths = bitwidth.(f.subspaces)
    bvecs = breakup(bitwidths, bvec)
    qs = get_quantum_number.(f.subspaces, bvecs)
    i = f.quantum_number_lookup[qs]
    j = 0
    # @show f.offsets[i]
    for (si, qi, bi) in reverse(collect(zip(f.subspaces, qs, bvecs)))
        j *= dimension(si.sectors[qi])
        k, a = get_basis_index_amplitude(si.sectors[qi], bi)
        @assert isone(a)
        j += (k-1)
        # @show k
    end
    return f.offsets[i] + j + 1
end

for b in bl
    i = breakup(foo, b)
    println(i, "\t", string(b, base=2, pad=6), "\t", string(bl[i], base=2, pad=6))
end


# struct ProductSpaceRep{QN}
#     subspaces::Vector{HilbertSpace{QN}}
# end

# function testfunc2(subspaces::Vararg{HilbertSpaceRepresentationSectorization{QN}, N}) where {QN, N}
#     qnd = combine_quantum_number(get_quantum_numbers.(subspaces)...)
#     offsets = Dict{QN, Vector{Int}}()
#     @show typeof(subspaces)
#     for (qtotal, qslist) in qnd
#         offsets_current = Int[]
#         off = 0
#         for qs in qslist
#             push!(offsets_current, off)
#             # x is a tuple of QNs
#             @show length(qs)
#             @show length(subspaces)
#             @assert N == length(qs) == length(subspaces)
#             d = prod(dimension(h.sectors[q]) for (h, q) in zip(subspaces, qs))
#             off += d
#         end
#         # @show qtotal
#         offsets[qtotal] = offsets_current
#     end
#     return FactorizedHilbertSpaceRepresentation(subspaces, qnd, offsets)
# end
# foo = testfunc2(ss1, ss1)


# -------------------------------

# # phs = ProductHilbertSpace(hs1, hs1, hs1)
# # phs |> factorize

# s = SectorizedHilbertSpaceRepresentation(hs1)
# sss = (s,s,s)

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
