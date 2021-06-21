
using QuantumHamiltonian

function subspace(hs::HilbertSpace{QN}, isites) where {QN}
    return HilbertSpace([hs.sites[isite] for isite in isites])
end

nsites = 16
hs, spin = QuantumHamiltonian.Toolkit.spin_system(nsites, 1//2)

hs1 = subspace(hs, 1:8)
hs2 = hs1
# hs2 = subspace(hs, 5:8)

reps1 = Dict(qn1 => represent(HilbertSpaceSector(hs1, qn1)) for qn1 in qns1)
# reps2 = Dict(qn2 => represent(HilbertSpaceSector(hs2, qn2)) for qn2 in qns2)
reps2 = reps1

qns1 = quantum_number_sectors(hs1)
qns2 = quantum_number_sectors(hs2)

out = Dict()
for  qn2 in qns2, qn1 in qns1
    qn = qn1 .+ qn2
    if !haskey(out, qn)
        out[qn] = []
    end
    push!(out[qn], (qn1, qn2))
end

newout = Dict(qn => Dict() for qn in keys(out))

for (qn, qnsplit) in out
    lowerdim = 0
    for (qn1, qn2) in qnsplit
        d1 = haskey(reps1, qn1) ? dimension(reps1[qn1]) : 0
        d2 = haskey(reps2, qn2) ? dimension(reps2[qn2]) : 0
        newout[qn][qn2] = (qn1=>d1, qn2=>d2, lowerdim)
        lowerdim += d1 * d2
    end
end

# qn_split = [(qn, sort(out[qn])) for qn in sort(collect(keys(out)))]


basis_list = UInt[]
for (qn1, qn2) in out[(2,)]
    hsr1 = reps1[qn1]
    hsr2 = reps2[qn2]
    for bvec2 in hsr2.basis, bvec1 in hsr1.basis
        bw1 = bitwidth(hsr1)
        bvec = bvec2 << bw1 | bvec1
        push!(basis_list, bvec)
    end
end
bm1 = get_bitmask(hs1)
bm2 = get_bitmask(hs2) << bitwidth(hs1)

#=
for (idx0, bvec) in enumerate(basis_list)
    bvec1 = bm1 & bvec
    bvec2 = (bm2 & bvec) >> bitwidth(hs1)
    qn = get_quantum_number(hs, bvec)
    qn1 = get_quantum_number(hs1, bvec1)
    qn2 = get_quantum_number(hs2, bvec2)
    # @show qn, qn1, qn2, bvec, bvec1, bvec2

    ((_, dim1), (_, dim2), index_offset) = newout[qn][qn2]
    idx = index_offset + reps1[qn1].basis_lookup[bvec1] + (reps2[qn2].basis_lookup[bvec2]-1) * dim1
    @show bvec, idx0, idx, idx0 == idx
end
=#


function join_by_quantum_number_sectors(
    hs1::HilbertSpace{QN}, hs2::HilbertSpace{QN},
    reps1, reps2
) where {QN}
    qns1 = quantum_number_sectors(hs1)
    qns2 = quantum_number_sectors(hs2)

    qn_splits = Dict{QN, Vector{Tuple{QN, QN}}}()
    for  qn2 in qns2, qn1 in qns1
        qn = qn1 .+ qn2
        if !haskey(qn_splits, qn)
            qn_splits[qn] = Tuple{QN, QN}[]
        end
        push!(qn_splits[qn], (qn1, qn2))
    end
    
    ItemType = NamedTuple{(:qn1, :qn2, :dimension1, :dimension2, :dimension_offset), Tuple{QN, QN, Int, Int, Int}}
    hs_splits = Dict(qn => Dict{QN, ItemType}() for qn in keys(qn_splits))
    for (qn, qnsplit) in qn_splits
        doff = 0
        for (qn1, qn2) in qnsplit
            d1 = haskey(reps1, qn1) ? dimension(reps1[qn1]) : 0
            d2 = haskey(reps2, qn2) ? dimension(reps2[qn2]) : 0
            hs_splits[qn][qn2] = (qn1=qn1, qn2=qn2, dimension1=d1, dimension2=d2, dimension_offset=doff)
            doff += d1 * d2
        end
    end
    return hs_splits
end

hs_splits = join_by_quantum_number_sectors(hs1, hs2, reps1, reps2)

function basis_lookup(hs, hs1, hs2, hs_splits, bvec)

    bvec1 = bm1 & bvec
    bvec2 = (bm2 & bvec) >> bitwidth(hs1)
    qn = get_quantum_number(hs, bvec)
    qn1 = get_quantum_number(hs1, bvec1)
    qn2 = get_quantum_number(hs2, bvec2)
    
    ((_, dim1), (_, dim2), index_offset) = hs_split[qn2]
    idx = index_offset + reps1[qn1].basis_lookup[bvec1] + (reps2[qn2].basis_lookup[bvec2]-1) * dim1
    @show bvec, idx0, idx, idx0 == idx
end

struct ConjoinedHilbertSpace{QN, LHS<:AbstractHilbertSpace{QN}, HHS<:AbstractHilbertSpace{QN}}
    lower::LHS
    higher::HHS
    
    function ConjoinedHilbertSpace(lower::LHS, higher::HHS) where {QN, LHS<:AbstractHilbertSpace{QN}, HHS<:AbstractHilbertSpace{QN}}
        return new{QN, LHS, HHS}(lower, higher)
    end
end

function represent(hs::ConjoinedHilbertSpace)
    qns1 = quantum_number_sectors(hs.lower)
    qns2 = quantum_number_sectors(hs.higher)
    qnsplits = Dict{QN, Vector{Tuple{QN, QN}}}()
    for  qn2 in qns2, qn1 in qns1
        qn = qn1 .+ qn2
        if !haskey(qnsplits, qn)
            qnsplits[qn] = Tuple{QN, QN}[]
        end
        push!(qnsplits[qn], (qn1, qn2))
    end
    ItemType = NamedTuple{(:qn1, :qn2, :dimension1, :dimension2, :dimension_offset), Tuple{QN, QN, Int, Int, Int}}
    hssplits = Dict(qn => Dict{QN, ItemType}() for qn in keys(qn_splits))
    for (qn, qnsplit) in qn_splits
        doff = 0
        for (qn1, qn2) in qnsplit
            d1 = haskey(reps1, qn1) ? dimension(reps1[qn1]) : 0
            d2 = haskey(reps2, qn2) ? dimension(reps2[qn2]) : 0
            hssplits[qn][qn2] = (qn1=qn1, qn2=qn2, dimension1=d1, dimension2=d2, dimension_offset=doff)
            doff += d1 * d2
        end
    end
end