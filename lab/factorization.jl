using QuantumHamiltonian
# export get_quantum_numbers

"""
Represents Hilbert space broken up into sum of sectors of Abelian quantum numbers.
H = H₁ ⊕ H₂ ⊕ H₃ ⊕ H₄
Note that this also affects the "ordering" of the basis states.
"""
struct HilbertSpaceRepresentationSectorization{QN, HS, BR, DictType}
    hilbertspace::HS
    sectors::Dict{QN, HilbertSpaceRepresentation{HS, BR, DictType}}

    function HilbertSpaceRepresentationSectorization(
        hilbertspace::HS,
        sectors::AbstractDict{QN, <:HilbertSpaceRepresentation{HS, BR, DictType}}
    ) where {QN, HS, BR, DictType}
        return new{QN, HS, BR, DictType}(hilbertspace, sectors)
    end
end


get_quantum_numbers(s::HilbertSpaceRepresentationSectorization) = sort(collect(keys(s.sectors)))

function represent_sectors_dict(hs::HilbertSpace{QN}, ::Type{BR}=UInt) where {QN, BR<:Unsigned}
    HS = HilbertSpace{QN}
    DictType = Dict{BR, Int}
    
    sector_basis_list = Dict{QN, Vector{UInt}}()
    for indexarray in keys(hs)
        bvec = compress(hs, indexarray, UInt)
        qn = reduce((x,y) -> x.+y, get_quantum_number(site, i)
                    for (site, i) in zip(hs.sites, indexarray.I))
        if haskey(sector_basis_list, qn)
            push!(sector_basis_list[qn], bvec)
        else
            sector_basis_list[qn] = [bvec]
        end
    end
    sectors = Dict{QN, HilbertSpaceRepresentation{HS, BR, DictType}}()
    for (qn, basis_list) in sector_basis_list
        basis_lookup = DictType(b => i for (i, b) in enumerate(basis_list))
        sectors[qn] = HilbertSpaceRepresentation(hs, basis_list, basis_lookup)
    end
    return HilbertSpaceRepresentationSectorization(hs, sectors)
end


QuantumHamiltonian.get_quantum_number(s::HilbertSpaceRepresentationSectorization, bvec::Unsigned) = get_quantum_number(s.hilbertspace, bvec)
QuantumHamiltonian.bitwidth(h::HilbertSpaceRepresentationSectorization) = bitwidth(h.hilbertspace)

function QuantumHamiltonian.get_bitmask(hs::HilbertSpaceRepresentationSectorization, ::Type{BR}=UInt)::BR where {BR<:Unsigned}
    return make_bitmask(bitwidth(hs), BR)
end


function combine_quantum_number(quantum_number_lists::Vararg{AbstractVector{QN}, N}) where {QN, N}
    qn_factorization = Dict{QN, Vector{NTuple{N, QN}}}()
    for qns in Iterators.product(quantum_number_lists...)
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

# qnf = combine_quantum_number(collect(keys(ss1.sectors)), collect(keys(ss1.sectors)))


struct FactorizedHilbertSpaceRepresentation{QN, N, T<:Tuple{Vararg{HilbertSpaceRepresentationSectorization{QN}, N}}}
    subspaces::T #NTuple{N, HilbertSpaceRepresentationSectorization{QN, HS, BR, DictType}}
    quantum_number_list::Vector{NTuple{N, QN}}
    quantum_number_lookup::Dict{NTuple{N, QN}, Int}
    offsets::Vector{Int}
    function FactorizedHilbertSpaceRepresentation(
        subspaces::NTuple{N, HilbertSpaceRepresentationSectorization{QN}},
        quantum_number_list::Vector{NTuple{N, QN}},
        quantum_number_lookup::Dict{NTuple{N, QN}, Int},
        offsets::Vector{Int},
        # quantum_number_decompositions::Dict{QN, Vector{NTuple{N, QN}}},
        # offsets::Dict{QN, Vector{Int}},
    ) where {QN, N}
        T = typeof(subspaces)
        return new{QN, N, T}(subspaces, quantum_number_list, quantum_number_lookup, offsets)
    end
end



function QuantumHamiltonian.get_basis_index_amplitude(f::FactorizedHilbertSpaceRepresentation, bvec::BR) where {BR<:Unsigned}
    bitwidths = bitwidth.(f.subspaces)
    bvecs = breakup(bitwidths, bvec)
    qs = get_quantum_number.(f.subspaces, bvecs)
    if !haskey(f.quantum_number_lookup, qs)
        return (-1, 0)
    end
    i = f.quantum_number_lookup[qs]
    j = 0
    aout = 1
    for (si, qi, bi) in reverse(collect(zip(f.subspaces, qs, bvecs)))
        j *= dimension(si.sectors[qi])
        k, a = get_basis_index_amplitude(si.sectors[qi], bi)
        @assert isone(a)
        aout *= a
        j += (k-1)
    end
    return (f.offsets[i] + j + 1, aout)
end



# function HilbertSpaceRepresentation(
#     hilbert_space::AbstractHilbertSpace,
#     basis_list::AbstractVector{BR},
#     basis_lookup::DictType


# function get_basis_state(hsr::HilbertSpaceRepresentation, index::Integer)
    # return hsr.basis_list[index]
# end

# function get_basis_index_amplitude(hsr::HilbertSpaceRepresentation, bvec::Unsigned)
    # index = get(hsr.basis_lookup, bvec, -1)
    # return (index=index, amplitude=(index <= 0) ? 0 : 1)
# end


# struct ProductHilbertSpace{QN, T<:Tuple{Vararg{<:AbstractHilbertSpace{QN}}}}
#     subspaces::T
#     ProductHilbertSpace(subspaces::AbstractHilbertSpace{QN}...) where {QN} = new{QN, typeof(subspaces)}(subspaces)
# end
