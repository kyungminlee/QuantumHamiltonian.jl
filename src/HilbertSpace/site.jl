export State, Site
export bitwidth, get_state, dimension
export qntype
export quantum_number_sectors
export get_quantum_number
export compress
export get_state_index


qntype(::T) where T = qntype(T)


"""
    State{QN}

State with quantum number type `QN`.

# Examples
```jldoctest
julia> using QuantumHamiltonian

julia> up = State("Up", 1)
State{Tuple{Int64}}("Up", (1,))

julia> State("Dn", (-1, 1))
State{Tuple{Int64, Int64}}("Dn", (-1, 1))
```
"""
struct State{QN<:Tuple{Vararg{<:AbstractQuantumNumber}}}
    name::String
    quantum_number::QN
    State(name::AbstractString) = new{Tuple{}}(name, ())

    function State(name::AbstractString, quantum_number::Integer)
        return new{Tuple{Int}}(name, (quantum_number,))
    end

    function State(name::AbstractString, quantum_number::QN) where {QN<:Tuple{Vararg{<:AbstractQuantumNumber}}}
        return new{QN}(name, quantum_number)
    end
end


function Base.:(==)(lhs::State{Q}, rhs::State{Q}) where Q
    return (lhs.name == rhs.name) && (lhs.quantum_number == rhs.quantum_number)
end


"""
    qntype(::Type{State{QN}})

Returns the quantum number type of the given state type.
"""
qntype(::Type{State{QN}}) where QN = QN


"""
    Site{QN}

A site with quantum number type `QN`.

# Examples
```jldoctest
julia> using QuantumHamiltonian

julia> site = Site([State("Up", 1), State("Dn", -1)]);
```
"""
struct Site{QN<:Tuple{Vararg{<:AbstractQuantumNumber}}}<:AbstractHilbertSpace{QN}
    states::Vector{State{QN}}
    Site(states::AbstractVector{State{QN}}) where QN = new{QN}(states)
end


"""
    qntype(::Type{Site{QN}})

Returns the quantum number type of the given site type.
"""
qntype(::Type{Site{QN}}) where QN = QN


function Base.:(==)(lhs::Site{Q1}, rhs::Site{Q2}) where {Q1, Q2}
    return lhs.states == rhs.states
end


"""
    bitwidth(site)

Number of bits necessary to represent the states of the given site.
"""
bitwidth(site::Site) = Int(ceil(log2(length(site.states))))


"""
    dimension(site)

Hilbert space dimension of a given site (= number of states).
"""
dimension(site::Site) = length(site.states)


"""
    get_state(site, binrep) where {QN, BR<:Unsigned}

Returns the state of `site` represented by the bits `binrep`.
"""
function get_state(site::Site, binrep::U) where {U<:Unsigned}
    return site.states[Int(binrep+1)]
end


"""
    compress(site, state_index, binary_type=UInt) -> binary_type

Get binary representation of the state specified by `state_index`.
Check bounds `1 <= state_index <= dimension(site)`, and returns binary representation of `state_index-1`.
"""
@inline function compress(site::Site, state_index::Integer, binary_type::Type{BR}=UInt) where {BR<:Unsigned}
    @boundscheck if !(1 <= state_index <= dimension(site))
        throw(BoundsError("attempt to access a $(dimension(site))-state site at index $state_index"))
    end
    return BR(state_index-1)
end


"""
    get_state_index(site, binrep)

Gets the state index of the binary representation. Returns `Int(binrep+1)`.
"""
@inline function get_state_index(site::Site, binrep::U) where {U<:Unsigned}
    i = Int(binrep+1)
    @boundscheck if !(1 <= i <= dimension(site))
        throw(BoundsError("attempt to access a $(dimension(site))-state site at index $i"))
    end
    return i
end


"""
    quantum_number_sectors(site) -> Vector{QN}

Gets a list of possible quantum numbers as a sorted vector of QN.
"""
function quantum_number_sectors(site::Site{QN})::Vector{QN} where QN
    return sort(collect(Set([state.quantum_number for state in site.states])))
end


"""
    get_quantum_number(site, state_index)

Gets the quantum number of state specified by state_index.
"""
function get_quantum_number(site::Site{QN}, state_index::Integer)::QN where QN
    return site.states[state_index].quantum_number
end


Base.keys(site::Site) = Base.OneTo(dimension(site))
