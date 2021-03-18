using Random
using StatsBase

"""
    estimate_density([rng], m; nsample=10, tol=√ϵ)

Estimate the density of nonzero elements in the abstract operator representation.
Returns mean and error of the density.
`nsample` should be greater than 2.
If the dimension of the operator representation is less than 10, exact counting of nonzero elements is performed.
"""
function estimate_density(
    rng::AbstractRNG,
    m::AbstractOperatorRepresentation{T};
    nsample::Integer=min(10, size(m, 2)),
    tol::Real=Base.rtoldefault(real(T))
) where {T}
    dim = size(m, 2)
    vals = Dict{Int, T}()
    function column_density(col::Integer)
        empty!(vals)
        for (i, v) in get_column_iterator(m, col)
            vals[i] = get(vals, i, zero(T)) + v
        end
        choptol!(vals, tol)
        return length(vals) / dim
    end

    if dim < 10
        md = Matrix(m)
        μ = mean(x -> !iszero(x), md)
        return μ, 0.0
    else
        nsample >= 2 || throw(ArgumentError("number of samples should be no less than 2"))
        sample_column_list = sample(rng, 1:dim, nsample; replace=false)
        μ, σ = mean_and_std(column_density(col) for col in sample_column_list)
        return μ, σ * sqrt((dim - nsample) / ((dim - 1) * nsample))
    end
end


function estimate_density(
    m::AbstractOperatorRepresentation{T};
    nsample::Integer=min(10, size(m, 2)),
    tol::Real=Base.rtoldefault(real(T))
) where {T}
    return estimate_density(Random.GLOBAL_RNG, m; nsample=nsample, tol=tol)
end
