#=
    Enable `lookback K steps | tiles` to utilize 1..K prior observations or 1..K prior updates
=#


# Count
mutable struct AccCount{T} <: Accumulator{T}
    n⬭::Int

    function AccCount(; n::T=0) where {T}
        new{T}(n)
    end
end

function (acc::AccCount{T})() where {T}
    acc.n⬭
end

function (acc::AccCount{T})(x) where {T}
    acc.n⬭ += 1
    acc
end

function (acc::AccCount{T})(xs::Seq) where {T}
    acc.n⬭ += length(xs)
    acc
end


# Short Memory, accumulate observations, spilling oldest once filled
mutable struct AccObs{T, N} <: Accumulator{T}
    n⬭::Int
    obs⬭::NTuple{N,T}
end

function AccObs{T}(nobs::Int=1) where {T}
     init_obs = ntuple(∅->zero(T), nobs)
     AccObs(0, init_obs)
end
     
function (acc::AccObs{T, N})() where {T, N}
    acc.obs⬭
end

function (acc::AccObs{T, 2})(x) where {T}
    acc.n⬭ += 1
    acc.obs⬭ = (acc.obs⬭[2], x)
    acc
end

function (acc::AccObs{T, N})(x) where {T, N}
    acc.n⬭ += 1
    acc.obs⬭ = (acc.obs⬭[1:N-1]..., x)
    acc
end

function (acc::AccObs{T, N})(xs::Seq) where {T, N}
    acc.n⬭ += length(xs)
    acc
end
