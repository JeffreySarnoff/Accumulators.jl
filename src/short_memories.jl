#=
    Enable `lookback K steps | tiles` to utilize 1..K prior observations or 1..K prior updates
=#


# Count
mutable struct AccCount{T} <: Accumulator{T}
    n⬭::Int

    function AccCount(n::T=0) where {T}
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
    obs⬭::Tuple{N,T}

    function AccObs(n::Int=0, obs::NTuple{N,T}=ntuple(∅->0, N)) where {T, N}
        new{T, N}(n, obs)
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
