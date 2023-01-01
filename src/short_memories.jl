#=
    Enable `lookback K steps | tiles` to utilize 1..K prior observations or 1..K prior updates
=#


# Count
mutable struct AccumCount{T,F} <: AccumCacher{T,F}
    nᵪ::Int

    function AccumCount(; n::T=0) where {T,F}
        new{T,F}(n)
    end
end

function (accum::AccumCount{T,F})() where {T,F}
    accum.nᵪ
end

function (accum::AccumCount{T,F})(x) where {T,F}
    accum.nᵪ += 1
    accum
end

function (accum::AccumCount{T,F})(xs::Seq) where {T,F}
    accum.nᵪ += length(xs)
    accum
end


# Short Memory, accumumulate observations, spilling oldest once filled
mutable struct AccumObs{T,F,N} <: AccumCacher{T,F}
    nᵪ::Int
    obsᵪ::NTuple{N,T}
end

function AccumObs{T,F}(nobs::Int=1) where {T,F}
    init_obs = ntuple(∅ -> zero(T), nobs)
    AccumObs(0, init_obs)
end

function (accum::AccumObs{T,F,N})() where {T,F,N}
    accum.obsᵪ
end

function (accum::AccumObs{T,2})(x) where {T,F}
    accum.nᵪ += 1
    accum.obsᵪ = (accum.obsᵪ[2], x)
    accum
end

function (accum::AccumObs{T,F,N})(x) where {T,F,N}
    accum.nᵪ += 1
    accum.obsᵪ = (accum.obsᵪ[1:N-1]..., x)
    accum
end

function (accum::AccumObs{T,F,N})(xs::Seq) where {T,F,N}
    accum.nᵪ += length(xs)
    accum
end


# Minimum

mutable struct AccumMinimum1{T,F} <: AccumCacher{T,F}
    obsᵪ::NTuple{N,T}   # cache most recent N obs
    const fn::F         # preapply to each element when observed
    nobsᵪ::Int          # count each observation
    nmin::Int           # count distinct minima
    min::T              # current minimum
end

function AccumMinimum1(::Type{T}=AccumNum; fn::F=identity, ncached=4) where {T,F}
    cache = ntuple(i -> zero(T), ncached)
    AccMinimum1{T,F}(cache, 0, 0, typemax(T), fn)
end

function (acc::AccMinimum{T,F})() where {T,F}
    acc.min
end

function (acc::AccMinimum{T,F})(x::T) where {T,F}
    xx = acc.fn(x)
    acc.nobs += 1
    if x < acc.min
        acc.nmin += 1
        acc.min = xx
    end
    accum
end

function (acc::AccMinimum{T,F})(xs::Seq{T}) where {T,F}
    xxs = map(acc.fn, xs)
    acc.nobs += length(xs)
    x = vminimum(xxs)
    if x < acc.min
        acc.nmin += 1
        acc.min = x
    end
    accum
end
