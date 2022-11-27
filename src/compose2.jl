# compose with two data streams

#=

def online_covariance(data1, data2):
    meanx = meany = C = n = 0
    for x, y in zip(data1, data2):
        n += 1
        dx = x - meanx
        meanx += dx / n
        meany += (y - meany) / n
        C += dx * (y - meany)

    population_covar = C / n
    # Bessel's correction for sample variance
    sample_covar = C / (n - 1)
=#

mutable struct AccCov{T} <: Accumulator{T}
    nobs::Int
    xmean::T
    ymean::T
    c::T
end

function AccCov(::Type{T}=Float64) where {T}
    AccCov{T}(0, zero(T), zero(T), zero(T))
end

function (acc::AccCov{T})() where {T}
    acc.c / (acc.nobs - 1)
end

function (acc::AccCov{T})(x, y) where {T}
    acc.nobs += 1
    dx = x - acc.xmean
    dy = y - acc.ymean
    acc.xmean += dx / acc.nobs
    acc.ymean += dy / acc.nobs
    acc.c += dx * dy
    acc
end

function (acc::AccCov{T})(xs::Seq{T}, ys::Seq{T}) where {T}
    acc.nobs += length(xs)
    xmean = vmean(xs)
    dx = xmean - acc.xmean
    acc.xmean += dx / acc.nobs
    ymean = vmean(ys)
    dy = ymean - acc.ymean
    acc.ymean += dy / acc.nobs
    acc.c += dx * dy
    acc
end

