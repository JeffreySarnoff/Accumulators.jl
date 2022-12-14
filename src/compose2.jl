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

# Cov
# see https://softwareengineering.stackexchange.com/questions/337617/algorithm-for-windowed-online-covariance

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
    ifelse(acc.nobs < 2, zero(T), acc.c / (acc.nobs - 1))
end

function (acc::AccCov{T})(x::T, y::T) where {T}
    acc.nobs += 1
    dx = x - acc.xmean
    dy = y - acc.ymean
    acc.xmean += dx / acc.nobs
    acc.ymean += dy / acc.nobs
    acc.c += ifelse(acc.nobs === 1, zero(T), dx * dy)
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
    acc.c += ifelse(acc.nobs === 1, zero(T), dx * dy)
    acc
end

# Cov

mutable struct AccCov{T} <: Accumulator{T}
    nobs::Int
    xmean::T
    ymean::T
    scov::T
end

function AccCov(::Type{T}=Float64) where {T}
    AccCov{T}(0, zero(T), zero(T), zero(T))
end

function (acc::AccCov{T})() where {T}
    ifelse(acc.nobs < 2, zero(T), acc.scov)
end

function (acc::AccCov{T})(x::T, y::T) where {T}
    acc.nobs += 1
    dx = x - acc.xmean
    dy = y - acc.ymean
    acc.xmean += dx / acc.nobs
    acc.ymean += dy / acc.nobs
    t = inv(acc.nobs)
    tm1 = inv(t - 1)
    tm2 = acc.nobs - 2
    prior_scov = acc.scov
    acc.scov = (tm2 * tm1) * prior_scov + t * dx * dy
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
    acc.c += ifelse(acc.nobs === 1, zero(T), dx * dy)
    acc
end


# Cov
# see https://softwareengineering.stackexchange.com/questions/337617/algorithm-for-windowed-online-covariance

mutable struct AccCov{T} <: Accumulator{T}
    nobs::Int
    x::T
    y::T
    dx::T
    dy::T
    xsum::T
    ysum::T
    xmean::T
    ymean::T
    dxdy::T
    xy::T
    c::T
end

function AccCov(::Type{T}=Float64) where {T}
    AccCov{T}(0, zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T))
end

function (acc::AccCov{T})() where {T}
    ifelse(acc.nobs < 2, zero(T), acc.c / (acc.nobs - 1))
end

function (acc::AccCov{T})(x::T, y::T) where {T}
    acc.nobs += 1
    acc.xsum += x
    acc.ysum += y
    acc.dx = ifelse(acc.nobs>1, x - acc.x, zero(T))
    acc.dy = ifelse(acc.nobs>1, y - acc.y, zero(T))
    acc.x = x
    acc.y = y
    acc.xmean += acc.dx / acc.nobs
    acc.ymean += acc.dy / acc.nobs
    acc.dxdy += acc.dx * acc.dy
    acc.dxdy = ifelse(!isfinite(acc.dxdy), zero(T), acc.dxdy)
    acc.xy += x * y
    acc.xy = ifelse(!isfinite(acc.xy), zero(T), acc.xy)
    acc.c += ifelse(acc.nobs === 1, zero(T), acc.dxdy)
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
    acc.c += ifelse(acc.nobs === 1, zero(T), dx * dy)
    acc
end

# Corr

mutable struct AccCor{T} <: Accumulator{T}
    nobs::Int
    xmean::T
    ymean::T
    xsvar::T
    ysvar::T
    c::T
end

function AccCor(::Type{T}=Float64) where {T}
    AccCor{T}(0, zero(T), zero(T), zero(T), zero(T), zero(T))
end

function (acc::AccCor{T})() where {T}
    if acc.nobs < 2
        return zero(T)
    end
    corr = acc.c / (acc.nobs - 1)
    unbiased_xvar = acc.xsvar / (acc.nobs - 1)
    unbiased_yvar = acc.ysvar / (acc.nobs - 1)
    corr / (unbiased_xvar * unbiased_yvar)
end

function (acc::AccCor{T})(x::T, y::T) where {T}
    acc.nobs += 1
    prior_xmean = acc.xmean
    dx = x - prior_xmean
    acc.xmean += dx / acc.nobs
    acc.xsvar = acc.xsvar + (x - prior_xmean) * (x - acc.xmean)
    prior_ymean = acc.ymean
    dy = y - acc.ymean
    acc.ymean += dy / acc.nobs
    acc.ysvar = acc.ysvar + (y - prior_ymean) * (y - acc.ymean)
    acc.c += dx * dy
    acc
end

function (acc::AccCorr{T})(xs::Seq{T}, ys::Seq{T}) where {T}
    if length(xs) !== length(ys) 
        throw(ArgumentError(string("lengths must be equal (", n, " != ", length(ys), ")")))
    end
    acc.nobs += length(xs)
    prior_xmean = acc.xmean
    xmean = vmean(xs)
    dx = xmean - acc.xmean
    acc.xmean += dx / acc.nobs
    acc.xsvar = acc.xsvar + (xmean - prior_xmean) * (xmean - acc.xmean)
    prior_ymean = acc.ymean
    ymean = vmean(ys)
    dy = ymean - acc.ymean
    acc.ymean += dy / acc.nobs
    acc.ysvar = acc.ysvar + (ymean - prior_ymean) * (ymean - acc.ymean)
    acc.c += dx * dy
    acc
end

