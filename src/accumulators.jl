#=
     AccCount, 
     AccMinimum, AccMaximum, AccExtrema, 
     AccSum, AccProd,
     AccMean, AccGeoMean, AccHarmMean, AccGenMean,
     AccMeanAndVar, AccMeanAndStd, AccStats,
     AccExpWtMean, AccExpWtMeanVar, AccExpWtMeanStd
=#

# Count

mutable struct AccCount{T,F} <: Accumulator{T,F}
    nobs::Int       # count each observation
    fn::F           # preapply to each element when observed
end

function AccCount(::Type{T}=Int64) where {T,F}
     AccCount{T,F}(zero(T))
end

function (acc::AccCount{T,F})() where {T,F}
    acc.nobs
end
     
function (acc::AccCount{T,F})(x::T) where {T,F}
    acc.nobs += one(I)
    accum
end

function (acc::AccCount{T,T})(xs::Seq{T}) where {T,F}
    acc.nobs += length(xs)
    accum
end

# Minimum

mutable struct AccMinimum{T,F} <: Accumulator{T,F}
    nobs::Int       # count each observation
    nmin::Int       # count distinct minima
    min::T          # current minimum
    fn::F           # preapply to each element when observed
end

function AccMinimum(::Type{T}=AccumNum; fn::F=identity) where {T,F}
     AccMinimum{T,F}(0, 0, typemax(T), fn)
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

# Maximum

mutable struct AccMaximum{T,F} <: Accumulator{T,F}
    nobs::Int       # count each observation
    nmax::Int       # count distinct maxima
    max::T          # current maximum
    fn::F
end

function AccMaximum(::Type{T}=AccumNum; fn::F=identity) where {T,F}
     AccMaximum{T,F}(0, 0, typemin(T), fn)
end

function (acc::AccMaximum{T,F})() where {T,F}
    acc.max
end

function (acc::AccMaximum{T,F})(x::T) where {T,F}
    xx = acc.fn(x)
    acc.nobs += 1
    if xx > acc.max
        acc.nmax += 1
        acc.max = xx
    end
    accum
end

function (acc::AccMaximum{T,F})(xs::Seq{T}) where {T,F}
    xxs = map(acc.fn, xs)
    acc.nobs += length(xs)
    x = vmaximum(xxs)
    if x > acc.max
        acc.nmax += 1
        acc.max = x
    end
    accum
end

# Extrema

mutable struct AccExtrema{T,F} <: Accumulator{T,F}
    nobs::Int       # count each observation
    nmin::Int       # count distinct minima
    nmax::Int       # count distinct maxima
    min::T          # current minimum
    max::T          # current maximum
    fn::F
end

function AccExtrema(::Type{T}=AccumNum; fn::F=identity) where {T,F}
     AccExtrema{T,F}(0, 0, 0, typemax(T), typemin(T), fn)
end

function (acc::AccExtrema{T,F})() where {T,F}
    (acc.min, acc.max)
end

function (acc::AccExtrema{T,F})(x::T) where {T,F}
    xx = acc.fn(x)
    acc.nobs += 1
    if xx < acc.min
        acc.nmin += 1
        acc.min = xx
    end
    if xx > acc.max
        acc.nmax += 1
        acc.max = xx
    end
    accum
end

function (acc::AccExtrema{T,F})(xs::Seq{T}) where {T,F}
    xxs = map(acc.fn, xs)
    acc.nobs += length(xs)     
    mn, mx = vextrema(xxs)
    if mn < acc.min
        acc.nmin += 1
        acc.min = mn
    end
    if mx > acc.max
        acc.nmax += 1
        acc.max = mx
    end
    accum
end

# Sum

mutable struct AccSum{T,F} <: Accumulator{T,F}
    nobs::Int       # count each observation
    sum::T          # current sum
    fn::F
end

function AccSum(::Type{T}=AccumNum; fn::F=identity) where {T,F}
     AccSum{T,F}(0, zero(T), fn)
end

function (acc::AccSum{T,F})() where {T,F}
    acc.sum
end

function (acc::AccSum{T,F})(x::T) where {T,F}
    xx = acc.fn(x)
    acc.nobs += 1
    acc.sum += xx
    accum
end

function (acc::AccSum{T,F})(xs::Seq{T}) where {T,F}
    xxs = map(acc.fn, xs)
    acc.nobs += length(xs)
    x = vsum(xxs)
    acc.sum += x
    accum
end

# Prod

mutable struct AccProd{T,F} <: Accumulator{T,F}
    nobs::Int       # count each observation
    prod::T         # current product
    fn::F
end

function AccProd(::Type{T}=AccumNum; fn::F=identity) where {T,F}
     AccProd{T,F}(0, one(T), fn)
end

function (acc::AccProd{T,F})() where {T,F}
    acc.prod
end

function (acc::AccProd{T,F})(x::T) where {T,F}
    xx = acc.fn(x)
    acc.nobs += 1
    acc.prod *= xx
    accum
end

function (acc::AccProd{T,F})(xs::Seq{T}) where {T,F}
    xxs = map(acc.fn, xs)
    acc.nobs += length(xs)
    x = vprod(xxs)
    acc.prod *= x
    accum
end

# Mean

mutable struct AccMean{T,F} <: Accumulator{T,F}
    nobs::Int       # count each observation
    mean::T         # current mean
    fn::F
end

function AccMean(::Type{T}=AccumNum; fn::F=identity) where {T,F}
     AccMean{T,F}(0, zero(T), fn)
end

function (acc::AccMean{T,F})() where {T,F}
    acc.mean
end

function (acc::AccMean{T,F})(x::T) where {T,F}
    xx = acc.fn(x)
    acc.nobs += 1
    acc.mean += (xx - acc.mean) / acc.nobs
    accum
end

function (acc::AccMean{T,F})(xs::Seq{T}) where {T,F}
    xxs = map(acc.fn, xs)
    acc.nobs += length(xs)
    xmean = vmean(xxs)
    acc.mean += (xmean - acc.mean) / acc.nobs
    accum
end

# GeoMean

mutable struct AccGeoMean{T,F} <: Accumulator{T,F}
    nobs::Int       # count each observation
    sumlog::T       # ∑(i=1:nobs) log(xᵢ)
    fn::F
end

function AccGeoMean(::Type{T}=AccumNum; fn::F=identity) where {T,F}
     AccGeoMean{T,F}(0, zero(T), fn)
end

function (acc::AccGeoMean{T,F})() where {T,F}
    n = ifelse(iszero(acc.nobs), 1, acc.nobs)
    exp(acc.sumlog / n)
end

function (acc::AccGeoMean{T,F})(x::T) where {T,F}
    xx = acc.fn(x)
    acc.nobs += 1
    acc.sumlog += logabs(xx)
    accum
end

function (acc::AccGeoMean{T,F})(xs::Seq{T}) where {T,F}
    xxs = map(acc.fn, xs)
    acc.nobs += length(xs)
    acc.sumlog += sum(map(logabs, xxs))
    accum
end

# Harmonic Mean
mutable struct AccHarmMean{T,F} <: Accumulator{T,F}
    nobs::Int       # count each observation
    invhmean::T     # 1 / current harmonic mean
    fn::F
end

function AccHarmMean(::Type{T}=AccumNum; fn::F=identity) where {T,F}
     AccHarmMean{T,F}(0, zero(T), fn)
end

function (acc::AccHarmMean{T,F})() where {T,F}
    n = ifelse(iszero(acc.nobs), 1, acc.nobs)
    n / acc.invhmean
end

function (acc::AccHarmMean{T,F})(x::T) where {T,F}
    xx = acc.fn(x)
    acc.nobs += 1
    acc.invhmean += one(T) / xx
    accum
end

function (acc::AccHarmMean{T,F})(xs::Seq{T}) where {T,F}
    xxs = map(acc.fn, xs)
    acc.nobs += length(xs)
    acc.invhmean += sum(map(inv, xxs))
    accum
end

# Generalized Mean (defaults to Quadratic Mean [root-mean-squared])

mutable struct AccGenMean{T,F} <: Accumulator{T,F}
    nobs::Int       # count each observation
    gmean::T        # current mean((xᵢ)ᵖʷʳ), i=1:nobs
    const pwr::T    # power
    const rpwr::T   # reciprocal of power
    fn::F
end

function AccGenMean(::Type{T}=AccumNum; fn::F=identity, power::Real=2.0) where {T,F}
    AccGenMean{T,F}(0, zero(T), T(power), T(1/power), fn)
end

function (acc::AccGenMean{T,F})() where {T,F}
    (acc.gmean)^acc.rpwr
end

function (acc::AccGenMean{T,F})(x::T) where {T,F}
    xx = acc.fn(x)
    acc.nobs += 1
    acc.gmean += ((xx^acc.pwr) - acc.gmean) / acc.nobs
    accum
end

function (acc::AccGenMean{T,F})(xs::Seq{T}) where {T,F}
    xxs = map(acc.fn, xs)
    acc.nobs += length(xs)     
    xmean = vmean(map(x->x^acc.pwr, xxs))
    acc.gmean += (xmean - acc.gmean) / acc.nobs
    accum
end

# Unbiased Sample Variation (with Mean)
# see https://www.johndcook.com/blog/standard_deviation/

mutable struct AccMeanAndVar{T,F} <: Accumulator{T,F}
    nobs::Int       # count each observation
    mean::T         # current mean
    svar::T         # sum of variances x[1:1=0,1:2,..,1:nobs]
    fn::F
end

function AccMeanAndVar(::Type{T}=Float64; fn::F=identity) where {T,F}
    AccMeanAndVar{T,F}(0, zero(T), zero(T), fn)
end

function (acc::AccMeanAndVar{T,F})() where {T,F}
    unbiased_var = acc.svar / (acc.nobs - 1)
    (mean=acc.mean, var=unbiased_var)
end

function (acc::AccMeanAndVar{T,F})(x::T) where {T,F}
    xx = acc.fn(x)
    acc.nobs += 1
    prior_mean = acc.mean
    acc.mean = prior_mean + (xx - prior_mean) / acc.nobs
    acc.svar = acc.svar + (xx - prior_mean) * (xx - acc.mean)
    accum
end

function (acc::AccMeanAndVar{T,F})(xs::Seq{T}) where {T,F}
    xxs = map(acc.fn, xs)
    acc.nobs += length(xs)
    prior_mean = acc.mean
    xmean = vmean(xxs)
    acc.mean += (xmean - prior_mean) / acc.nobs
    acc.svar = acc.svar + (xmean - prior_mean) * (xmean - acc.mean)
    accum
end

mutable struct AccMeanAndStd{T,F} <: Accumulator{T,F}
    nobs::Int       # count each observation
    mean::T         # current mean
    svar::T         # sum of variances x[1:1=0,1:2,..,1:nobs]
    fn::F
end

function AccMeanAndStd(::Type{T}=Float64; fn::F=identity) where {T,F}
    AccMeanAndStd{T,F}(0, zero(T), zero(T), fn)
end

function (acc::AccMeanAndStd{T,F})() where {T,F}
    unbiased_std = sqrt(acc.svar / (acc.nobs - 1))
    (mean=acc.mean, std=unbiased_std)
end

function (acc::AccMeanAndStd{T,F})(x::T) where {T,F}
    xx = acc.fn(x)
    acc.nobs += 1
    prior_mean = acc.mean
    acc.mean = prior_mean + (xx - prior_mean) / acc.nobs
    acc.svar = acc.svar + (xx - prior_mean) * (xx - acc.mean)
    accum
end

function (acc::AccMeanAndStd{T,F})(xs::Seq{T}) where {T,F}
    xxs = map(acc.fn, xs)
    acc.nobs += length(xs)
    prior_mean = acc.mean
    xmean = vmean(xxs)
    acc.mean += (xmean - prior_mean) / acc.nobs
    acc.svar = acc.svar + (xmean - prior_mean) * (acc.mean - xmean)
    accum
end

# see https://www.johndcook.com/blog/skewness_kurtosis/

mutable struct AccStats{T,F} <: Accumulator{T,F}
    nobs::Int
    m1::T
    m2::T
    m3::T
    m4::T
    fn::F
end

AccumStats(::Type{T}=Float64; fn::F=identity) where {T,F} =
    AccStats(0, zero(T), zero(T), zero(T), zero(T), fn)

function (acc::AccStats{T,F})() where {T,F}
    (nobs=nobs(accum), mean=mean(accum), var=var(accum), std=std(accum), skewness=skewness(accum), kurtosis=kurtosis(accum))
end

function (acc::AccStats{T,F})(x) where {T,F}
    xx = acc.fn(x)
    n1 = acc.nobs
    acc.nobs += 1
    delta = xx - acc.m1
    delta_n = delta / acc.nobs
    delta_n2 = delta_n^2
    term1 = delta * delta_n * n1
    acc.m1 += delta_n
    acc.m4 += term1 * delta_n2 * (acc.nobs^2 - 3*acc.nobs + 3) + 
              6 * delta_n2 * acc.m2 - 
              4 * delta_n * acc.m3
    acc.m3 += term1 * delta_n * (acc.nobs - 2) - 3 * delta_n * acc.m2
    acc.m2 += term1
end

function (acc::AccStats{T,F})(xs::Seq{T}) where {T,F}
    for i in eachindex(xs)
        accum(xs[i])
    end
    accum
end

#=
reference for AccumExpWtMean, AccExpWtMeanVar

Incremental calculation of weighted mean and variance
by Tony Finch
=#

mutable struct AccExpWtMean{T,F} <: Accumulator{T,F}
    nobs::Int
    alpha::T
    expwtmean::T
    fn::F
end

AccumExpWtMean(::Type{T}=Float64; alpha::T=T(0.5), fn::F=identity) where {T,F} =
    AccExpWtMean{T,F}(0, T(alpha), zero(T), fn)

(accum::AccExpWtMean{T,F})() where {T,F} = acc.expwtmean

function (acc::AccExpWtMean{T,F})(x) where {T,F}
    xx = acc.fn(x)
    acc.nobs += 1
    acc.expwtmean += acc.alpha * (xx - acc.expwtmean)
    accum
end

function (acc::AccExpWtMean{T,F})(xs::Seq{T}) where {T,F}
    for x in eachindex(xs)
        accum(xs[i])
    end
    accum
end

mutable struct AccExpWtMeanVar{T,F} <: Accumulator{T,F}
    nobs::Int
    alpha::T
    expwtmean::T
    expwtsvar::T
    fn::F
end

AccumExpWtMeanVar(::Type{T}=Float64; alpha::T=T(0.5), fn::F=identity) where {T,F} =
    AccExpWtMeanVar(0, alpha, zero(T), zero(T), fn)

function(accum::AccExpWtMeanVar{T,F})() where {T,F}
    unbiased_expwtvar = acc.expwtsvar / (acc.nobs - 1)
    (expwt_mean=acc.expwtmean, expwt_var=unbiased_expwtvar)
end

function (acc::AccExpWtMeanVar{T,F})(x) where {T,F}
    xx = acc.fn(x)
    acc.nobs += 1
    diff = xx - acc.expwtmean
    incr = acc.alpha * diff
    acc.expwtmean += acc.alpha * (xx - acc.expwtmean)
    acc.expwtsvar = (one(T) - acc.alpha) * (acc.expwtsvar + diff * incr)
    accum
end

function (acc::AccExpWtMeanVar{T,F})(xs::Seq{T}) where {T,F}
    for i in eachindex(xs)
        accum(xs[i])
    end
    accum
end

mutable struct AccExpWtMeanStd{T,F} <: Accumulator{T,F}
    nobs::Int
    alpha::T
    expwtmean::T
    expwtsvar::T
    fn::F
end

AccumExpWtMeanStd(::Type{T}=Float64; alpha::T=T(0.5), fn::F=identity) where {T,F} =
    AccExpWtMeanStd(0, alpha, zero(T), zero(T), fn)

function(accum::AccExpWtMeanStd{T,F})() where {T,F}
    unbiased_expwtstd = sqrt(acc.expwtsvar / (acc.nobs - 1))
    (expwt_mean=acc.expwtmean, expwt_std=unbiased_expwtstd)
end

function (acc::AccExpWtMeanStd{T,F})(x) where {T,F}
    xx = acc.fn(x)
    acc.nobs += 1
    diff = xx - acc.expwtmean
    incr = acc.alpha * diff
    acc.expwtmean += acc.alpha * (xx - acc.expwtmean)
    acc.expwtsvar = (one(T) - acc.alpha) * (acc.expwtsvar + diff * incr)
    accum
end

function (acc::AccExpWtMeanStd{T,F})(xs::Seq{T}) where {T,F}
    for i in eachindex(xs)
        accum(xs[i])
    end
    accum
end

Base.length(@nospecialize acc::Accumulator) = acc.nobs
StatsBase.nobs(@nospecialize acc::Accumulator) = acc.nobs

for (F,A) in ((:(Base.minimum), :AccMinimum), (:(Base.maximum), :AccMaximum), (:(Base.extrema), :AccExtrema),
              (:(Base.sum), :AccSum), (:(Base.prod), :AccProd),
              (:(StatsBase.mean), :AccMean), (:(StatsBase.geomean), :AccGeoMean), (:(StatsBase.harmmean), :AccHarmMean))
     @eval $F(accum::$A) = Accum()
end

Base.minimum(accum::AccExtrema) = acc.min
Base.maximum(accum::AccExtrema) = acc.max
midrange(accum::AccExtrema) = (acc.min / 2) + (acc.max / 2)
proportionalrange(accum::AccExtrema, proportion) = (acc.min * proportion) + (acc.max * (1 - proportion))

nminima(accum::AccMinimum) = acc.nmin
nminima(accum::AccExtrema) = acc.nmin
nmaxima(accum::AccMinimum) = acc.nmax
nmaxima(accum::AccExtrema) = acc.nmax

StatsBase.mean(accum::AccMeanAndVar) = acc.mean
StatsBase.var(accum::AccMeanAndVar) = acc.svar / (acc.nobs - 1)
StatsBase.std(accum::AccMeanAndVar) = sqrt(acc.svar / (acc.nobs - 1))
StatsBase.mean(accum::AccMeanAndStd) = acc.mean
StatsBase.var(accum::AccMeanAndStd) = acc.svar / (acc.nobs - 1)
StatsBase.std(accum::AccMeanAndStd) = sqrt(acc.svar / (acc.nobs - 1))

StatsBase.mean(accum::AccStats{T}) where {T,F} = T(acc.m1)
StatsBase.var(accum::AccStats{T}) where {T,F} = T(acc.m2 / (acc.nobs - 1))
StatsBase.std(accum::AccStats{T}) where {T,F} = T(sqrt(var(accum)))
StatsBase.skewness(accum::AccStats{T}) where {T,F} = T(sqrt(acc.nobs) * acc.m3 / (acc.m2 * sqrt(acc.m2)))
StatsBase.kurtosis(accum::AccStats{T}) where {T,F} = T( ((acc.nobs * acc.m4) / (acc.m2^2)) - 3)

