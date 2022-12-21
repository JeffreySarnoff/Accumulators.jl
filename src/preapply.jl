#=
     AccumCount, 
     AccumMinimum, AccumMaximum, AccumExtrema, 
     AccumSum, AccumProd,
     AccumMean, AccumGeoMean, AccumHarmMean, AccumGenMean,
     AccumMeanAndVar, AccumMeanAndStd, AccumStats,
     AccumExpWtMean, AccumExpWtMeanVar, AccumExpWtMeanStd
=#

# Count

mutable struct AccumCount{I,F} <: Accumulator{I}
    nobs::Int       # count each observation
    fn::F           # preapply to each element when observed
end

function AccumCount(::Type{I}=Int64) where {I,F}
     AccumCount{I,F}(zero(I))
end

function (accum::AccumCount{I})() where {I}
    accum.nobs
end
     
function (accum::AccumCount{I})(x::T) where {I,T}
    accum.nobs += one(I)
    accum
end

function (accum::AccumCount{I})(xs::Seq{T}) where {I,T}
    accum.nobs += length(xs)
    accum
end

# Minimum

mutable struct AccumMinimum{T,F} <: AccumLensed{T,F}
    nobs::Int       # count each observation
    nmin::Int       # count distinct minima
    min::T          # current minimum
    fn::F           # preapply to each element when observed
end

function AccumMinimum(::Type{T}=AccumNum; fn::F=identity) where {T,F}
     AccumMinimum{T,F}(0, 0, typemax(T), fn)
end

function (accum::AccumMinimum{T,F})() where {T,F}
    accum.min
end

function (accum::AccumMinimum{T,F})(x::T) where {T,F}
    xx = accum.fn(x)
    accum.nobs += 1
    if x < accum.min
        accum.nmin += 1
        accum.min = xx
    end
    accum
end

function (accum::AccumMinimum{T,F})(xs::Seq{T}) where {T,F}
    xxs = map(accum.fn, xs)
    accum.nobs += length(xs)
    x = vminimum(xxs)
    if x < accum.min
        accum.nmin += 1
        accum.min = x
    end
    accum
end

# Maximum

mutable struct AccumMaximum{T,F} <: AccumLensed{T,F}
    nobs::Int       # count each observation
    nmax::Int       # count distinct maxima
    max::T          # current maximum
    fn::F
end

function AccumMaximum(::Type{T}=AccumNum; fn::F=identity) where {T,F}
     AccumMaximum{T,F}(0, 0, typemin(T), fn)
end

function (accum::AccumMaximum{T,F})() where {T,F}
    accum.max
end

function (accum::AccumMaximum{T,F})(x::T) where {T,F}
    xx = accum.fn(x)
    accum.nobs += 1
    if xx > accum.max
        accum.nmax += 1
        accum.max = xx
    end
    accum
end

function (accum::AccumMaximum{T,F})(xs::Seq{T}) where {T,F}
    accum.nobs += length(xs)     
    x = vmaximum(xs)
    if x > accum.max
        accum.nmax += 1
        accum.max = x
    end
    accum
end

# Extrema

mutable struct AccumExtrema{T,F} <: AccumLensed{T,F}
    nobs::Int       # count each observation
    nmin::Int       # count distinct minima
    nmax::Int       # count distinct maxima
    min::T          # current minimum
    max::T          # current maximum
    fn::F
end

function AccumExtrema(::Type{T}=AccumNum; fn::F=identity) where {T,F}
     AccumExtrema{T,F}(0, 0, 0, typemax(T), typemin(T), fn)
end

function (accum::AccumExtrema{T,F})() where {T,F}
    (accum.min, accum.max)
end

function (accum::AccumExtrema{T,F})(x::T) where {T,F}
    xx = accum.fn(x)
    accum.nobs += 1
    if xx < accum.min
        accum.nmin += 1
        accum.min = xx
    end
    if xx > accum.max
        accum.nmax += 1
        accum.max = xx
    end
    accum
end

function (accum::AccumExtrema{T,F})(xs::Seq{T}) where {T,F}
    accum.nobs += length(xs)     
    mn, mx = vextrema(xs)
    if mn < accum.min
        accum.nmin += 1
        accum.min = mn
    end
    if mx > accum.max
        accum.nmax += 1
        accum.max = mx
    end
    accum
end

# Sum

mutable struct AccumSum{T,F} <: AccumLensed{T,F}
    nobs::Int       # count each observation
    sum::T          # current sum
    fn::F
end

function AccumSum(::Type{T}=AccumNum; fn::F=identity) where {T,F}
     AccumSum{T,F}(0, zero(T), fn)
end

function (accum::AccumSum{T,F})() where {T,F}
    accum.sum
end

function (accum::AccumSum{T,F})(x::T) where {T,F}
    xx = accum.fn(x)
    accum.nobs += 1
    accum.sum += xx
    accum
end

function (accum::AccumSum{T,F})(xs::Seq{T}) where {T,F}
    accum.nobs += length(xs)     
    x = vsum(xs)
    accum.sum += x
    accum
end

# Prod

mutable struct AccumProd{T,F} <: AccumLensed{T,F}
    nobs::Int       # count each observation
    prod::T         # current product
    fn::F
end

function AccumProd(::Type{T}=AccumNum; fn::F=identity) where {T,F}
     AccumProd{T,F}(0, one(T), fn)
end

function (accum::AccumProd{T,F})() where {T,F}
    accum.prod
end

function (accum::AccumProd{T,F})(x::T) where {T,F}
    xx = accum.fn(x)
    accum.nobs += 1
    accum.prod *= xx
    accum
end

function (accum::AccumProd{T,F})(xs::Seq{T}) where {T,F}
    accum.nobs += length(xs)     
    x = vprod(xs)
    accum.prod *= x
    accum
end

# Mean

mutable struct AccumMean{T,F} <: AccumLensed{T,F}
    nobs::Int       # count each observation
    mean::T         # current mean
    fn::F
end

function AccumMean(::Type{T}=AccumNum; fn::F=identity) where {T,F}
     AccumMean{T,F}(0, zero(T), fn)
end

function (accum::AccumMean{T,F})() where {T,F}
    accum.mean
end

function (accum::AccumMean{T,F})(x::T) where {T,F}
    xx = accum.fn(x)
    accum.nobs += 1
    accum.mean += (xx - accum.mean) / accum.nobs
    accum
end

function (accum::AccumMean{T,F})(xs::Seq{T}) where {T,F}
    accum.nobs += length(xs)     
    xmean = vmean(xs)
    accum.mean += (xmean - accum.mean) / accum.nobs
    accum
end

# GeoMean

mutable struct AccumGeoMean{T,F} <: AccumLensed{T,F}
    nobs::Int       # count each observation
    sumlog::T       # ∑(i=1:nobs) log(xᵢ)
    fn::F
end

function AccumGeoMean(::Type{T}=AccumNum; fn::F=identity) where {T,F}
     AccumGeoMean{T,F}(0, zero(T), fn)
end

function (accum::AccumGeoMean{T,F})() where {T,F}
    n = ifelse(iszero(accum.nobs), 1, accum.nobs)
    exp(accum.sumlog / n)
end

function (accum::AccumGeoMean{T,F})(x::T) where {T,F}
    xx = accum.fn(x)
    accum.nobs += 1
    accum.sumlog += logabs(xx)
    accum
end

function (accum::AccumGeoMean{T,F})(xs::Seq{T}) where {T,F}
    accum.nobs += length(xs)
    accum.sumlog += sum(map(logabs, xs))
    accum
end

# Harmonic Mean
mutable struct AccumHarmMean{T,F} <: AccumLensed{T,F}
    nobs::Int       # count each observation
    invhmean::T     # 1 / current harmonic mean
    fn::F
end

function AccumHarmMean(::Type{T}=AccumNum; fn::F=identity) where {T,F}
     AccumHarmMean{T,F}(0, zero(T), fn)
end

function (accum::AccumHarmMean{T,F})() where {T,F}
    n = ifelse(iszero(accum.nobs), 1, accum.nobs)
    n / accum.invhmean
end

function (accum::AccumHarmMean{T,F})(x::T) where {T,F}
    xx = accum.fn(x)
    accum.nobs += 1
    accum.invhmean += one(T) / xx
    accum
end

function (accum::AccumHarmMean{T,F})(xs::Seq{T}) where {T,F}
    accum.nobs += length(xs)
    accum.invhmean += sum(map(inv, xs))
    accum
end

# Generalized Mean (defaults to Quadratic Mean [root-mean-squared])

mutable struct AccumGenMean{T,F} <: AccumLensed{T,F}
    nobs::Int       # count each observation
    gmean::T        # current mean((xᵢ)ᵖʷʳ), i=1:nobs
    const pwr::T    # power
    const rpwr::T   # reciprocal of power
    fn::F
end

function AccumGenMean(::Type{T}=AccumNum; ; fn::F=identity, power::Real=2.0) where {T,F}
    AccumGenMean{T,F}(0, zero(T), T(power), T(1/power), fn)
end
                                        
function (accum::AccumGenMean{T,F})() where {T,F}
    (accum.gmean)^accum.rpwr
end

function (accum::AccumGenMean{T,F})(x::T) where {T,F}
    xx = accum.fn(x)
    accum.nobs += 1
    accum.gmean += ((xx^accum.pwr) - accum.gmean) / accum.nobs
    accum
end

function (accum::AccumGenMean{T,F})(xs::Seq{T}) where {T,F}
    accum.nobs += length(xs)     
    xmean = vmean(map(x->x^accum.pwr, xs))
    accum.gmean += (xmean - accum.gmean) / accum.nobs
    accum
end
              
# Unbiased Sample Variation (with Mean)
# see https://www.johndcook.com/blog/standard_deviation/

mutable struct AccumMeanAndVar{T,F} <: AccumLensed{T,F}
    nobs::Int       # count each observation
    mean::T         # current mean
    svar::T         # sum of variances x[1:1=0,1:2,..,1:nobs]
    fn::F
end

function AccumMeanAndVar(::Type{T}=Float64; fn::F=identity) where {T,F}
    AccumMeanAndVar{T,F}(0, zero(T), zero(T), fn)
end

function (accum::AccumMeanAndVar{T,F})() where {T,F}
    unbiased_var = accum.svar / (accum.nobs - 1)
    (mean=accum.mean, var=unbiased_var)
end

function (accum::AccumMeanAndVar{T,F})(x::T) where {T,F}
    xx = accum.fn(x)
    accum.nobs += 1
    prior_mean = accum.mean
    accum.mean = prior_mean + (xx - prior_mean) / accum.nobs
    accum.svar = accum.svar + (xx - prior_mean) * (xx - accum.mean)
    accum
end

function (accum::AccumMeanAndVar{T,F})(xs::Seq{T}) where {T,F}
    accum.nobs += length(xs)
    prior_mean = accum.mean
    xmean = vmean(xs)
    accum.mean += (xmean - prior_mean) / accum.nobs
    accum.svar = accum.svar + (x - prior_mean) * (x - xmean)
    accum
end

mutable struct AccumMeanAndStd{T,F} <: AccumLensed{T,F}
    nobs::Int       # count each observation
    mean::T         # current mean
    svar::T         # sum of variances x[1:1=0,1:2,..,1:nobs]
    fn::F
end

function AccumMeanAndStd(::Type{T}=Float64; fn::F=identity) where {T,F}
    AccumMeanAndStd{T,F}(0, zero(T), zero(T), fn)
end

function (accum::AccumMeanAndStd{T,F})() where {T,F}
    unbiased_std = sqrt(accum.svar / (accum.nobs - 1))
    (mean=accum.mean, std=unbiased_std)
end

function (accum::AccumMeanAndStd{T,F})(x::T) where {T,F}
    xx = accum.fn(x)
    accum.nobs += 1
    prior_mean = accum.mean
    accum.mean = prior_mean + (xx - prior_mean) / accum.nobs
    accum.svar = accum.svar + (xx - prior_mean) * (x - accum.mean)
    accum
end

function (accum::AccumMeanAndStd{T,F})(xs::Seq{T}) where {T,F}
    accum.nobs += length(xs)
    prior_mean = accum.mean
    xmean = vmean(xs)
    accum.mean += (xmean - prior_mean) / accum.nobs
    accum.svar = accum.svar + (x - prior_mean) * (x - xmean)
    accum
end

# see https://www.johndcook.com/blog/skewness_kurtosis/

mutable struct AccumStats{T,F} <: AccumLensed{T,F}
    nobs::Int
    m1::T
    m2::T
    m3::T
    m4::T
    fn::F
end

AccumStats(::Type{T}=Float64; fn::F=identity) where {T,F} =
    AccumStats(0, zero(T), zero(T), zero(T), zero(T), fn)

function (accum::AccumStats{T,F})() where {T,F}
    (nobs=nobs(accum), mean=mean(accum), var=var(accum), std=std(accum), skewness=skewness(accum), kurtosis=kurtosis(accum))
end

#=
                                                  
function (accum::AccumMeanVar{T,F})(x::T) where {T,F}
    accum.nobs += 1
    prior_mean = accum.mean
    accum.mean = prior_mean + (x - prior_mean) / accum.nobs
    accum.svar = accum.svar + (x - prior_mean) * (x - accum.mean)
    accum
end
=#
                                                  
function (accum::AccumStats{T,F})(x) where {T,F}
    xx = accum.fn(x)
    n1 = accum.nobs
    accum.nobs += 1
    delta = xx - accum.m1
    delta_n = delta / accum.nobs
    delta_n2 = delta_n^2
    term1 = delta * delta_n * n1
    accum.m1 += delta_n
    accum.m4 += term1 * delta_n2 * (accum.nobs^2 - 3*accum.nobs + 3) + 
              6 * delta_n2 * accum.m2 - 
              4 * delta_n * accum.m3
    accum.m3 += term1 * delta_n * (accum.nobs - 2) - 3 * delta_n * accum.m2
    accum.m2 += term1
end

function (accum::AccumStats{T,F})(xs::Seq{T}) where {T,F}
    for x in xs
        Accum(x)
    end
    accum
end

#=
reference for AccumExpWtMean, AccumExpWtMeanVar

Incremental calculation of weighted mean and variance
by Tony Finch
=#

mutable struct AccumExpWtMean{T,F} <: AccumLensed{T,F}
    nobs::Int
    alpha::T
    expwtmean::T
    fn::F
end

AccumExpWtMean(::Type{T}=Float64; alpha::T=T(0.5); fn::F=identity) where {T,F} =
    AccumExpWtMean{T,F}(0, T(alpha), zero(T), fn)

(accum::AccumExpWtMean{T,F})() where {T,F} = accum.expwtmean

function (accum::AccumExpWtMean{T,F})(x) where {T,F}
    xx = accum.fn(x)
    accum.nobs += 1
    accum.expwtmean += accum.alpha * (xx - accum.expwtmean)
    accum
end

function (accum::AccumExpWtMean{T,F})(xs::Seq{T}) where {T,F}
    for x in xs
        Accum(x)
    end
    accum
end

mutable struct AccumExpWtMeanVar{T,F} <: AccumLensed{T,F}
    nobs::Int
    alpha::T
    expwtmean::T
    expwtsvar::T
    fn::F
end

AccumExpWtMeanVar(::Type{T}=Float64; alpha::T=T(0.5); fn::F=identity) where {T,F} =
    AccumExpWtMeanVar(0, alpha, zero(T), zero(T), fn)

function(accum::AccumExpWtMeanVar{T,F})() where {T,F}
    unbiased_expwtvar = accum.expwtsvar / (accum.nobs - 1)
    (expwt_mean=accum.expwtmean, expwt_var=unbiased_expwtvar)
end

function (accum::AccumExpWtMeanVar{T,F})(x) where {T,F}
    xx = accum.fn(x)
    accum.nobs += 1
    diff = xx - accum.expwtmean
    incr = accum.alpha * diff
    accum.expwtmean += accum.alpha * (xx - accum.expwtmean)
    accum.expwtsvar = (one(T) - accum.alpha) * (accum.expwtsvar + diff * incr)
    accum
end

function (accum::AccumExpWtMeanVar{T,F})(xs::Seq{T}) where {T,F}
    for x in xs
        accum(x)
    end
    accum
end

mutable struct AccumExpWtMeanStd{T,F} <: AccumLensed{T,F}
    nobs::Int
    alpha::T
    expwtmean::T
    expwtsvar::T
    fn::F
end

AccumExpWtMeanStd(::Type{T}=Float64; alpha::T=T(0.5); fn::F=identity) where {T,F} =
    AccumExpWtMeanStd(0, alpha, zero(T), zero(T), fn)

function(accum::AccumExpWtMeanStd{T,F})() where {T,F}
    unbiased_expwtstd = sqrt(accum.expwtsvar / (accum.nobs - 1))
    (expwt_mean=accum.expwtmean, expwt_std=unbiased_expwtstd)
end

function (accum::AccumExpWtMeanStd{T,F})(x) where {T,F}
    xx = accum.fn(x)
    accum.nobs += 1
    diff = xx - accum.expwtmean
    incr = accum.alpha * diff
    accum.expwtmean += accum.alpha * (xx - accum.expwtmean)
    accum.expwtsvar = (one(T) - accum.alpha) * (accum.expwtsvar + diff * incr)
    accum
end

function (accum::AccumExpWtMeanStd{T,F})(xs::Seq{T}) where {T,F}
    for x in xs
        Accum(x)
    end
    accum
end

Base.length(@nospecialize Accum::Accumulator) = accum.nobs
StatsBase.nobs(@nospecialize Accum::Accumulator) = accum.nobs

for (F,A) in ((:(Base.minimum), :AccumMinimum), (:(Base.maximum), :AccumMaximum), (:(Base.extrema), :AccumExtrema),
              (:(Base.sum), :AccumSum), (:(Base.prod), :AccumProd),
              (:(StatsBase.mean), :AccumMean), (:(StatsBase.geomean), :AccumGeoMean), (:(StatsBase.harmmean), :AccumHarmMean))
     @eval $F(accum::$A) = Accum()
end

Base.minimum(accum::AccumExtrema) = accum.min
Base.maximum(accum::AccumExtrema) = accum.max
midrange(accum::AccumExtrema) = (accum.min / 2) + (accum.max / 2)
proportionalrange(accum::AccumExtrema, proportion) = (accum.min * proportion) + (accum.max * (1 - proportion))

nminima(accum::AccumMinimum) = accum.nmin
nminima(accum::AccumExtrema) = accum.nmin
nmaxima(accum::AccumMinimum) = accum.nmax
nmaxima(accum::AccumExtrema) = accum.nmax

StatsBase.mean(accum::AccumMeanAndVar) = accum.mean
StatsBase.var(accum::AccumMeanAndVar) = accum.svar / (accum.nobs - 1)
StatsBase.std(accum::AccumMeanAndVar) = sqrt(accum.svar / (accum.nobs - 1))
StatsBase.mean(accum::AccumMeanAndStd) = accum.mean
StatsBase.var(accum::AccumMeanAndStd) = accum.svar / (accum.nobs - 1)
StatsBase.std(accum::AccumMeanAndStd) = sqrt(accum.svar / (accum.nobs - 1))

StatsBase.mean(accum::AccumStats{T}) where {T,F} = T(accum.m1)
StatsBase.var(accum::AccumStats{T}) where {T,F} = T(accum.m2 / (accum.nobs - 1))
StatsBase.std(accum::AccumStats{T}) where {T,F} = T(sqrt(var(accum)))
StatsBase.skewness(accum::AccumStats{T}) where {T,F} = T(sqrt(accum.nobs) * accum.m3 / (accum.m2 * sqrt(accum.m2)))
StatsBase.kurtosis(accum::AccumStats{T}) where {T,F} = T( ((accum.nobs * accum.m4) / (accum.m2^2)) - 3)




