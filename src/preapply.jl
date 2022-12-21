#=
     AccCount, 
     AccMinimum, AccMaximum, AccExtrema, 
     AccSum, AccProd,
     AccMean, AccGeoMean, AccHarmMean, AccGenMean,
     AccMeanAndVar, AccMeanAndStd, AccStats,
     AccExpWtMean, AccExpWtMeanVar, AccExpWtMeanStd
=#
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

function (Accum::AccumCount{I})() where {I}
    Accum.nobs
end
     
function (Accum::AccumCount{I})(x::T) where {I,T}
    Accum.nobs += one(I)
    Accum
end

function (Accum::AccumCount{I})(xs::Seq{T}) where {I,T}
    Accum.nobs += length(xs)
    Accum
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

function (Accum::AccumMinimum{T,F})() where {T,F}
    Accum.min
end

function (Accum::AccumMinimum{T})(x::T) where {T,F}
    Accum.nobs += 1
    if x < Accum.min
        Accum.nmin += 1
        Accum.min = x
    end
    Accum
end

function (Accum::AccumMinimum{T})(xs::Seq{T}) where {T,F}
    Accum.nobs += length(xs)     
    x = vminimum(xs)
    if x < Accum.min
        Accum.nmin += 1
        Accum.min = x
    end
    Accum
end

# Maximum

mutable struct AccumMaximum{T,F} <: AccumLensed{T,F}
    nobs::Int       # count each observation
    nmax::Int       # count distinct maxima
    max::T          # current maximum
end

function AccumMaximum(::Type{T}=AccumNum; fn::F=identity) where {T,F}
     AccumMaximum{T,F}(0, 0, typemin(T), fn)
end

function (Accum::AccumMaximum{T})() where {T,F}
    Accum.max
end

function (Accum::AccumMaximum{T})(x::T) where {T,F}
    Accum.nobs += 1
    if x > Accum.max
        Accum.nmax += 1
        Accum.max = x
    end
    Accum
end

function (Accum::AccumMaximum{T})(xs::Seq{T}) where {T,F}
    Accum.nobs += length(xs)     
    x = vmaximum(xs)
    if x > Accum.max
        Accum.nmax += 1
        Accum.max = x
    end
    Accum
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

function (Accum::AccumExtrema{T})() where {T,F}
    (Accum.min, Accum.max)
end

function (Accum::AccumExtrema{T})(x::T) where {T,F}
    Accum.nobs += 1
    if x < Accum.min
        Accum.nmin += 1
        Accum.min = x
    end
    if x > Accum.max
        Accum.nmax += 1
        Accum.max = x
    end
    Accum
end

function (Accum::AccumExtrema{T})(xs::Seq{T}) where {T,F}
    Accum.nobs += length(xs)     
    mn, mx = vextrema(xs)
    if mn < Accum.min
        Accum.nmin += 1
        Accum.min = mn
    end
    if mx > Accum.max
        Accum.nmax += 1
        Accum.max = mx
    end
    Accum
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

function (Accum::AccumSum{T})() where {T,F}
    Accum.sum
end

function (Accum::AccumSum{T})(x::T) where {T,F}
    Accum.nobs += 1
    Accum.sum += x
    Accum
end

function (Accum::AccumSum{T})(xs::Seq{T}) where {T,F}
    Accum.nobs += length(xs)     
    x = vsum(xs)
    Accum.sum += x
    Accum
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

function (Accum::AccumProd{T})() where {T,F}
    Accum.prod
end

function (Accum::AccumProd{T})(x::T) where {T,F}
    Accum.nobs += 1
    Accum.prod *= x
    Accum
end

function (Accum::AccumProd{T})(xs::Seq{T}) where {T,F}
    Accum.nobs += length(xs)     
    x = vprod(xs)
    Accum.prod *= x
    Accum
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

function (Accum::AccumMean{T})() where {T,F}
    Accum.mean
end

function (Accum::AccumMean{T})(x::T) where {T,F}
    Accum.nobs += 1
    Accum.mean += (x - Accum.mean) / Accum.nobs
    Accum
end

function (Accum::AccumMean{T})(xs::Seq{T}) where {T,F}
    Accum.nobs += length(xs)     
    xmean = vmean(xs)
    Accum.mean += (xmean - Accum.mean) / Accum.nobs
    Accum
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

function (Accum::AccumGeoMean{T})() where {T,F}
    n = ifelse(iszero(Accum.nobs), 1, Accum.nobs)
    exp(Accum.sumlog / n)
end

function (Accum::AccumGeoMean{T})(x::T) where {T,F}
    Accum.nobs += 1
    Accum.sumlog += logabs(x::T)
    Accum
end

function (Accum::AccumGeoMean{T})(xs::Seq{T}) where {T,F}
    Accum.nobs += length(xs)
    Accum.sumlog += sum(map(logabs, xs))
    Accum
end

# HarmMean
                                   
mutable struct AccumHarmMean{T,F} <: AccumLensed{T,F}
    nobs::Int       # count each observation
    invhmean::T     # 1 / current harmonic mean
    fn::F
end

function AccumHarmMean(::Type{T}=AccumNum; fn::F=identity) where {T,F}
     AccumHarmMean{T,F}(0, zero(T), fn)
end

function (Accum::AccumHarmMean{T})() where {T,F}
    n = ifelse(iszero(Accum.nobs), 1, Accum.nobs)
    n / Accum.invhmean
end

function (Accum::AccumHarmMean{T})(x::T) where {T,F}
    Accum.nobs += 1
    Accum.invhmean += one(T) / x
    Accum
end

function (Accum::AccumHarmMean{T})(xs::Seq{T}) where {T,F}
    Accum.nobs += length(xs)
    Accum.invhmean += sum(map(inv, xs))
    Accum
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
                                        
function (Accum::AccumGenMean{T})() where {T,F}
    (Accum.gmean)^Accum.rpwr
end

function (Accum::AccumGenMean{T})(x::T) where {T,F}
    Accum.nobs += 1
    Accum.gmean += ((x^Accum.pwr) - Accum.gmean) / Accum.nobs
    Accum
end

function (Accum::AccumGenMean{T})(xs::Seq{T}) where {T,F}
    Accum.nobs += length(xs)     
    xmean = vmean(map(x->x^Accum.pwr, xs))
    Accum.gmean += (xmean - Accum.gmean) / Accum.nobs
    Accum
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

function (Accum::AccumMeanAndVar{T})() where {T,F}
    unbiased_var = Accum.svar / (Accum.nobs - 1)
    (mean=Accum.mean, var=unbiased_var)
end

function (Accum::AccumMeanAndVar{T})(x::T) where {T,F}
    Accum.nobs += 1
    prior_mean = Accum.mean
    Accum.mean = prior_mean + (x - prior_mean) / Accum.nobs
    Accum.svar = Accum.svar + (x - prior_mean) * (x - Accum.mean)
    Accum
end

function (Accum::AccumMeanAndVar{T})(xs::Seq{T}) where {T,F}
    Accum.nobs += length(xs)
    prior_mean = Accum.mean
    xmean = vmean(xs)
    Accum.mean += (xmean - prior_mean) / Accum.nobs
    Accum.svar = Accum.svar + (x - prior_mean) * (x - xmean)
    Accum
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

function (Accum::AccumMeanAndStd{T})() where {T,F}
    unbiased_std = sqrt(Accum.svar / (Accum.nobs - 1))
    (mean=Accum.mean, std=unbiased_std)
end

function (Accum::AccumMeanAndStd{T})(x::T) where {T,F}
    Accum.nobs += 1
    prior_mean = Accum.mean
    Accum.mean = prior_mean + (x - prior_mean) / Accum.nobs
    Accum.svar = Accum.svar + (x - prior_mean) * (x - Accum.mean)
    Accum
end

function (Accum::AccumMeanAndStd{T})(xs::Seq{T}) where {T,F}
    Accum.nobs += length(xs)
    prior_mean = Accum.mean
    xmean = vmean(xs)
    Accum.mean += (xmean - prior_mean) / Accum.nobs
    Accum.svar = Accum.svar + (x - prior_mean) * (x - xmean)
    Accum
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

function (Accum::AccumStats{T})() where {T,F}
    (nobs=nobs(Accum), mean=mean(Accum), var=var(Accum), std=std(Accum), skewness=skewness(Accum), kurtosis=kurtosis(Accum))
end

#=
                                                  
function (Accum::AccumMeanVar{T})(x::T) where {T,F}
    Accum.nobs += 1
    prior_mean = Accum.mean
    Accum.mean = prior_mean + (x - prior_mean) / Accum.nobs
    Accum.svar = Accum.svar + (x - prior_mean) * (x - Accum.mean)
    Accum
end
=#
                                                  
function (Accum::AccumStats{T})(x) where {T,F}
    n1 = Accum.nobs
    Accum.nobs += 1
    delta = x - Accum.m1
    delta_n = delta / Accum.nobs
    delta_n2 = delta_n^2
    term1 = delta * delta_n * n1
    Accum.m1 += delta_n
    Accum.m4 += term1 * delta_n2 * (Accum.nobs^2 - 3*Accum.nobs + 3) + 
              6 * delta_n2 * Accum.m2 - 
              4 * delta_n * Accum.m3
    Accum.m3 += term1 * delta_n * (Accum.nobs - 2) - 3 * delta_n * Accum.m2
    Accum.m2 += term1
end

function (Accum::AccumStats{T})(xs::Seq{T}) where {T,F}
    for x in xs
        Accum(x)
    end
    Accum
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

(Accum::AccumExpWtMean{T})() where {T,F} = Accum.expwtmean

function (Accum::AccumExpWtMean{T})(x) where {T,F}
    Accum.nobs += 1
    Accum.expwtmean += Accum.alpha * (x - Accum.expwtmean)
    Accum
end

function (Accum::AccumExpWtMean{T})(xs::Seq{T}) where {T,F}
    for x in xs
        Accum(x)
    end
    Accum
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

function(Accum::AccumExpWtMeanVar{T})() where {T,F}
    unbiased_expwtvar = Accum.expwtsvar / (Accum.nobs - 1)
    (expwt_mean=Accum.expwtmean, expwt_var=unbiased_expwtvar)
end

function (Accum::AccumExpWtMeanVar{T})(x) where {T,F}
    Accum.nobs += 1
    diff = x - Accum.expwtmean
    incr = Accum.alpha * diff
    Accum.expwtmean += Accum.alpha * (x - Accum.expwtmean)
    Accum.expwtsvar = (one(T) - Accum.alpha) * (Accum.expwtsvar + diff * incr)
    Accum
end

function (Accum::AccumExpWtMeanVar{T})(xs::Seq{T}) where {T,F}
    for x in xs
        Accum(x)
    end
    Accum
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

function(Accum::AccumExpWtMeanStd{T})() where {T,F}
    unbiased_expwtstd = sqrt(Accum.expwtsvar / (Accum.nobs - 1))
    (expwt_mean=Accum.expwtmean, expwt_std=unbiased_expwtstd)
end

function (Accum::AccumExpWtMeanStd{T})(x) where {T,F}
    Accum.nobs += 1
    diff = x - Accum.expwtmean
    incr = Accum.alpha * diff
    Accum.expwtmean += Accum.alpha * (x - Accum.expwtmean)
    Accum.expwtsvar = (one(T) - Accum.alpha) * (Accum.expwtsvar + diff * incr)
    Accum
end

function (Accum::AccumExpWtMeanStd{T})(xs::Seq{T}) where {T,F}
    for x in xs
        Accum(x)
    end
    Accum
end

Base.length(@nospecialize Accum::Accumulator) = Accum.nobs
StatsBase.nobs(@nospecialize Accum::Accumulator) = Accum.nobs

for (F,A) in ((:(Base.minimum), :AccumMinimum), (:(Base.maximum), :AccumMaximum), (:(Base.extrema), :AccumExtrema),
              (:(Base.sum), :AccumSum), (:(Base.prod), :AccumProd),
              (:(StatsBase.mean), :AccumMean), (:(StatsBase.geomean), :AccumGeoMean), (:(StatsBase.harmmean), :AccumHarmMean))
     @eval $F(Accum::$A) = Accum()
end

Base.minimum(Accum::AccumExtrema) = Accum.min
Base.maximum(Accum::AccumExtrema) = Accum.max
midrange(Accum::AccumExtrema) = (Accum.min / 2) + (Accum.max / 2)
proportionalrange(Accum::AccumExtrema, proportion) = (Accum.min * proportion) + (Accum.max * (1 - proportion))

nminima(Accum::AccumMinimum) = Accum.nmin
nminima(Accum::AccumExtrema) = Accum.nmin
nmaxima(Accum::AccumMinimum) = Accum.nmax
nmaxima(Accum::AccumExtrema) = Accum.nmax

StatsBase.mean(Accum::AccumMeanAndVar) = Accum.mean
StatsBase.var(Accum::AccumMeanAndVar) = Accum.svar / (Accum.nobs - 1)
StatsBase.std(Accum::AccumMeanAndVar) = sqrt(Accum.svar / (Accum.nobs - 1))
StatsBase.mean(Accum::AccumMeanAndStd) = Accum.mean
StatsBase.var(Accum::AccumMeanAndStd) = Accum.svar / (Accum.nobs - 1)
StatsBase.std(Accum::AccumMeanAndStd) = sqrt(Accum.svar / (Accum.nobs - 1))

StatsBase.mean(Accum::AccumStats{T}) where {T,F} = T(Accum.m1)
StatsBase.var(Accum::AccumStats{T}) where {T,F} = T(Accum.m2 / (Accum.nobs - 1))
StatsBase.std(Accum::AccumStats{T}) where {T,F} = T(sqrt(var(Accum)))
StatsBase.skewness(Accum::AccumStats{T}) where {T,F} = T(sqrt(Accum.nobs) * Accum.m3 / (Accum.m2 * sqrt(Accum.m2)))
StatsBase.kurtosis(Accum::AccumStats{T}) where {T,F} = T( ((Accum.nobs * Accum.m4) / (Accum.m2^2)) - 3)


