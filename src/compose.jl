#=
     AccCount, 
     AccMinimum, AccMaximum, AccExtrema, 
     AccSum, AccProd,
     AccMean, AccGeoMean, AccHarmMean, AccGenMean,
     AccMeanVar, AccMeanStd, AccStats,
     AccExpWtMean, AccExpWtMeanVar, AccExpWtMeanStd
=#
#=
names used in StatsBase
:adjr2, :adjr², :aic, :aicc, :autocor, :autocov, :aweights, 
:bic, 
:coef, :coefnames, :coeftable, :competerank, :confint, :cooksdistance, :cor, :cor2cov, :corkendall, :corspearman,
:counteq, :countmap, :countne, :counts, :cov, :cov2cor, :cronbachalpha, 
:crosscor, :crosscov, :crossentropy, :crossmodelmatrix, :cumulant, 
:denserank, :describe, :deviance, :dof, :dof_residual, 
:ecdf, :entropy, :eweights, 
:fit, :fitted, :fweights, 
:genmean, :genvar, :geomean, :gkldiv, 
:harmmean, 
:indexmap, :indicatormat, :informationmatrix, :inverse_rle, :iqr, :isfitted, :islinear, 
:kldivergence, :kurtosis, 
:levelsmap, :leverage, :loglikelihood, 
:mad, :maxad, :mean, :mean_and_cov, :mean_and_std, :mean_and_var, :meanad, :meanresponse, 
:median, :middle, :midpoints, :mode, :model_response, :modelmatrix, :modes, :moment, :msd, :mss, 
:nobs, :norepeats, :nquantile, :nulldeviance, :nullloglikelihood, 
:ordinalrank, 
:pacf, :pairwise, :partialcor, :percentile, :percentilerank, :predict, 
:proportionmap, :proportions, :psnr, :pweights, 
:quantile, :quantilerank, 
:r2, :renyientropy, :residuals, :response, :responsename, :rle, :rmsd, :rss, :r², 
:sample, :samplepair, :scattermat, :scattermat_zm, :scattermatm, :score, :sem, 
:skewness, :span, :sqL2dist, :standardize, :std, :stderror, :sum, :summarystats, 
:tiedrank, :totalvar, :trim, :trimvar, 
:uweights, 
:values, :var, :variation, :vcov, 
:weights, :winsor, :wmedian, :wquantile, :wsample, :wsum, 
:zscore
=#

# Count

mutable struct AccCount{T} <: Accumulator{T}
    nobs::Int       # count each observation
end

function AccCount(::Type{T}=Int64) where {T}
     AccCount{T}(zero(T))
end

function (acc::AccCount{T})() where {T}
    acc.nobs
end
     
function (acc::AccCount{T})(x::T) where {T}
    acc.nobs += one(T)
    acc
end

function (acc::AccCount{T})(xs::Seq{T}) where {T}
    acc.nobs += length(xs)
    acc
end

# Min

mutable struct AccMinimum{T} <: Accumulator{T}
    nobs::Int       # count each observation
    nmin::Int       # count distinct minima
    min::T          # current minimum
end

function AccMinimum(::Type{T}=AccNum) where {T}
     AccMinimum{T}(0, 0, typemax(T))
end

function (acc::AccMinimum{T})() where {T}
    acc.min
end

function (acc::AccMinimum{T})(x::T) where {T}
    acc.nobs += 1
    if x < acc.min
        acc.nmin += 1
        acc.min = x
    end
    acc
end

function (acc::AccMinimum{T})(xs::Seq{T}) where {T}
    acc.nobs += length(xs)     
    x = vminimum(xs)
    if x < acc.min
        acc.nmin += 1
        acc.min = x
    end
    acc
end

# Max

mutable struct AccMaximum{T} <: Accumulator{T}
    nobs::Int       # count each observation
    nmax::Int       # count distinct maxima
    max::T          # current maximum
end

function AccMaximum(::Type{T}=AccNum) where {T}
     AccMaximum{T}(0, 0, typemin(T))
end

function (acc::AccMaximum{T})() where {T}
    acc.max
end

function (acc::AccMaximum{T})(x::T) where {T}
    acc.nobs += 1
    if x > acc.max
        acc.nmax += 1
        acc.max = x
    end
    acc
end

function (acc::AccMaximum{T})(xs::Seq{T}) where {T}
    acc.nobs += length(xs)     
    x = vmaximum(xs)
    if x > acc.max
        acc.nmax += 1
        acc.max = x
    end
    acc
end

# Extrema

mutable struct AccExtrema{T} <: Accumulator{T}
    nobs::Int       # count each observation
    nmin::Int       # count distinct minima
    nmax::Int       # count distinct maxima
    min::T          # current minimum
    max::T          # current maximum
end

function AccExtrema(::Type{T}=AccNum) where {T}
     AccExtrema{T}(0, 0, 0, typemax(T), typemin(T))
end

function (acc::AccExtrema{T})() where {T}
    (acc.min, acc.max)
end

function (acc::AccExtrema{T})(x::T) where {T}
    acc.nobs += 1
    if x < acc.min
        acc.nmin += 1
        acc.min = x
    end
    acc
end

function (acc::AccExtrema{T})(xs::Seq{T}) where {T}
    acc.nobs += length(xs)     
    mn, mx = vextrema(xs)
    if mn < acc.min
        acc.nmin += 1
        acc.min = mn
    end
    if mx > acc.max
        acc.nmax += 1
        acc.max = mx
    end
    acc
end

# Sum

mutable struct AccSum{T} <: Accumulator{T}
    nobs::Int       # count each observation
    sum::T          # current sum
end

function AccSum(::Type{T}=AccNum) where {T}
     AccSum{T}(0, zero(T))
end

function (acc::AccSum{T})() where {T}
    acc.sum
end

function (acc::AccSum{T})(x::T) where {T}
    acc.nobs += 1
    acc.sum += x
    acc
end

function (acc::AccSum{T})(xs::Seq{T}) where {T}
    acc.nobs += length(xs)     
    x = vsum(xs)
    acc.sum += x
    acc
end

# Prod

mutable struct AccProd{T} <: Accumulator{T}
    nobs::Int       # count each observation
    prod::T         # current product
end

function AccProd(::Type{T}=AccNum) where {T}
     AccProd{T}(0, one(T))
end

function (acc::AccProd{T})() where {T}
    acc.prod
end

function (acc::AccProd{T})(x::T) where {T}
    acc.nobs += 1
    acc.prod *= x
    acc
end

function (acc::AccProd{T})(xs::Seq{T}) where {T}
    acc.nobs += length(xs)     
    x = vprod(xs)
    acc.prod *= x
    acc
end

# Mean

mutable struct AccMean{T} <: Accumulator{T}
    nobs::Int       # count each observation
    mean::T         # current mean
end

function AccMean(::Type{T}=AccNum) where {T}
     AccMean{T}(0, zero(T))
end

function (acc::AccMean{T})() where {T}
    acc.mean
end

function (acc::AccMean{T})(x::T) where {T}
    acc.nobs += 1
    acc.mean += (x - acc.mean) / acc.nobs
    acc
end

function (acc::AccMean{T})(xs::Seq{T}) where {T}
    acc.nobs += length(xs)     
    xmean = vmean(xs)
    acc.mean += (xmean - acc.mean) / acc.nobs
    acc
end

# GeoMean

mutable struct AccGeoMean{T} <: Accumulator{T}
    nobs::Int       # count each observation
    sumlog::T       # ∑(i=1:nobs) log(xᵢ)
end

function AccGeoMean(::Type{T}=AccNum) where {T}
     AccGeoMean{T}(0, one(T))
end

function (acc::AccGeoMean{T})() where {T}
    n = ifelse(acc.nobs === 0 ? 1 : acc.nobs)
    exp(acc.sumlog / n)
end

function (acc::AccGeoMean{T})(x::T) where {T}
    acc.nobs += 1
    acc.sumlog += logabs(x::T)
    acc
end

function (acc::AccGeoMean{T})(xs::Seq{T}) where {T}
    acc.nobs += length(xs)
    acc.sumlog += sum(map(logabs, xs))
    acc
end

# HarmMean
                                   
mutable struct AccHarmMean{T} <: Accumulator{T}
    nobs::Int       # count each observation
    invhmean::T     # 1 / current harmonic mean
end

function AccHarmMean(::Type{T}=AccNum) where {T}
     AccHarmMean{T}(0, zero(T))
end

function (acc::AccHarmMean{T})() where {T}
    n = ifelse(acc.nobs === 0 ? 1 : acc.nobs)
    n / acc.invhmean
end

function (acc::AccHarmMean{T})(x::T) where {T}
    acc.nobs += 1
    acc.invhmean += one(T) / x
    acc
end

function (acc::AccHarmMean{T})(xs::Seq{T}) where {T}
    acc.nobs += length(xs)
    acc.invhmean += sum(map(inv, xs))
    acc
end

# Generalized Mean

mutable struct AccGenMean{T} <: Accumulator{T}
    nobs::Int       # count each observation
    gmean::T        # current mean((xᵢ)ᵖʷʳ), i=1:nobs
    const pwr::T    # power
    const rpwr::T   # reciprocal of power
end

function AccGenMean(::Type{T}=AccNum; pwr::Real) where {T}
    AccGenMean{T}(0, zero(T), T(pwr), T(1/pwr))
end
                                        
function (acc::AccGenMean{T})() where {T}
    (acc.gmean)^acc.rpwr
end

function (acc::AccGenMean{T})(x::T) where {T}
    acc.nobs += 1
    acc.gmean += ((x^acc.pwr) - acc.gmean) / acc.nobs
    acc
end

function (acc::AccGenMean{T})(xs::Seq{T}) where {T}
    acc.nobs += length(xs)     
    xmean = vmean(map(x->x^acc.pwr, xs))
    acc.gmean += (xmean - acc.gmean) / acc.nobs
    acc
end
              
# Unbiased Sample Variation (with Mean)
# see https://www.johndcook.com/blog/standard_deviation/

mutable struct AccMeanVar{T} <: Accumulator{T}
    nobs::Int       # count each observation
    mean::T         # current mean
    svar::T         # sum of variances x[1:1=0,1:2,..,1:nobs]
end

function AccMeanVar(::Type{T}=Float64) where {T}
    AccMeanVar{T}(0, zero(T), zero(T))
end

function (acc::AccMeanVar{T})() where {T}
    unbiased_var = acc.svar / (acc.nobs - 1)
    (mean=acc.mean, var=unbiased_var)
end

function (acc::AccMeanVar{T})(x::T) where {T}
    acc.nobs += 1
    prior_mean = acc.mean
    acc.mean = prior_mean + (x - prior_mean) / acc.nobs
    acc.svar = acc.svar + (x - prior_mean) * (x - acc.mean)
    acc
end

function (acc::AccMeanVar{T})(xs::Seq{T}) where {T}
    acc.nobs += length(xs)
    prior_mean = acc.mean
    xmean = vmean(xs)
    acc.mean += (xmean - prior_mean) / acc.nobs
    acc.svar = acc.svar + (x - prior_mean) * (x - xmean)
    acc
end

mutable struct AccMeanStd{T} <: Accumulator{T}
    nobs::Int       # count each observation
    mean::T         # current mean
    svar::T         # sum of variances x[1:1=0,1:2,..,1:nobs]
end

function AccMeanStd(::Type{T}=Float64) where {T}
    AccMeanVar{T}(0, zero(T), zero(T))
end

function (acc::AccMeanStd{T})() where {T}
    unbiased_std = sqrt(acc.svar / (acc.nobs - 1))
    (mean=acc.mean, std=unbiased_std)
end

function (acc::AccMeanStd{T})(x::T) where {T}
    acc.nobs += 1
    prior_mean = acc.mean
    acc.mean = prior_mean + (x - prior_mean) / acc.nobs
    acc.svar = acc.svar + (x - prior_mean) * (x - acc.mean)
    acc
end

function (acc::AccMeanStd{T})(xs::Seq{T}) where {T}
    acc.nobs += length(xs)
    prior_mean = acc.mean
    xmean = vmean(xs)
    acc.mean += (xmean - prior_mean) / acc.nobs
    acc.svar = acc.svar + (x - prior_mean) * (x - xmean)
    acc
end

# see https://www.johndcook.com/blog/skewness_kurtosis/

mutable struct AccStats{T} <: Accumulator{T}
    nobs::Int
    m1::T
    m2::T
    m3::T
    m4::T
end

AccStats(::Type{T}=Float64) where {T} = AccStats(0, zero(T), zero(T), zero(T), zero(T))

function (acc::AccStats{T})() where {T}
    (count=count(acc), mean=mean(acc), var=var(acc), std=std(acc), skewness=skewness(acc), kurtosis=kurtosis(acc))
end

#=
                                                  
function (acc::AccMeanVar{T})(x::T) where {T}
    acc.nobs += 1
    prior_mean = acc.mean
    acc.mean = prior_mean + (x - prior_mean) / acc.nobs
    acc.svar = acc.svar + (x - prior_mean) * (x - acc.mean)
    acc
end
=#
                                                  
function (acc::AccStats{T})(x) where {T}
    n1 = acc.nobs
    acc.nobs += 1
    delta = x - acc.m1
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

function (acc::AccStats{T})(xs::Seq{T}) where {T}
    for x in xs
        acc(x)
    end
    acc
end

#=
reference for AccExpWtMean, AccExpWtMeanVar

Incremental calculation of weighted mean and variance
by Tony Finch
=#

mutable struct AccExpWtMean{T} <: Accumulator{T}
    nobs::Int
    alpha::T
    expwtmean::T
end

AccExpWtMean(::Type{T}=Float64; alpha::T=T(0.5)) where {T} = new{T}(0, T(alpha), zero(T))

(acc::AccExpWtMean{T})() where {T} = acc.expwtmean

function (acc::AccExpWtMean{T})(x) where {T}
    acc.nobs += 1
    acc.expwtmean += acc.alpha * (x - acc.expwtmean)
    acc
end

function (acc::AccExpWtMean{T})(xs::Seq{T}) where {T}
    for x in xs
        acc(x)
    end
    acc
end

mutable struct AccExpWtMeanVar{T} <: Accumulator{T}
    nobs::Int
    alpha::T
    expwtmean::T
    expwtsvar::T
end

AccExpWtMeanVar(::Type{T}=Float64; alpha::T=T(0.5)) where {T} = AccExpWtMeanVar(0, alpha, zero(T), zero(T))

function(acc::AccExpWtMeanVar{T})() where {T}
    unbiased_expwtvar = acc.expwtsvar / (acc.nobs - 1)
    (expwt_mean=acc.expwtmean, expwt_var=unbiased_expwtvar)
end

function (acc::AccExpWtMeanVar{T})(x) where {T}
    acc.nobs += 1
    diff = x - acc.expwtmean
    incr = acc.alpha * diff
    acc.expwtmean += acc.alpha * (x - acc.expwtmean)
    acc.expwtsvar = (one(T) - acc.alpha) * (acc.expwtsvar + diff * incr)
    acc
end

function (acc::AccExpWtMeanVar{T})(xs::Seq{T}) where {T}
    for x in xs
        acc(x)
    end
    acc
end

mutable struct AccExpWtMeanStd{T} <: Accumulator{T}
    nobs::Int
    alpha::T
    expwtmean::T
    expwtsvar::T
end

AccExpWtMeanStd(::Type{T}=Float64; alpha::T=T(0.5)) where {T} = AccExpWtMeanStd(0, alpha, zero(T), zero(T))

function(acc::AccExpWtMeanStd{T})() where {T}
    unbiased_expwtstd = sqrt(acc.expwtsvar / (acc.nobs - 1))
    (expwt_mean=acc.expwtmean, expwt_std=unbiased_expwtstd)
end

function (acc::AccExpWtMeanStd{T})(x) where {T}
    acc.nobs += 1
    diff = x - acc.expwtmean
    incr = acc.alpha * diff
    acc.expwtmean += acc.alpha * (x - acc.expwtmean)
    acc.expwtsvar = (one(T) - acc.alpha) * (acc.expwtsvar + diff * incr)
    acc
end

function (acc::AccExpWtMeanStd{T})(xs::Seq{T}) where {T}
    for x in xs
        acc(x)
    end
    acc
end

Base.length(@nospecialize acc::Accumulator) = acc.nobs
StatsBase.nobs(@nospecialize acc::Accumulator) = acc.nobs

for (F,A) in ((:(Base.minimum), :AccMinimum), (:(Base.maximum), :AccMaximum), (:(Base.extrema), :AccExtrema),
              (:(Base.sum), :AccSum), (:(Base.prod), :AccProd),
              (:(StatsBase.mean), :AccMean), (:(StatsBase.geomean), :AccGeoMean), (:(StatsBase.harmmean), :AccHarmMean))
     @eval $F(acc::$A) = acc()
end

Base.minimum(acc::AccExtrema) = acc.min
Base.maximum(acc::AccExtrema) = acc.max
midrange(acc::AccExtrema) = (acc.min / 2) + (acc.max / 2)
proportionalrange(acc::AccExtrema, proportion) = (acc.min * proportion) + (acc.max * (1 - proportion))

nminima(acc::AccMinimum) = acc.nmin
nminima(acc::AccExtrema) = acc.nmin
nmaxima(acc::AccMinimum) = acc.nmax
nmaxima(acc::AccExtrema) = acc.nmax

StatsBase.mean(acc::AccMeanVar) = acc.mean
StatsBase.var(acc::AccMeanVar) = acc.svar / (acc.nobs - 1)
StatsBase.std(acc::AccMeanVar) = sqrt(acc.svar / (acc.nobs - 1))

StatsBase.mean(acc::AccStats{T}) where {T} = T(acc.m1)
StatsBase.var(acc::AccStats{T}) where {T} = T(acc.m2 / (acc.nobs - 1))
StatsBase.std(acc::AccStats{T}) where {T} = T(sqrt(var(acc)))
StatsBase.skewness(acc::AccStats{T}) where {T} = T(sqrt(acc.nobs) * acc.m3 / (acc.m2 * sqrt(acc.m2)))
StatsBase.kurtosis(acc::AccStats{T}) where {T} = T( ((acc.nobs * acc.m4) / (acc.m2^2)) - 3)
                                        
           
