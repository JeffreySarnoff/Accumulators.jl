#=
     AccCount, 
     AccMinimum, AccMaximum, AccExtrema, 
     AccSum, AccProd,
     AccMean, AccGeoMean, AccHarmoMean,
     AccMeanVar,
     AccExpWtMean, AccExpWtMeanVar
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

DefaultFloat = Float64

StatsBase.nobs(acc::A) where {T, A<:Accumulator{T}} = acc.nobs

Base.minimum(acc::AccMinimum) = acc.min
Base.minimum(acc::AccExtrema) = acc.min
Base.maximum(acc::AccMaximum) = acc.max
Base.maximum(acc::AccExtrema) = acc.max

midrange(acc::AccExtrema) = (acc.min / 2) + (acc.max / 2)
proportionalrange(acc::AccExtrema, proportion) = (acc.min * proportion) + (acc.max * (1 - proportion))

nminima(acc::AccMinimum) = acc.nmin
nminima(acc::AccExtrema) = acc.nmin
nmaxima(acc::AccMinimum) = acc.nmax
nmaxima(acc::AccExtrema) = acc.nmax

Base.sum(acc::AccSum) = acc.sum
Base.prod(acc::AccProd) = acc.prod

StatsBase.mean(acc::AccMean) = acc.mean

# Count

mutable struct AccCount{T} <: Accumulator{T}
    nobs::T
end

function AccCount(::Type{T}=Int64) where {T}
     AccCount{T}(zero(T))
end

function (acc::AccCount{T})() where {T}
    acc.nobs
end
     
function (acc::AccCount{T})(x) where {T}
    acc.nobs += one(T)
    acc
end

function (acc::AccCount{T})(xs::A) where {T, A<:AbstractVector{T}}
    acc.nobs += length(xs)
    acc
end

function (acc::AccCount{T})(xs::NTuple{N,T}) where {T, N}
    acc.nobs += N
    acc
end

# Min

mutable struct AccMinimum{T} <: Accumulator{T}
    nobs::Int
    nmin::Int
    min::T
end

function AccMinimum(::Type{T}=DefaultFloat) where {T}
     AccMinimum{T}(0, 0, typemax(T))
end

function (acc::AccMinimum{T})() where {T}
    acc.min
end

function (acc::AccMinimum{T})(x) where {T}
    acc.nobs += 1
    if x < acc.min
        acc.nmin += 1
        acc.min = x
    end
    acc
end

function (acc::AccMinimum{T})(xs::A) where {T, A<:AbstractVector{T}}
    acc.nobs += length(xs)     
    x = vminimum(xs)
    if x < acc.min
        acc.nmin += 1
        acc.min = x
    end
    acc
end

function (acc::AccMinimum{T})(xs::NTuple{N,T}) where {T, N}
    acc.nobs += N
    x = minimum(xs)
    if x < acc.min
        acc.nmin += 1
        acc.min = x
    end
    acc
end

# Max

mutable struct AccMaximum{T} <: Accumulator{T}
    nobs::Int
    nmax::Int
    max::T
end

function AccMaximum(::Type{T}=DefaultFloat) where {T}
     AccMaximum{T}(0, 0, typemin(T))
end

function (acc::AccMaximum{T})() where {T}
    acc.max
end

function (acc::AccMaximum{T})(x) where {T}
    acc.nobs += 1
    if x > acc.max
        acc.nmax += 1
        acc.max = x
    end
    acc
end

function (acc::AccMaximum{T})(xs::A) where {T, A<:AbstractVector{T}}
    acc.nobs += length(xs)     
    x = vmaximum(xs)
    if x > acc.max
        acc.nmax += 1
        acc.max = x
    end
    acc
end

function (acc::AccMaximum{T})(xs::NTuple{N,T}) where {T, N}
    acc.nobs += N
    x = maximum(xs)
    if x > acc.max
        acc.nmax += 1
        acc.max = x
    end
    acc
end

# Extrema

mutable struct AccExtrema{T} <: Accumulator{T}
    nobs::Int
    nmin::Int
    nmax::Int
    min::T
    max::T
end

function AccExtrema(::Type{T}=DefaultFloat) where {T}
     AccExtrema{T}(0, 0, 0, typemax(T), typemin(T))
end

function (acc::AccExtrema{T})() where {T}
    (acc.min, acc.max)
end

function (acc::AccExtrema{T})(x) where {T}
    acc.nobs += 1
    if x < acc.min
        acc.nmin += 1
        acc.min = x
    end
    acc
end

function (acc::AccExtrema{T})(xs::A) where {T, A<:AbstractVector{T}}
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

function (acc::AccExtrema{T})(xs::NTuple{N,T}) where {T, N}
    acc.nobs += N
    mn, mx = extrema(xs)
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
    nobs::Int
    sum::T
end

function AccSum(::Type{T}=DefaultFloat) where {T}
     AccSum{T}(0, zero(T))
end

function (acc::AccSum{T})() where {T}
    acc.sum
end

function (acc::AccSum{T})(x) where {T}
    acc.nobs += 1
    acc.sum += x
    acc
end

function (acc::AccSum{T})(xs::A) where {T, A<:AbstractVector{T}}
    acc.nobs += length(xs)     
    x = vsum(xs)
    acc.sum += x
    acc
end

function (acc::AccSum{T})(xs::NTuple{N,T}) where {T, N}
    acc.nobs += N
    x = sum(xs)
    acc.sum += x
    acc
end

# Prod

mutable struct AccProd{T} <: Accumulator{T}
    nobs::Int
    prod::T
end

function AccProd(::Type{T}=DefaultFloat) where {T}
     AccProd{T}(0, one(T))
end

function (acc::AccProd{T})() where {T}
    acc.prod
end

function (acc::AccProd{T})(x) where {T}
    acc.nobs += 1
    acc.prod *= x
    acc
end

function (acc::AccProd{T})(xs::A) where {T, A<:AbstractVector{T}}
    acc.nobs += length(xs)     
    x = vprod(xs)
    acc.prod *= x
    acc
end

function (acc::AccProd{T})(xs::NTuple{N,T}) where {T, N}
    acc.nobs += N
    x = prod(xs)
    acc.prod *= x
    acc
end

# Mean

mutable struct AccMean{T} <: Accumulator{T}
    nobs::Int
    mean::T
end

function AccMean(::Type{T}=DefaultFloat) where {T}
     AccMean{T}(0, zero(T))
end

function (acc::AccMean{T})() where {T}
    acc.mean
end

function (acc::AccMean{T})(x) where {T}
    acc.nobs += 1
    acc.mean += (x - acc.mean) / acc.nobs
    acc
end

function (acc::AccMean{T})(xs::A) where {T, A<:AbstractVector{T}}
    acc.nobs += length(xs)     
    xmean = vmean(xs)
    acc.mean += (xmean - acc.mean) / acc.nobs
    acc
end

function (acc::AccMean{T})(xs::NTuple{N,T}) where {T, N}
    acc.nobs += N
    xmean = mean(xs)
    acc.mean += (xmean - acc.mean) / acc.nobs
    acc
end

# GeoMean

logabs(x) = log(abs(x))
                                   
mutable struct AccGeoMean{T} <: Accumulator{T}
    nobs::Int
    sumlog::T
end

function AccGeoMean(::Type{T}=DefaultFloat) where {T}
     AccGeoMean{T}(0, one(T))
end

function (acc::AccGeoMean{T})() where {T}
    n = ifelse(acc.nobs === 0 ? 1 : acc.nobs)
    exp(acc.sumlog / n)
end

function (acc::AccGeoMean{T})(x) where {T}
    acc.nobs += 1
    acc.sumlog += logabs(x)
    acc
end

function (acc::AccGeoMean{T})(xs::A) where {T, A<:AbstractVector{T}}
    acc.nobs += length(xs)
    acc.sumlog += sum(map(logabs, xs))
    acc
end

function (acc::AccGeoMean{T})(xs::NTuple{N,T}) where {T, N}
    acc.nobs += N
    acc.sumlog += sum(map(logabs, xs))
    acc
end


#
                                                                      
acc_mean(acc::AccStats{T}) where {T} = T(acc.m1)
acc_var(acc::AccStats{T}) where {T} = T(acc.m2 / (acc.n - 1))
acc_std(acc::AccStats{T}) where {T} = T(sqrt(acc_var(acc)))
acc_skew(acc::AccStats{T}) where {T} = T(sqrt(acc.n) * acc.m3 / (acc.m2 * sqrt(acc.m2)))
acc_kurt(acc::AccStats{T}) where {T} = T((acc.n * acc.m4) / (acc.m2^2) - 3)
                                                                      
acc_mean(acc::AccMean{T}) where {T} = acc.mean
acc_mean(acc::AccGeometricMean{T}) where {T} = acc()
acc_mean(acc::AccHarmonicMean{T}) where {T} = acc()

acc_mean(acc::AccMeanVar{T}) where {T} = acc.mean
acc_var(acc::AccMeanVar{T}) where {T} = acc.svar / (acc.n - 1)
acc_std(acc::AccMeanVar{T}) where {T} = sqrt(acc.svar / (acc.n - 1))

acc_mean(acc::AccExpWtMean{T}) where {T} = acc.mean

acc_mean(acc::AccExpWtMeanVar{T}) where {T} = acc.mean
acc_var(acc::AccExpWtMeanVar{T}) where {T} = acc.svar / (acc.n - 1)
acc_std(acc::AccExpWtMeanVar{T}) where {T} = sqrt(acc_var(acc))

