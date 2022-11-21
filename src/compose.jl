#=
     AccCount, 
     AccMin, AccMax, AccExtrema, 
     AccSum, AccProd,
     AccMean, AccGeometricMean, AccHarmonicMean,
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

nminima(acc::AccMinimum) = acc.nmin
nminima(acc::AccExtrema) = acc.nmin
nmaxima(acc::AccMinimum) = acc.nmax
nmaxima(acc::AccExtrema) = acc.nmax

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


            
          
          
          
          
          
          
          
          
          
          
          
          
          
          
mutable struct AccMax{T} <: Accumulator{T}
    max::T
end

function AccMax(::Type{T}=DefaultFloat) where {T}
     AccMax{T}(typemin(T))
end

function (acc::AccMax{T})() where {T}
    acc.max
end
     
function (acc::AccMax{T})(x) where {T}
    if x > acc.max
        acc.max = x
    end
    acc
end

function (acc::AccMax{T})(xs::A) where {T, A<:AbstractVector{T}}
    acc(vmaximum(xs))
end

function (acc::AccCount{T})(xs::NTuple{N,T}) where {T, N}
    acc(maximum(xs))
end



# Max
     
mutable struct AccMin{T} <: Accumulator{T}
    n::Int
    min::T
    AccMin(::Type{T}=Float64) where {T} =
        (T <: Integer) ? new{T}(typemax(T)) : new{T}(0, floatmax(T))
end

(acc::AccMin{T})() where {T} = (acc.min)
(acc::AccMin{T})(x) where {T} = (acc.n +=1; acc.min = ifelse(x < acc.min, T(x), acc.min); acc)
(acc::AccMin{T})(xs::Seq) where {T} = (acc.n += length(xs); x = T(vminimum(xs)); acc.min = ifelse(x < acc.min, x, acc.min); acc)

mutable struct AccMax{T} <: Accumu
acc_nmin(acc::AccExtrema{T}) where {T} = acc.nmin
acc_nmax(acc::AccExtrema{T}) where {T} = acc.nmax
acc_midrange(acc::AccExtrema{T}) where {T} = (acc.max / 2) + (acc.min / 2)
                                                                      
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

