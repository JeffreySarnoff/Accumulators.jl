#=
   const AccStruct =
   [ AccCount, AccMin, AccMax, AccMinAbs, AccMaxAbs,
     AccExtrema, AccExtremaAbs, 
     AccSum, AccProd,
     AccMean, AccMeanAbs, AccMeanAbs2, AccGeometricMean, AccHarmonicMean,
     AccMeanVar,
     AccExpWMean, AccExpWMeanVar]
=#

mutable struct AccCount{T} <: Accumulator{T}
    n::T
    AccCount(::Type{T}=Int64) where {T} = new{T}(zero(T))
end

(acc::AccCount{T})() where {T} = (acc.n)
(acc::AccCount{T})(x) where {T} = (acc.n += one(T); nothing)

mutable struct AccMin{T} <: Accumulator{T}
    min::T
    AccMin(::Type{T}=Float64) where {T} =
        (T <: Integer) ? new{T}(typemax(T)) : new{T}(floatmax(T))
end

(acc::AccMin{T})() where {T} = (acc.min)
(acc::AccMin{T})(x) where {T} = (acc.min = ifelse(x < acc.min, T(x), acc.min); nothing)

mutable struct AccMax{T} <: Accumulator{T}
    max::T
    AccMax(::Type{T}=Float64) where {T} =
        (T <: Integer) ? new{T}(typemin(T)) : new{T}(floatmin(T))
end

(acc::AccMax{T})() where {T} = (acc.max)
(acc::AccMax{T})(x) where {T} = (acc.max = ifelse(acc.max < x, T(x), acc.max); nothing)

mutable struct AccExtrema{T} <: Accumulator{T}
    min::T
    max::T
    AccExtrema(::Type{T}=Float64) where {T} =
        (T <: Integer) ? new{T}(typemax(T), typemin(T)) : new{T}(floatmax(T), floatmin(T))
end

(acc::AccExtrema{T})() where {T} = (acc.min, acc.max)
(acc::AccExtrema{T})(x) where {T} = 
    (acc.min = ifelse(x < acc.min, T(x), acc.min);
     acc.max = ifelse(acc.max < x, T(x), acc.max); nothing)

acc_max(acc::AccExtrema{T}) where {T} = acc.max
acc_min(acc::AccExtrema{T}) where {T} = acc.min
acc_midrange(acc::AccExtrema{T}) where {T} = (acc.max / 2) + (acc.min / 2)

mutable struct AccSum{T} <: Accumulator{T}
    sum::T
    AccSum(::Type{T}=Float64) where {T} = new{T}(zero(T),)
end

(acc::AccSum{T})() where {T} = (acc.sum)
(acc::AccSum{T})(x) where {T} = (acc.sum += x; nothing)

mutable struct AccProd{T} <: Accumulator{T}
    prod::T
    AccProd(::Type{T}=Float64) where {T} = new{T}(one(T),)
end

(acc::AccProd{T})() where {T} = (acc.prod)
(acc::AccProd{T})(x) where {T} = (acc.prod *= x; nothing)

mutable struct AccMean{T} <: Accumulator{T}
    n::Int
    mean::T
    AccMean(::Type{T}=Float64) where {T} = new{T}(0, zero(T))
end

(acc::AccMean{T})() where {T} = (acc.mean)
(acc::AccMean{T})(x) where {T} =
    (acc.n += 1; acc.mean += (x - acc.mean) / acc.n; nothing)

# geometric mean (of abs(xs))
# see https://github.com/stdlib-js/stats/blob/main/incr/gmean/lib/main.js
mutable struct AccGeometricMean{T} <: Accumulator{T}
    n::Int
    sumlog::T
    AccGeometricMean(::Type{T}=Float64) where {T} = new{T}(0, zero(T))
end

(acc::AccGeometricMean{T})() where {T} = (iszero(acc.n) ? one(T) : exp(acc.sumlog / acc.n))
(acc::AccGeometricMean{T})(x) where {T} = (acc.n += 1; acc.sumlog += log(abs(x)); nothing)

# harmonic mean
# see https://github.com/stdlib-js/stats/blob/main/incr/hmean/lib/main.js
mutable struct AccHarmonicMean{T} <: Accumulator{T}
    n::Int
    hmean::T
    AccHarmonicMean(::Type{T}=Float64) where {T} = new{T}(0, zero(T))
end

(acc::AccHarmonicMean{T})() where {T} = (iszero(acc.n) ? one(T) : acc.n / acc.hmean)
(acc::AccHarmonicMean{T})(x) where {T} = (acc.n += 1; acc.hmean += (one(T) / x); nothing)

# Unbiased Sample Variation (with Mean)
# see https://www.johndcook.com/blog/standard_deviation/

mutable struct AccMeanVar{T} <: Accumulator{T}
    n::Int
    mean::T
    svar::T
    AccMeanVar(::Type{T}=Float64) where {T} = new{T}(0, zero(T), zero(T))
end

function (acc::AccMeanVar{T})() where {T}
    unbiased_var = acc.svar / (acc.n - 1)
    (acc.mean, unbiased_var)
end

function (acc::AccMeanVar{T})(x) where {T}
    acc.n += 1
    if !iszero(acc.n)
        oldmean = acc.mean
        acc.mean = oldmean + (x - oldmean) / acc.n
        acc.svar = acc.svar + (x - oldmean) * (x - acc.mean)
    else
        acc.mean = x  # svar is zero already
    end
    nothing
end

# see https://www.johndcook.com/blog/skewness_kurtosis/

mutable struct AccStats{T} <: Accumulator{T}
    n::Int
    m1::T
    m2::T
    m3::T
    m4::T
    AccStats(::Type{T}=Float64) where {T} = new{T}(0, zero(T), zero(T), zero(T), zero(T))
end

function (acc::AccStats{T})() where {T}
    (acc.acc_count(), acc.acc_mean(), acc.acc_var(), acc.acc_std(), acc.acc_skew(), acc.acc_kurt())
end

function (acc::AccStats{T})(x) where {T}
    n1 = acc.n
    acc.n += 1
    delta = x - acc.m1
    delta_n = delta / acc.n
    delta_n2 = delta_n^2
    term1 = delta * delta_n * n1
    acc.m1 += delta_n
    acc.m4 += term1 * delta_n2 * (acc.n^2 - 3*acc.n + 3) + 
              6 * delta_2 * acc.m2 - 
              4 * delta_n * acc.m3
    acc.m3 += term1 * delta_n * (n - 2) - 3 * delta_n * acc.m2
    acc.m2 += term1
    nothing
end

acc_count(acc::AccStats{T}) where {T} = acc.n
acc_mean(acc::AccStats{T}) where {T} = T(acc.m1)
acc_var(acc::AccStats{T}) where {T} = T(acc.m2 / (acc.n - 1))
acc_std(acc::AccStats{T}) where {T} = T(sqrt(var(x)))
acc_skew(acc::AccStats{T}) where {T} = T(sqrt(acc.n) * acc.m3 / (acc.m2 * sqrt(acc.m2)))
acc_kurt(acc::AccStats{T}) where {T} = T((acc.n * acc.m4) / (acc.m2^2) - 3)

#=
reference for AccExpWtMean, AccExpWtMeanVar

Incremental calculation of weighted mean and variance
by Tony Finch
=#
                           
mutable struct AccExpWtMean{T} <: Accumulator{T}
    n::Int
    alpha::T
    mean::T
    AccExpWtMean(::Type{T}=Float64; alpha::T=0.5) where {T} = new{T}(0, alpha, zero(T))
end

(acc::AccExpWtMean{T})() where {T} = (acc.mean)
(acc::AccExpWtMean{T})(x) where {T} = (acc.n += 1; acc.mean += acc.alpha * (x - acc.mean))

mutable struct AccExpWtMeanVar{T} <: Accumulator{T}
    n::Int
    alpha::T
    mean::T
    svar::T
    AccExpWtMeanVar(::Type{T}=Float64; alpha::T=0.5) where {T} = new{T}(0, alpha, zero(T), zero(T))
end

function(acc::AccExpWtMeanVar{T})() where {T}
    unbiased_var = acc.svar / (acc.n = 1)
    (acc.mean, unbiased_var)
end

function (acc::AccExpWtMeanVar{T})(x) where {T}
    acc.n += 1
    diff = x - acc.mean
    incr = acc.alpha * diff
    acc.mean += acc.alpha * (x - acc.mean)
    acc.svar = (one(T) - acc.alpha) * (acc.svar + diff * incr)
end

# other derived

acc_count(acc::Accumulator) = acc.n

acc_mean(acc::AccMeanVar{T}) where {T} = acc.mean
acc_var(acc::AccMeanVar{T}) where {T} = acc.svar / (acc.n - 1)
acc_std(acc::AccMeanVar{T}) where {T} = sqrt(acc.svar / (acc.n - 1))

acc_mean(acc::AccExpWtMean{T}) where {T} = acc.mean

acc_mean(acc::AccExpWtMeanVar{T}) where {T} = acc.mean
acc_var(acc::AccExpWtMeanVar{T}) where {T} = acc.svar / (acc.n - 1)
acc_std(acc::AccExpWtMeanVar{T}) where {T} = sqrt(acc_var(acc))


