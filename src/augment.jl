mutable struct AccumCount{T,FN} <:  Accumulator{T}
    n::T
    const fn::FN
    AccumCount(::Type{T}=Int64, fn::FN=identity) where {T,FN} = new{T,FN}(zero(T), fn)
end

(accum::AccumCount{T,FN})() where {T,FN} = (accum.n)
(accum::AccumCount{T,FN})(x) where {T,FN} = (accum.n += accum.fn(one(T)); accum.n)

mutable struct AccumMin{T,FN} <:  Accumulator{T}
    min::T
    const fn::FN
    AccumMin(::Type{T}=Float64, fn::FN=identity) where {T,FN} =
        (T <: Integer) ? new{T,FN}(typemax(T), fn) : new{T,FN}(floatmax(T), fn)
end

(accum::AccumMin{T,FN})() where {T,FN} = (accum.min)
(accum::AccumMin{T,FN})(x) where {T,FN} = (accum.min = ifelse(x < accum.min, accum.fn(T(x)), accum.min); accum.min)

mutable struct AccumMax{T,FN} <:  Accumulator{T}
    max::T
    const fn::FN
    AccumMax(::Type{T}=Float64, fn::FN=identity) where {T,FN} =
        (T <: Integer) ? new{T,FN}(typemin(T), fn) : new{T,FN}(floatmin(T), fn)
end

(accum::AccumMax{T,FN})() where {T,FN} = (accum.max)
(accum::AccumMax{T,FN})(x) where {T,FN} = (accum.max = ifelse(accum.max < x, accum.fn(T(x)), accum.max); accum.max)

mutable struct AccumExtrema{T,FN} <:  Accumulator{T}
    min::T
    max::T
    const fn::FN
    AccumExtrema(::Type{T}=Float64, fn::FN=identity) where {T,FN} =
        (T <: Integer) ? new{T,FN}(typemax(T), typemin(T), fn) : new{T,FN}(floatmax(T), floatmin(T), fn)
end

(accum::AccumExtrema{T,FN})() where {T,FN} = (accum.min, accum.max)
(accum::AccumExtrema{T,FN})(x) where {T,FN} = 
    (accum.min = ifelse(x < accum.min, accum.fn(T(x)), accum.min); 
     accum.max = ifelse(accum.max < x, fn(T(x)), accum.max); (accum.min, accum.max))

acc_min(accum::AccumExtrema{T,FN}) where {T,FN} = accum.min
acc_max(accum::AccumExtrema{T,FN}) where {T,FN} = accum.max
acc_midrange(accum::AccumExtrema{T,FN}) where {T,FN} = (accum.max / 2) + (accum.min / 2)

mutable struct AccumSum{T,FN} <:  Accumulator{T}
    sum::T
    const fn::FN
    AccumSum(::Type{T}=Float64, fn::FN=identity) where {T,FN} = new{T,FN}(zero(T), fn)
end

(accum::AccumSum{T,FN})() where {T,FN} = (accum.sum)
(accum::AccumSum{T,FN})(x) where {T,FN} = (accum.sum += accum.fn(x); accum.sum)

mutable struct AccumProd{T,FN} <:  Accumulator{T}
    prod::T
    const fn::FN
    AccumProd(::Type{T}=Float64, fn::FN=identity) where {T,FN} = new{T,FN}(one(T), fn)
end

(accum::AccumProd{T,FN})() where {T,FN} = (accum.prod)
(accum::AccumProd{T,FN})(x) where {T,FN} = (accum.prod *= accum.fn(x); accum.prod)

mutable struct AccumMean{T,FN} <:  Accumulator{T}
    n::Int
    mean::T
    const fn::FN
    AccumMean(::Type{T}=Float64, fn::FN=identity) where {T,FN} = new{T,FN}(0, zero(T), fn)
end

(accum::AccumMean{T,FN})() where {T,FN} = (accum.mean)
(accum::AccumMean{T,FN})(x) where {T,FN} =
    (accum.n += 1; accum.mean += (accum.fn(x) - accum.mean) / accum.n; accum.mean)

# geometric mean (of abs(xs))
# see https://github.com/stdlib-js/stats/blob/main/incr/gmean/lib/main.js
mutable struct AccumGeometricMean{T,FN} <:  Accumulator{T}
    n::Int
    sumlog::T
    const fn::FN
    AccumGeometricMean(::Type{T}=Float64, fn::FN=identity) where {T,FN} = new{T,FN}(0, zero(T), fn)
end

(accum::AccumGeometricMean{T,FN})() where {T,FN} = (iszero(accum.n) ? one(T) : exp(accum.sumlog / accum.n))
(accum::AccumGeometricMean{T,FN})(x) where {T,FN} = (accum.n += 1; accum.sumlog += log(abs(accum.fn(x))); accum())

# harmonic mean
# see https://github.com/stdlib-js/stats/blob/main/incr/hmean/lib/main.js
mutable struct AccumHarmonicMean{T,FN} <:  Accumulator{T}
    n::Int
    hmean::T
    const fn::FN
    AccumHarmonicMean(::Type{T}=Float64, fn::FN=identity) where {T,FN} = new{T,FN}(0, zero(T), fn)
end

(accum::AccumHarmonicMean{T})() where {T} = (iszero(accum.n) ? one(T) : accum.n / accum.hmean)
(accum::AccumHarmonicMean{T})(x) where {T} = (accum.n += 1; accum.hmean += (one(T) / accum.fn(x)); accum())

# Unbiased Sample Variation (with Mean)
# see https://www.johndcook.com/blog/standard_deviation/

mutable struct AccumMeanVar{T,FN} <: Accumulator{T}
    n::Int
    mean::T
    svar::T
    const fn::FN
    AccumMeanVar(::Type{T}=Float64, fn::FN=identity) where {T,FN} = new{T,FN}(0, zero(T), zero(T), fn)
end

function (accum::AccumMeanVar{T,FN})() where {T,FN}
    unbiased_var = accum.svar / (accum.n - 1)
    (accum.mean, unbiased_var)
end

function (accum::AccumMeanVar{T,FN})(x) where {T,FN}
    if !iszero(accum.n)
        x = fn(x)
        oldmean = accum.mean
        accum.mean = accum.fn(oldmean + (x - oldmean) / accum.n)
        accum.svar = accum.svar + (x - oldmean) * (x - accum.mean)
    else
        accum.mean = accum.fn(x)  # svar is zero already
    end
    accum()
end

# see https://www.johndcook.com/blog/skewness_kurtosis/

mutable struct AccumStats{T,FN} <: Accumulator{T}
    n::Int
    m1::T
    m2::T
    m3::T
    m4::T
    const fn::FN
    AccumStats(::Type{T}=Float64, fn::FN=identity) where {T} = new{T}(0, zero(T), zero(T), zero(T), zero(T), fn)
end

function (acc::AccumStats{T})() where {T}
    (acc.acc_count(), acc.acc_mean(), acc.acc_var(), acc.acc_std(), acc.acc_skew(), acc.acc_kurt())
end

function (acc::AccumStats{T})(x) where {T}
    n1 = acc.n
    acc.n += 1
    delta = acc.fn(x) - acc.m1
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
