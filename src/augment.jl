mutable struct AccumCount{T,F} <: Accumulator{T}
    n::T
    const fn::F
    AccumCount(::Type{T}=Int64, fn::F=identity) where {T,F} = new{T,F}(zero(T), fn)
end

(acc::AccumCount{T,F})() where {T,F} = acc.n
(acc::AccumCount{T,F})(x) where {T,F} = (acc.n += fn(one(T)); acc)
(acc::AccumCount{T,F})(xs::Seq) where {T,F} = (acc.n += fn(T(length(xs))); acc)

mutable struct AccumMin{T,F} <: Accumulator{T}
    n::Int
    min::T
    const fn::F
    AccumMin(::Type{T}=Float64, fn::F=identity) where {T,F} =
        (T <: Integer) ? new{T,F}(0, typemax(T), fn) : new{T,F}(0, floatmax(T), fn)
end

(acc::AccumMin{T,F})() where {T,F} = acc.min
(acc::AccumMin{T,F})(x) where {T,F} = (acc.n += 1; acc.min = ifelse(x < acc.min, fn(T(x)), acc.min); acc)
(acc::AccumMin{T,F})(xs::Seq) where {T,F} = (acc.n += length(xs); x = T(vminimum(xs)); acc.min = ifelse(x < acc.min, fn(x), acc.min); acc)

mutable struct AccumMax{T,F} <: Accumulator{T}
    n::Int
    max::T
    const fn::F
    AccumMax(::Type{T}=Float64, fn::F=identity) where {T,F} =
        (T <: Integer) ? new{T,F}(0, typemin(T), fn) : new{T,F}(0, floatmin(T), fn)
end

(acc::AccumMax{T,F})() where {T,F} = acc.max
(acc::AccumMax{T,F})(x) where {T,F} = (acc.n += 1; acc.min = ifelse(x < acc.min, fn(T(x)), acc.max); acc)
(acc::AccumMax{T,F})(xs::Seq) where {T,F} = (acc.n += length(xs); x = T(vmaximum(xs)); acc.max = ifelse(x > acc.max, fn(x), acc.max); acc)

mutable struct AccumExtrema{T,F} <: Accumulator{T}
    n::Int
    nmin::Int
    nmax::Int
    min::T
    max::T
    const fn::F
    AccumExtrema(::Type{T}=Float64, fn::F=identity) where {T,F} =
        (T <: Integer) ? new{T, F}(0, 0, 0, typemax(T), typemin(T), fn) : new{T, F}(0, 0, 0, floatmax(T), floatmin(T), fn)
end

(acc::AccumExtrema{T,F})() where {T,F} = (acc.min, acc.max)

function (acc::AccumExtrema{T,F})(x) where {T,F}
    acc.n += 1
    xx = fn(T(x))
    if xx < acc.min
       acc.nmin += 1
       acc.min = xx
    end
    if xx > acc.max
       acc.nmax += 1
       acc.max = xx
    end
    acc
 end
 
function (acc::AccumExtrema{T,F})(xs::Seq) where {T,F}
    acc.n += length(xs)
    xxs = map(fn, xs)
    mn, mx = map(T, vextrema(xxs))
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

mutable struct AccumSum{T,FN} <:  Accumulator{T}
    n::Int
    sum::T
    const fn::FN
    AccumSum(::Type{T}=Float64, fn::FN=identity) where {T,FN} = new{T, FN}(0, zero(T), fn)
end

(accum::AccumSum{T,FN})() where {T,FN} = accum.sum
(accum::AccumSum{T,FN})(x) where {T,FN} = (accum.n += 1; accum.sum += accum.fn(x); accum)
(accum::AccumSum{T,FN})(xs::Seq) where {T,FN} = (accum.n += length(xs); xxs = map(accum.fn, xs); x = T(vsum(xxs)); accum.sum += x; accum)

mutable struct AccumProd{T,FN} <:  Accumulator{T}
    n::Int
    prod::T
    const fn::FN
    AccumProd(::Type{T}=Float64, fn::FN=identity) where {T,FN} = new{T, FN}(0, one(T), fn)
end

(accum::AccumProd{T,FN})() where {T,FN} = accum.prod
(accum::AccumProd{T,FN})(x) where {T,FN} = (accum.n += 1; accum.prod *= accum.fn(x); accum)
(accum::AccumProd{T,FN})(xs::Seq) where {T,FN} = (accum.n += length(xs); xxs = map(accum.fn, xs); x = T(prod(xxs)); accum.prod += x; accum)

mutable struct AccumMean{T,FN} <:  Accumulator{T}
    n::Int
    mean::T
    const fn::FN
    AccumMean(::Type{T}=Float64, fn::FN=identity) where {T,FN} = new{T,FN}(0, zero(T), fn)
end

(accum::AccumMean{T,FN})() where {T,FN} = (accum.mean)
(accum::AccumMean{T,FN})(x) where {T,FN} =
    (accum.n += 1; accum.mean += (accum.fn(x) - accum.mean) / accum.n; accum)
(accum::AccumMean{T,FN})(xs::Seq) where {T,FN} = (accum.n += length(xs); xxs = map(accum.fn, xs); x = T(vmean(xxs)); accum.sum += x; accum)

# geometric mean (of abs(xs))
# see https://github.com/stdlib-js/stats/blob/main/incr/gmean/lib/main.js
mutable struct AccumGeometricMean{T,FN} <:  Accumulator{T}
    n::Int
    sumlog::T
    const fn::FN
    AccumGeometricMean(::Type{T}=Float64, fn::FN=identity) where {T,FN} = new{T,FN}(0, zero(T), fn)
end

(accum::AccumGeometricMean{T,FN})() where {T,FN} = (iszero(accum.n) ? one(T) : exp(accum.sumlog / accum.n))
(accum::AccumGeometricMean{T,FN})(x) where {T,FN} = (accum.n += 1; accum.sumlog += log(abs(accum.fn(x))); accum)
(accum::AccumGeometricMean{T,FN})(xs::Seq) where {T,FN} = (accum.n += length(xs); xxs = map(accum, xs); accum)

# harmonic mean
# see https://github.com/stdlib-js/stats/blob/main/incr/hmean/lib/main.js
mutable struct AccumHarmonicMean{T,FN} <:  Accumulator{T}
    n::Int
    hmean::T
    const fn::FN
    AccumHarmonicMean(::Type{T}=Float64, fn::FN=identity) where {T,FN} = new{T,FN}(0, zero(T), fn)
end

(accum::AccumHarmonicMean{T})() where {T} = (iszero(accum.n) ? one(T) : accum.n / accum.hmean)
(accum::AccumHarmonicMean{T})(x) where {T} = (accum.n += 1; accum.hmean += (one(T) / accum.fn(x)); accum)
(accum::AccumHarmonicMean{T,FN})(xs::Seq) where {T,FN} = (accum.n += length(xs); xxs = map(accum, xs); accum)

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
    accum
end

(accum::AccumMeanVar{T,FN})(xs::Seq) where {T,FN} = (accum.n += length(xs); xxs = map(accum, xs); accum)

# see https://www.johndcook.com/blog/skewness_kurtosis/

mutable struct AccumStats{T,FN} <: Accumulator{T}
    n::Int
    m1::T
    m2::T
    m3::T
    m4::T
    const fn::FN
    AccumStats(::Type{T}=Float64, fn::FN=identity) where {T,FN} = new{T,FN}(0, zero(T), zero(T), zero(T), zero(T), fn)
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
    acc
end

(accum::AccumStats{T,FN})(xs::Seq) where {T,FN} = (accum.n += length(xs); xxs = map(accum, xs); accum)

#=
reference for AccExpWtMean, AccExpWtMeanVar
Incremental calculation of weighted mean and variance
by Tony Finch
=#

mutable struct AccumExpWtMean{T, F} <: Accumulator{T}
    n::Int
    alpha::T
    mean::T
    fn::F
    AccumExpWtMean(::Type{T}=Float64, fn::F=identity; alpha::T=0.5) where {T, F} = new{T, F}(0, alpha, zero(T), fn)
end

(acc::AccumExpWtMean{T, F})() where {T, F} = acc.mean
(acc::AccumExpWtMean{T, F})(x) where {T, F} = (acc.n += 1; acc.mean += acc.alpha * (x - acc.mean); acc)
(acc::AccumExpWtMean{T,FN})(xs::Seq) where {T,FN} = (acc.n += length(xs); xxs = map(acc, xs); acc)

function (acc::AccumExpWtMean{T, F})(xs::Seq) where {T, F}
    @inbounds for i ∈ eachindex(xs)
        acc(xs[i])
    end
    acc
end

mutable struct AccumExpWtMeanVar{T, F} <: Accumulator{T}
    n::Int
    alpha::T
    mean::T
    svar::T
    fn::F
    AccumExpWtMeanVar(::Type{T}=Float64, fn::F=identity; alpha::T=0.5) where {T, F} = new{T, F}(0, alpha, zero(T), zero(T), fn)
end

function(acc::AccumExpWtMeanVar{T, F})() where {T, F}
    unbiased_var = acc.svar / (acc.n = 1)
    (acc.mean, unbiased_var)
end

function (acc::AccumExpWtMeanVar{T, F})(x) where {T, F}
    xx = fn(T(x))
    acc.n += 1
    diff = xx - acc.mean
    incr = acc.alpha * diff
    acc.mean += acc.alpha * (xx - acc.mean)
    acc.svar = (one(T) - acc.alpha) * (acc.svar + diff * incr)
    acc
end

function (acc::AccumExpWtMeanVar{T, F})(xs::Seq) where {T, F}
    @inbounds for i ∈ eachindex(xxs)
        acc(xxs[i])
    end
    acc
end

acc_min(acc::AccumMin{T, F}) where {T, F} = acc.min
acc_max(acc::AccumMax{T, F}) where {T, F} = acc.max

acc_min(acc::AccumExtrema{T, F}) where {T, F} = acc.min
acc_max(acc::AccumExtrema{T, F}) where {T, F} = acc.max
acc_nmin(acc::AccumExtrema{T, F}) where {T, F} = acc.nmin
acc_nmax(acc::AccumExtrema{T, F}) where {T, F} = acc.nmax
acc_midrange(acc::AccumExtrema{T, F}) where {T, F} = (acc.max / 2) + (acc.min / 2)
                                                                      
acc_mean(acc::AccumStats{T, F}) where {T, F} = T(acc.m1)
acc_var(acc::AccumStats{T, F}) where {T, F} = T(acc.m2 / (acc.n - 1))
acc_std(acc::AccumStats{T, F}) where {T, F} = T(sqrt(var(x)))
acc_skew(acc::AccumStats{T, F}) where {T, F} = T(sqrt(acc.n) * acc.m3 / (acc.m2 * sqrt(acc.m2)))
acc_kurt(acc::AccumStats{T, F}) where {T, F} = T((acc.n * acc.m4) / (acc.m2^2) - 3)
                                                                      
acc_mean(acc::AccumMean{T, F}) where {T, F} = acc.mean
acc_mean(acc::AccumGeometricMean{T, F}) where {T, F} = acc()
acc_mean(acc::AccumHarmonicMean{T, F}) where {T, F} = acc()

acc_mean(acc::AccumMeanVar{T, F}) where {T, F} = acc.mean
acc_var(acc::AccumMeanVar{T, F}) where {T, F} = acc.svar / (acc.n - 1)
acc_std(acc::AccumMeanVar{T, F}) where {T, F} = sqrt(acc.svar / (acc.n - 1))

acc_mean(acc::AccumExpWtMean{T, F}) where {T, F} = acc.mean

acc_mean(acc::AccumExpWtMeanVar{T, F}) where {T, F} = acc.mean
acc_var(acc::AccumExpWtMeanVar{T, F}) where {T, F} = acc.svar / (acc.n - 1)
acc_std(acc::AccumExpWtMeanVar{T, F}) where {T, F} = sqrt(acc_var(acc))
