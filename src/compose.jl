#=
     AccCount, 
     AccMin, AccMax, AccExtrema, 
     AccSum, AccProd,
     AccMean, AccGeometricMean, AccHarmonicMean,
     AccMeanVar,
     AccExpWtMean, AccExpWtMeanVar
=#

mutable struct AccCount{T} <: Accumulator{T}
    n::T
    AccCount(::Type{T}=Int64) where {T} = new{T}(zero(T))
end

(acc::AccCount{T})() where {T} = (acc.n)
(acc::AccCount{T})(x) where {T} = (acc.n += one(T); acc)
(acc::AccCount{T})(xs::Seq) where {T} = (acc.n += T(length(xs)); acc)

mutable struct AccMin{T} <: Accumulator{T}
    min::T
    AccMin(::Type{T}=Float64) where {T} =
        (T <: Integer) ? new{T}(typemax(T)) : new{T}(floatmax(T))
end

(acc::AccMin{T})() where {T} = (acc.min)
(acc::AccMin{T})(x) where {T} = (acc.min = ifelse(x < acc.min, T(x), acc.min); acc)
(acc::AccMin{T})(xs::Seq) where {T} = (x = T(vminimum(xs)); acc.min = ifelse(x < acc.min, x, acc.min); acc)

mutable struct AccMax{T} <: Accumulator{T}
    max::T
    AccMax(::Type{T}=Float64) where {T} =
        (T <: Integer) ? new{T}(typemin(T)) : new{T}(floatmin(T))
end

(acc::AccMax{T})() where {T} = (acc.max)
(acc::AccMax{T})(x) where {T} = (acc.max = ifelse(acc.max < x, T(x), acc.max); acc)
(acc::AccMax{T})(xs::Seq) where {T} = (x = T(vmaximum(xs)); acc.max = ifelse(x > acc.max, x, acc.max); acc)

mutable struct AccExtrema{T} <: Accumulator{T}
    min::T
    max::T
    AccExtrema(::Type{T}=Float64) where {T} =
        (T <: Integer) ? new{T}(typemax(T), typemin(T)) : new{T}(floatmax(T), floatmin(T))
end

(acc::AccExtrema{T})() where {T} = (acc.min, acc.max)
(acc::AccExtrema{T})(x) where {T} = 
    (acc.min = ifelse(x < acc.min, T(x), acc.min);
     acc.max = ifelse(acc.max < x, T(x), acc.max); acc)

function (acc::AccExtrema{T})(xs::Seq) where {T}
   mn, mx = map(T, vminimum(xs))
   acc.min = ifelse(mn < acc.min, mn, acc.min)
   acc.max = ifelse(mx > acc.max, mx, acc.max)
   acc
end

acc_max(acc::AccExtrema{T}) where {T} = acc.max
acc_min(acc::AccExtrema{T}) where {T} = acc.min
acc_midrange(acc::AccExtrema{T}) where {T} = (acc.max / 2) + (acc.min / 2)

mutable struct AccSum{T} <: Accumulator{T}
    sum::T
    AccSum(::Type{T}=Float64) where {T} = new{T}(zero(T),)
end

(acc::AccSum{T})() where {T} = (acc.sum)
(acc::AccSum{T})(x) where {T} = (acc.sum += x; acc)
(acc::AccSum{T})(xs::Seq) where {T} = (x = T(vsum(xs)); acc.sum += x; acc)

mutable struct AccProd{T} <: Accumulator{T}
    prod::T
    AccProd(::Type{T}=Float64) where {T} = new{T}(one(T),)
end

(acc::AccProd{T})() where {T} = (acc.prod)
(acc::AccProd{T})(x) where {T} = (acc.prod *= x; acc)

function (acc::AccProd{T})(xs::Seq) where {T}
    Σ = one(T)
    @turbo for i ∈ eachindex(xs)
        Σ *= xs[i]
    end
    acc.prod *= Σ
    acc
end
             
mutable struct AccMean{T} <: Accumulator{T}
    n::Int
    mean::T
    AccMean(::Type{T}=Float64) where {T} = new{T}(0, zero(T))
end

(acc::AccMean{T})() where {T} = (acc.mean)
(acc::AccMean{T})(x) where {T} =
    (acc.n += 1; acc.mean += (x - acc.mean) / acc.n; acc)
(acc::AccMean{T})(xs::Seq) where {T} =
     (acc.n += length(xs); x = T(vmean(xs)); acc.mean += (x - acc.mean) / acc.n; acc)

# geometric mean (of abs(xs))
# see https://github.com/stdlib-js/stats/blob/main/incr/gmean/lib/main.js
mutable struct AccGeometricMean{T} <: Accumulator{T}
    n::Int
    sumlog::T
    AccGeometricMean(::Type{T}=Float64) where {T} = new{T}(0, zero(T))
end

(acc::AccGeometricMean{T})() where {T} = (iszero(acc.n) ? one(T) : exp(acc.sumlog / acc.n))
(acc::AccGeometricMean{T})(x) where {T} = (acc.n += 1; acc.sumlog += logabs(x); acc)
(acc::AccGeometricMean{T})(xs::Seq) where {T} =
     (acc.n += length(xs); s = sumlogabs(xs); acc.sumlog += s; acc)

# harmonic mean
# see https://github.com/stdlib-js/stats/blob/main/incr/hmean/lib/main.js
mutable struct AccHarmonicMean{T} <: Accumulator{T}
    n::Int
    hmean::T
    AccHarmonicMean(::Type{T}=Float64) where {T} = new{T}(0, zero(T))
end

(acc::AccHarmonicMean{T})() where {T} = (iszero(acc.n) ? one(T) : acc.n / acc.hmean)
(acc::AccHarmonicMean{T})(x) where {T} = (acc.n += 1; acc.hmean += (one(T) / x); acc)
(acc::AccHarmonicMean{T})(xs::Seq) where {T} =
     (acc.n += length(xs); s = vsum(map(x->one(T)/x, xs)); acc.hmean += s; acc)

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
    acc
end

function (acc::AccMeanVar{T})(xs::Seq) where {T}
    n = length(xs)
    m = vmean(xs)
    # v = vvar(xs)
    # s = v * (n - 1)
    acc.n += n
    if !iszero(acc.n)
        oldmean = acc.mean
        acc.mean = oldmean + (m - oldmean) / acc.n
        acc.svar = acc.svar + (m - oldmean) * (m - acc.mean)
    else
        acc.mean = m
    end                                              
    acc
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
    accn = acc.n
    accm1 = acc.m1
    accm2 = acc.m2
    accm3 = acc.m3
    accm4 = acc.m4
                                                  
    n1 = accn
    accn += 1
    delta = x - accm1
    delta_n = delta / accn
    delta_n2 = delta_n * delta_n
    term1 = delta * delta_n * n1
    accm1 += delta_n
    accm4 += term1 * delta_n2 * (accn * accn - 3*accn + 3) + 
              6 * delta_n2 * accm2 - 
              4 * delta_n * accm3
    accm3 += term1 * delta_n * (n - 2) - 3 * delta_n * accm2
    accm2 += term1

    acc.n = accn
    acc.m1 = accm1
    acc.m2 = accm2
    acc.m3 = accm3
    acc.m4 = accm4
    acc
end

function (acc::AccStats{T})(xs::Seq) where {T}
    Σ = one(T)
    @inbounds for i ∈ eachindex(xs)
        acc(xs[i])
    end
    acc
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
(acc::AccExpWtMean{T})(x) where {T} = (acc.n += 1; acc.mean += acc.alpha * (x - acc.mean); acc)

function (acc::AccExpWtMean{T})(xs::Seq) where {T}
    @inbounds for i ∈ eachindex(xs)
        acc(xs[i])
    end
    acc
end
                                                  
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
    acc
end

function (acc::AccExpWtMeanVar{T})(xs::Seq) where {T}
    @inbounds for i ∈ eachindex(xs)
        acc(xs[i])
    end
    acc
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

