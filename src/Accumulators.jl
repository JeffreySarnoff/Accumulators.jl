module Accumulators

export Accumulator,
    # Acc is a typed accumlator
    AccCount,
    AccMin, AccMax, AccExtrema,
    AccSum, AccProd,
    AccMean, AccGeometricMean, AccHarmonicMean,
    AccMeanVar,
    AccStats,
    AccExpWtMean, AccExpWtMeanVar,
    AccMinAbs, AccMinAbs2, AccMaxAbs, AccMaxAbs2,
    AccSumAbs, AccSumAbs2, AccProdAbs, AccProdAbs2, 
    AccMeanAbs, AccMeanAbs2, AccMeanVarAbs, AccMeanVarAbs2,
    # Accum is Acc with a prespecified function applied to each new 'x'
    AccumCount,
    AccumMin, AccumMax, AccumExtrema,
    AccumSum, AccumProd,
    AccumMean, AccumGeometricMean, AccumHarmonicMean,
    AccumMeanVar,
    AccumStats,
    AccumExpWtMean, AccumExpWtMeanVar
    # acc_ are convienice getters for appropriate accumulators
    acc_count, acc_min, acc_max, acc_nmin, acc_nmax, acc_midrange,
    acc_mean, acc_var, acc_std, acc_skew, acc_kurt

using LoopVectorization: @turbo, @tturbo
using VectorizedStatistics

const Seq = Union{AbstractVector{T}, NTuple{N,T}} where {N,T}

logabs(x) = log(abs(x))
sumlogabs(xs::Seq) = vsum(map(logabs, xs))

abstract type Accumulator{T} end

include("compose.jl")
include("augment.jl")
include("absvals.jl")


end  # Accumulators
