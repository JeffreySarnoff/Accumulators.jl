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
    # acc_ are convienice getters for appropriate accumulators
    acc_count, acc_mean, acc_var, acc_std, acc_skew, acc_kurt,
    acc_midrange,
    # Accum is Acc with a specified function applied to each new 'x'
    AccumCount,
    AccumMin, AccumMax, AccumExtrema,
    AccumSum, AccumProd,
    AccumMean, AccumGeometricMean, AccumHarmonicMean,
    AccumMeanVar,
    AccumStats,
    AccumExpWtMean, AccumExpWtMeanVar

using LoopVectorization: @turbo, @tturbo
using VectorizedStatistics

const Seq = Union{AbstractVector{T}, NTuple{N,T}} where {N,T}

logabs(x) = log(abs(x))

abstract type Accumulator{T} end

include("compose.jl")
include("augment.jl")
include("absvals.jl")


end  # Accumulators
