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


abstract type Accumulator{T} end

include("compose.jl")
include("augment.jl")
include("absvals.jl")


end  # Accumulators
