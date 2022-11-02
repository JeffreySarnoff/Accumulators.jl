module Accumulators

export Accumulator,
    AccumCount,
    AccumMin, AccumMax, 
    AccumExtrema, acc_midrange,
    AccumSum, AccumProd,
    AccumMean, AccumGeometricMean, AccumHarmonicMean,
    AccumMeanVar,
    AccMinAbs, AccMaxAbs, AccMeanAbs, AccMeanAbs2,
    AccSumAbs, AccProdAbs, 
    AccMeanVarAbs, AccMeanVarAbs2, 
    AccStats, acc_count, acc_mean, acc_var, acc_std, acc_skew, acc_kurt,
    AccExpWtMean, AccExpWtMeanVar

abstract type Accumulator{T} end

include("compose.jl")
include("augment.jl")
include("absvals.jl")


end  # Accumulators
