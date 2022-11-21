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
    AccumMin, AccumMax, AccumExtrema,
    AccumSum, AccumProd,
    AccumMean, AccumGeometricMean, AccumHarmonicMean,
    AccumMeanVar,
    AccumStats,
    AccumExpWtMean, AccumExpWtMeanVar,
    # acc_ are convienice getters for appropriate accumulators
    acc_count, acc_min, acc_max, acc_nmin, acc_nmax, acc_midrange,
    acc_mean, acc_var, acc_std, acc_skew, acc_kurt

using LoopVectorization: @turbo, @tturbo
using VectorizedStatistics

#

accumulator_type = Float64
if isdefined(Main, :AccumulatorNumType)
    accumulator_type = AccumulatorNumType
elseif haskey(ENV, "AccumulatorNumType")
    defaultnumtype = parse(ENV["AccumulatorNumType"])
    if (isa(defaultnumtype, Type) && isa(defaultnumtype, Number))
        accumulator_type = defaultnumtype
    end
end
const AccNum = accumulator_type    

#

abstract type Accumulator{T} <:Function end

logabs(x) = log(abs(x))
sumlogabs(xs::A) where {T, A<:AbstractVector{T}} = vsum(map(logabs, xs))
sumlogabs(xs::NTuple{N,T}) where {T, N} = sum(map(logabs, xs))
    
include("compose.jl")
include("augment.jl")
include("absvals.jl")


end  # Accumulators
