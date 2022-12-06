module Accumulators

export Accumulator,
    # Acc is a typed accumlator
    AccCount,
    AccMinimum, AccMaximum, AccExtrema,
    AccSum, AccProd,
    AccMean, AccGeoMean, AccHarmMean,
    AccMeanVar, AccMeanStd, AccStats,
    AccExpWtMean, AccExpWtMeanVar, AccExpWtMeanStd,
#=
    AccMinimumAbs, AccMaximumAbs, AccExtremaAbs,
    AccSumAbs, AccProdAbs,

    AccMeanAbs, AccGeoMeanAbs, AccHarmMeanAbs,
    AccMeanVarMag, AccMeanStdMag, AccStatsMag,
    
    AccMagMinimum, AccMagMaximum, AccMagExtrema,
    AccMagSum, AccMagProd,
    AccAbsMinimum, AccMaximum, AccMagExtrema,
    AccMagSum, AccMagProd,

    AccMeanMag, AccGeoMeanMag, AccHarmMeanMag,
    AccMeanVarMag, AccMeanStdMag, AccStatsMag,

    AccExpWtMean, AccExpWtMeanVar,
    AccMinAbs, AccMinAbs2, AccMaxAbs, AccMaxAbs2,
    AccSumAbs, AccSumAbs2, AccProdAbs, AccProdAbs2, 
    AccMeanAbs, AccMeanAbs2, AccMeanVarAbs, AccMeanVarAbs2,
    # for paired data streams
    AccCov, AccCor,
    # Accum is Acc with a prespecified function applied to each new 'x'
    AccumMin, AccumMax, AccumExtrema,
    AccumSum, AccumProd,
    AccumMean, AccumGeometricMean, AccumHarmonicMean,
    AccumMeanVar,
    AccumStats,
    AccumExpWtMean, AccumExpWtMeanVar,
=#

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

const Seq = Union{AbstractVector{T}, Tuple{Vararg{T}}} where {T}

seq(x::AbstractVector{T}) where {T} = x
seq(x::NTuple{N,T}) where {N,T} = x

logabs(x) = log(abs(x))
sumlogabs(xs::A) where {T, A<:AbstractVector{T}} = vsum(map(logabs, xs))
sumlogabs(xs::NTuple{N,T}) where {T, N} = sum(map(logabs, xs))
    
include("compose.jl")
include("compose2.jl")
include("augment.jl")
include("absvals.jl")


end  # Accumulators
