module Accumulators

export Accumulator,
    # Acc is a typed accumlator
    AccCount,
    AccMinimum, AccMaximum, AccExtrema,
    AccSum, AccProd,
    AccMean, AccGeoMean, AccHarmMean, AccGenMean,
    AccMeanAndVar, AccMeanAndStd, AccStats,
    AccExpWtMean, AccExpWtMeanVar, AccExpWtMeanStd

using StatsBase
using LoopVectorization: @turbo, @tturbo
using VectorizedStatistics
using Chain

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

abstract type Accumulator{T} <: Function end

const Seq = Union{AbstractVector{T}, Tuple{Vararg{T}}} where {T}

seq(x::AbstractVector{T}) where {T} = x
seq(x::NTuple{N,T}) where {N,T} = x

logabs(x) = log(abs(x))
sumlogabs(xs::A) where {T, A<:AbstractVector{T}} = vsum(map(logabs, xs))
sumlogabs(xs::NTuple{N,T}) where {T, N} = sum(map(logabs, xs))
    
include("compose.jl")
# include("compose2.jl")
# include("augment.jl")
# include("absvals.jl")

end  # Accumulators
