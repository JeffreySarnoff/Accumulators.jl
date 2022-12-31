module Accumulators

export Seq, seq,
    Accumulator,  
    # Acc_ is a typed accumlator
    AccCount,
    AccMinimum, AccMaximum, AccExtrema,
    AccSum, AccProd,
    AccMean, AccGeoMean, AccHarmMean, AccGenMean,
    AccMeanAndVar, AccMeanAndStd, AccStats,
    AccExpWtMean, AccExpWtMeanVar, AccExpWtMeanStd,
    midrange, proportionalrange, 
    nminima, nmaxima

abstract type Accumulator{T,F} <: Function end

const Seq = Union{AbstractArray{T,N}, NTuple{N,T}}} where {N,T}

seq(x::AbstractArray{T,N}) where {N,T} = x
seq(x::NTuple{N,T}) where {N,T} = x

include("support.jl")

include("accum_one.jl")

include("accum_two.jl")
include("accum_three.jl")
include("accum_many.jl")


end  # Accumulators
