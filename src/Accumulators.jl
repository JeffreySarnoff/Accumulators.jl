module Accumulators

export Seq, seq,
    Accumulator, AccumCacher,
    # Acc_ is a typed accumlator
    AccCount,
    AccMinimum, AccMaximum, AccExtrema,
    AccSum, AccProd,
    AccMean, AccGeoMean, AccHarmMean, AccGenMean,
    AccMeanAndVar, AccMeanAndStd, AccStats,
    AccExpWtMean, AccExpWtMeanVar, AccExpWtMeanStd,
    midrange, proportionalrange,
    nminima, nmaxima

using StatsBase, VectorizedStatistics

#= ---- abstract types ---- =#

"""
    Accumulator{T,F}
"""
Accumulator

abstract type Accumulator{T,F} <: Function end

"""
    AccumCacher{T,F}
"""
AccumCacher

abstract type AccumCacher{T,F} <: Accumulator{T,F} end

"""
    AccumCaches{T,F,N}
"""
AccumCaches

abstract type AccumCaches{T,F,N} <: Accumulator{T,F} end

#= ---- generalized sequence ---- =#

"""
    Seq
"""
Seq

const Seq = Union{AbstractArray{T,N},NTuple{N,T}} where {N,T}

"""
    seq
"""
seq

seq(x::AbstractArray{T,N}) where {N,T} = x
seq(x::NTuple{N,T}) where {N,T} = x

const AccNum = Float64
const AccumNum = Float64

#= ---- code organization ---- =#

include("support.jl")

include("accum_one.jl")

include("accum_two.jl")
include("accum_three.jl")
include("accum_many.jl")

#= ---- notes ---- =#

end  # Accumulators
