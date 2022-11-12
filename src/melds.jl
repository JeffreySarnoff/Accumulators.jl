abstract type Melds <: Function end

function meld(meldwith::Function, values::Vararg{N}) where {N}
    map(meldwith, values)
end

function meld(meldwith::Function, values::T) where {T<:Tuple}
  map(meldwith, values)
end

function meld(meldwith::Function, values::T) where {T<:AbstractVector}
  map(meldwith, values)
end

function meld(meldwith::Function, values::T) where {T<:AbstractArray}
  map(meldwith, values)
end

  

#=
julia> xs'
1×6 adjoint(::Vector{Int64}) with eltype Int64:
 1  2  3  4  5  6

julia> ys'
1×6 adjoint(::Vector{Int64}) with eltype Int64:
 11  12  13  14  15  16

julia> map(mean,mapreduce(a->map(abs,a), zip, [xs, ys]))
6-element Vector{Float64}:
  5.5
  6.5
  7.5
  8.5
  9.5
 10.5

=#

#=


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

=#
