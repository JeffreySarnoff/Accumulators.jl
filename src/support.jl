const AccNum = Float64

# shift new obs into tuple, drop oldest obs
#=

    # leftshift

    # newest(cache) = last(cache)
    # oldest(cache) = first(cache)
    #
    # newest obs at cache[end], oldest obs at cache[1]
    # oldest obs at cache[1], newest obs at cache[end]

    init_cache   = (x₀, x₀, x₀, x₀)
    new_obs = x₁
    update_cache = (x₀, x₀, x₀, x₁)
    new_obs = x₂
    update_cache = (x₀, x₀, x₁, x₂)
    new_obs = x₃
    update_cache = (x₀, x₁, x₂, x₃)
    new_obs = x₄
    update_cache = (x₁, x₂, x₃, x₄)
    new_obs = x₅
    update_cache = (x₂, x₃, x₄, x₅)
    new_obs = x₆
    update_cache = (x₃, x₄, x₅, x₆)
    new_obs = x₇
    update_cache = (x₄, x₅, x₆, x₇)
    new_obs = x₈
    update_cache = (x₅, x₆, x₇, x₈)


    # rightshift

    # newest(cache) = first(cache)
    # oldest(cache) = last(cache)
    #
    # newest obs at cache[1], oldest obs at cache[end]
    # oldest obs at cache[end], newest obs at cache[1]

    init_cache   = (x₀, x₀, x₀, x₀)
    new_obs = x₁
    update_cache = (x₁, x₀, x₀, x₀)
    new_obs = x₂
    update_cache = (x₂, x₁, x₀, x₀)
    new_obs = x₃
    update_cache = (x₃, x₂, x₁, x₀)
    new_obs = x₄
    update_cache = (x₄, x₃, x₂, x₁)
    new_obs = x₅
    update_cache = (x₅, x₄, x₃, x₂)
    new_obs = x₆
    update_cache = (x₆, x₅, x₄, x₃)
    new_obs = x₇
    update_cache = (x₇, x₆, x₅, x₄)
    new_obs = x₈
    update_cache = (x₈, x₇, x₆, x₅)


=#

@inline cache_eltype(accum::AccumCaches{T,F,N}) where {T,F,N} = T
@inline cache_length(accum::AccumCaches{T,F,N}) where {T,F,N} = N
@inline cache_size(accum::AccumCaches{T,F,N}) where {T,F,N} = sizeof(NTuple{N,T})
@inline cache_bits(accum::AccumCaches{T,F,N}) where {T,F,N} = cache_size(accum) * 8

#=
adaptation of code borrowed from 
https://github.com/BioJulia/Kmers.jl/blob/master/src/tuple_bitflipping.jl


rightshift & leftshift


These methods are micro-optimised for use with an NTuple{N,T}.
- shift out the oldest observed value, drops the [first, last]
- shift or rewrite the youngest..penultimate_oldest values
  -- overwriting or rewriting each value with the next younger
  -- this leaves the youngest twice given
  -- overwrite or replace the most extremal of the two youngest
      with newest observed value

=#

@inline function dropfirst(x::NTuple{N,T}) where {N,T}
    ntuple(i->x[i+1], N-1)
end

@inline function droplast(x::NTuple{N,T}) where {N,T}
    ntuple(i->x[i], N-1)
end

@inline function popfirst(x::NTuple{N,T}) where {N,T}
    x[2:end]
end

@inline function poplast(x::NTuple{N,T}) where {N,T}
    x[1:end-1]
end

@inline function newfirst(x::NTuple{N,T}, newvalue::T) where {N,T}
    ntuple(i->ifelse(isone(i), newvalue, x[i]), N)
end

@inline function newlast(x::NTuple{N,T}, newvalue::T) where {N,T}
    ntuple(i->ifelse((i==N), newvalue, x[i]), N)
end

@inline function pushfirst(x::NTuple{N,T}, newvalue::T) where {N,T}
    (newvalue, x...)
end

@inline function pushlast(x::NTuple{N,T}, newvalue::T) where {N,T}
    (x..., newvalue)
end



@inline function rightshift1(x::NTuple{N,T}, newvalue::T) where {N,T}
    return _rightshift1(newvalue, @view(x,2,N))
end

@inline function _rightshift1(hewvalue::T, xs::Vararg{T,N}) where {N,T}

@inline function _rightshift1(xs::Vararg{T,N}) where {N,T}
    return (xs[(head >> nbits) | carry, _rightshift_carry(nbits, (head & ((one(UInt64) << nbits) - 1)) << (64 - nbits), tail...)...)
end

@inline _rightshift1(carry::T) where {T} = ()


#=
These methods are micro-optimised (or should be!!!) 
for shifting the bits in an NTuple of unsigned integers, 
carrying the bits "shifted off" one word over to the next word.
The carry can also be "seeded" so as other methods like
pushfirst and pushlast can be efficiently implemented without
duplication of code or less efficient implementations
that first shift and then insert an element.


rightshift_carry & leftshift_carry

These methods are micro-optimised (or should be!!!) for shifting the bits in 
an NTuple of unsigned integers, carrying the bits "shifted off" one word 
over to the next word. The carry can also be "seeded" so as other methods like
pushfirst and pushlast can be efficiently implemented without duplication of code
or less efficient implementations that first shift and then insert an element.
=#

@inline function rightshift_carry(x::NTuple{N,UInt64}, nbits::Integer, prevcarry=zero(UInt64)) where {N}
    return _rightshift_carry(nbits, prevcarry, x...)
end

@inline function _rightshift_carry(nbits::Integer, carry::UInt64, head::UInt64, tail...)
    return ((head >> nbits) | carry, _rightshift_carry(nbits, (head & ((one(UInt64) << nbits) - 1)) << (64 - nbits), tail...)...)
end

@inline _rightshift_carry(nbits::Integer, carry::UInt64) = ()

@inline function leftshift_carry(x::NTuple{N,UInt64}, nbits::Integer, prevcarry::UInt64=zero(UInt64)) where {N}
    _, newbits = _leftshift_carry(nbits, prevcarry, x...)
    return newbits
end

@inline function _leftshift_carry(nbits::Integer, prevcarry::UInt64, head::UInt64, tail...)
    carry, newtail = _leftshift_carry(nbits, prevcarry, tail...)
    return head >> (64 - nbits), ((head << nbits) | carry, newtail...)
end

@inline _leftshift_carry(nbits::Integer, prevcarry::UInt64) = prevcarry, ()


# math functions for internal use

logabs(x) = log(abs(x))
sumlogabs(xs::A) where {T,A<:AbstractVector{T}} = vsum(map(logabs, xs))
sumlogabs(xs::NTuple{N,T}) where {T,N} = sum(map(logabs, xs))

