const AccNum = Float64

# math functions for internal use

logabs(x) = log(abs(x))
sumlogabs(xs::A) where {T, A<:AbstractVector{T}} = vsum(map(logabs, xs))
sumlogabs(xs::NTuple{N,T}) where {T, N} = sum(map(logabs, xs))

