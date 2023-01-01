const AccNum = F

# shift new obs into tuple, drop oldest obs
#=
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

=#



# math functions for internal use

logabs(x) = log(abs(x))
sumlogabs(xs::A) where {T, A<:AbstractVector{T}} = vsum(map(logabs, xs))
sumlogabs(xs::NTuple{N,T}) where {T, N} = sum(map(logabs, xs))

