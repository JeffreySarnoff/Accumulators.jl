using Accumulators, StatsBase, Test

≐(a, b) = isapprox(a, b; rtol=eps(a)^0.75)

@inline function accvals(acc, vals)
    for x in vals
        acc(x)
    end
end

include("constants.jl")

@testset "compose" begin
    include("compose.jl")
end

@testset "augment" begin
    include("augment.jl")
end
