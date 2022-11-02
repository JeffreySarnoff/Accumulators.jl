using Accumulators, StatsBase, Test

≐(a, b) = isapprox(a, b; rtol=eps(a)^0.875)

@inline function accvals(acc, vals)
    for x in vals
        acc(x)
    end
end

include("testvals.jl")

@testset "compose" begin
    include("compose.jl")
end

@testset "absvals" begin
    include("absvals.jl")
end

@testset "augment" begin
    include("augment.jl")
end
