using Accumulators, StatsBase, Test

include("testvals.jl")

@inline function accvals(acc, vals)
    for x in vals
        acc(x)
    end
end

@testset "compose" begin
    include("compose.jl")
end

@testset "absvals" begin
    include("absvals.jl")
end

@testset "augment" begin
    include("augment.jl")
end
