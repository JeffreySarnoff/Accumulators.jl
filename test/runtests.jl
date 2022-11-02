using Accumulators, StatsBase, Test

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
