using Accumulators, StatsBase, RandomNumbers, Test

rng = Xoshiro(1618);
ns = (4,16,63,255);
urands = rand.(rng, ns);
nrands = randn.(rng, ns);

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

#testset "augment" begin
#   include("augment.jl")
#nd
