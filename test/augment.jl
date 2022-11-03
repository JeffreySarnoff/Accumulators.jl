@testset "simple accumulators" begin
   acc = AccumCount(fn=abs2)
   accvals(acc, one2five)
   @test acc() == 5*5
  
   acc = AccumMin(fn=abs2)
   accvals(acc, sortn8)
   @test acc() == minimum(map(abs2, sortn8))
  
   acc = AccumMax(fn=abs2)
   accvals(acc, sortn8)
   @test acc() == maximum(map(abs2, sortn8))

   acc = AccumExtrema(fn=abs2)
   accvals(acc, sortn8)
   @test acc() == (minimum(map(abs2, sortn8)), maximum(map(abs2, sortn8)))
end

@testset "sum, prod" begin
    acc = AccumSum(fn=abs2)
    accvals(acc, sortn8)
    @test acc() == foldl(+, map(abs2,sortn8); init=0.0)

    acc = AccumProd(fn=abs2)
    accvals(acc, sortn8)
    @test acc() == foldl(*, map(abs2,sortn8); init=1.0)
end

@testset "means" begin
    acc = AccumMean(fn=abs2)
    accvals(acc, sortn8)
    @test acc() ≐ mean(map(abs2,sortn8))
   
    acc = AccumGeometricMean(fn=abs2)
    accvals(acc, sortn8)
    @test acc() ≐ geomean(map(abs2, sortn8))

    acc = AccumHarmonicMean(fn=abs2)
    accvals(acc, sortn8)
    @test acc() ≐ harmmean(map(abs2,sortn8))
end
   
@testset "stats" begin
   
end

@testset "expwt" begin
end

#=
 AccMeanVar,
    AccStats,
    AccExpWtMean, AccExpWtMeanVar,
    acc_count, acc_mean, acc_var, acc_std, acc_skew, acc_kurt,
    acc_midrange,
=#
