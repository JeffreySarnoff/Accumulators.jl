@testset "simple accumulators" begin
   acc = AccumCount()
   accvals(acc, one2five)
   @test acc() == 5
  
   acc = AccumMin()
   accvals(acc, sortn8)
   @test acc() == sortn8[1]
  
   acc = AccumMax()
   accvals(acc, sortn8)
   @test acc() == sortn8[end]

   acc = AccumExtrema()
   accvals(acc, sortn8)
   @test acc() == (sortn8[1], sortn8[end])
end

@testset "sum, prod" begin
    acc = AccumSum()
    accvals(acc, sortn8)
    @test acc() == foldl(+, sortn8; init=0.0)

    acc = AccumProd()
    accvals(acc, sortn8)
    @test acc() == foldl(*, sortn8; init=1.0)
end

@testset "means" begin
    acc = AccumMean()
    accvals(acc, sortn8)
    @test acc() ≐ mean(sortn8)
   
    acc = AccumGeometricMean()
    accvals(acc, sortn8)
    @test acc() ≐ geomean(map(abs, sortn8))

    acc = AccumHarmonicMean()
    accvals(acc, sortn8)
    @test acc() ≐ harmmean(sortn8)
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
