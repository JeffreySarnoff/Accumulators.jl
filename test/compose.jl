#=
     AccCount, 
     AccMinimum, AccMaximum, AccExtrema,

     AccSum, AccProd,

     AccMean, AccGeoMean, AccHarmMean, AccGenMean,

     AccMeanVar, AccMeanStd, AccStats,
     AccExpWtMean, AccExpWtMeanVar, AccExpWtMeanStd
=#

@testset "simple accumulators" begin
   acc = AccCount()
   accvals(acc, one2five)
   @test acc() == 5
  
   acc = AccMinimum()
   accvals(acc, sortn8)
   @test acc() == sortn8[1]
  
   acc = AccMaximum()
   accvals(acc, sortn8)
   @test acc() == sortn8[end]

   acc = AccExtrema()
   accvals(acc, sortn8)
   @test acc() == (sortn8[1], sortn8[end])
end

@testset "sum, prod" begin
    acc = AccSum()
    accvals(acc, sortn8)
    @test acc() == foldl(+, sortn8; init=0.0)

    acc = AccProd()
    accvals(acc, sortn8)
    @test acc() == foldl(*, sortn8; init=1.0)
end

@testset "means" begin
    acc = AccMean()
    accvals(acc, sortn8)
    @test acc() ≐ mean(sortn8)
   
    acc = AccGeoMean()
    accvals(acc, sortn8)
    @test acc() ≐ geomean(map(abs, sortn8))

    acc = AccHarmMean()
    accvals(acc, sortn8)
    @test acc() ≐ harmmean(sortn8)

    acc = AccGenMean(power=2)
    accvals(acc, sortn8)
    @test acc() ≐ genmean(sortn8, 2)
end
   
@testset "stats" begin
    acc = AccMeanAndVar()
    accvals(acc, sortn8)
    @test all(map(≐, Tuple(acc()), mean_and_var(sortn8))

    acc = AccMeanAndStd()
    accvals(acc, sortn8)
    @test all(map(≐, Tuple(acc()), mean_and_std(sortn8))

    acc = AccStats()
    accvals(acc, sortn8)
    tst = (nobs = length(sortn8), mean = mean(sortn8), var = var(sortn8), std = std(sortn8), skewness = skewness(sortn8), kurtosis = kurtosis(sortn8)) 
    @test all(map(≐, Tuple(acc()), Tuple(tst))
end

@testset "expwt" begin
#     AccExpWtMean, AccExpWtMeanVar, AccExpWtMeanStd   
end

