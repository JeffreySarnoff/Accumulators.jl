@testset "simple accumulators" begin
   acc = AccCount()
   accvals(acc, one2five)
   @test acc() == 5
  
   acc = AccMin()
   accvals(acc, sortn8)
   @test acc() == sortn8[1]
  
   acc = AccMax()
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
    @test acc() == mean(sortn8)
   
    acc = AccGeometricMean()
    accvals(acc, sortn8)
    @test acc() == geomean(map(abs, sortn8))

    acc = AccHarmonicMean()
    accvals(acc, sortn8)
    @test acc() == harmmean(sortn8)
end
   
   
