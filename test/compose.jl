@inline function accvals(acc, vals)
    for x in vals
      acc(x)
    end
end

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
    accvals(acc, randn8)
    @test acc() == mean(randn8)
   
end
   
   
