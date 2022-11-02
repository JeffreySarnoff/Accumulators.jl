@testset "simple accumulators" begin
   acc = AccCount()
   for x in one2five
       acc(x)
   end
   @test acc() == 5
  
   acc = AccMin()
   for x in sortn8
       acc(x)
   end
   @test acc() == sortn8[1]
  
   acc = AccMax()
   for x in sortn8
       acc(x)
   end
   @test acc() == sortn8[end]

   acc = AccExtrema()
   for x in sortn8
       acc(x)
   end
   @test acc() == (sortn8[1], sortn8[end])
end

