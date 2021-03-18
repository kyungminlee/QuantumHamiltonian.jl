using Test
using QuantumHamiltonian

@testset "FrozenSortedArray" begin
  @test_throws ArgumentError FrozenSortedArrayIndex(['b', 'a', 'c'])
  @test_throws ArgumentError FrozenSortedArrayIndex(['a', 'a', 'b'])

  arr = FrozenSortedArrayIndex(['a', 'b', 'c', 'd', 'e'])
  @test eltype(arr) == Pair{Char, Int}
  @test length(arr) == 5
  @test keys(arr) == ['a', 'b', 'c', 'd', 'e']
  @test collect(arr) == ['a' => 1, 'b' => 2, 'c' => 3, 'd' => 4, 'e' => 5]
end
