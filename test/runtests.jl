using Random
using StringAlgorithms
using Test

@testset "RadixSort" begin
    b = zeros(Int32, 3)
    a = [0,1,2]

    Random.seed!( 0 )
    
    prefix_result = randstring(40961)

    str = "yabbadabbado"

    #str = prefix_result

    T = convert(Vector{Int8}, collect(str))
    n = length(T)
    push!(T, 0)
    push!(T, 0)
    push!(T, 0)

    SA = Vector{Int32}(undef, n)

    StringAlgorithms.suffix_array(T, SA, n, 257)

    for i in 2:n
        @test str[begin + SA[i-1]:end]< str[begin + SA[i]:end]
    end

#    radix_pass!
    @test 1 == 1
end
