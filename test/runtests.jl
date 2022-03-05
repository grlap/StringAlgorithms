using Random
using StringAlgorithms
using Test

@testset "RadixSort" begin
    Random.seed!(0)

    prefix_result = randstring(8 * 1024 * 1024)

    str = "yabbadabbado"

    str = prefix_result

    T = convert(Vector{Int8}, collect(str))
    n = length(T)
    push!(T, 0)
    push!(T, 0)
    push!(T, 0)

    SA = Vector{Int32}(undef, n)

    @time StringAlgorithms.suffix_array(T, SA, n, 257)

    for i in 2:n
        @test view(str, (SA[i - 1] + 1):n) < view(str, (SA[i] + 1):n)
    end

    #    radix_pass!
    @test 1 == 1
end
