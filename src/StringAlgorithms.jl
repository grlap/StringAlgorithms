module StringAlgorithms

greet() = print("Hello World!")

@inline leq(a1, a2, b1, b2) = a1 < b1 || a1 == b1 && a2 <= b2

@inline leq(a1, a2, a3, b1, b2, b3) = a1 < b1 || a1 == b1 && leq(a2, a3, b2, b3)

function radix_pass(
    a::Vector{Int32},
    b::Vector{Int32},
    r::AbstractVector{T},
    n::Int32,
    K::Int64) where {T}

    c = zeros(Int32, K + 1)

    for i in 1:n
        c[begin + r[begin + a[i]]] += 1
    end

    sum::Int32 = 0
    for i in 0:K
        t = c[begin + i]
        c[begin + i] = sum
        sum += t
    end

    for i in 1:n
        b[begin + c[begin + r[begin + a[i]]]] = a[i]
        c[begin + r[begin + a[i]]] += 1
    end

    @show "radix_pass_result", a[1:n], b[1:n]
end

function suffix_array(
    T::AbstractVector{V},
    SA::Vector{Int32},
    n,
    K) where {V}

    @show "[] => suffix_array", T, n, K

    n0::Int32 = (n + 2) ÷ 3
    n1::Int32 = (n + 1) ÷ 3
    n2::Int32 = n ÷ 3
    n02 = n0 + n2

    @show n02

    R = Vector{Int32}(undef, n02 + 3)
    R[begin + n02:begin + n02 + 2] = [0,0,0]

    @show R

    SA12 = Vector{Int32}(undef, n02 + 3)
    SA12[begin + n02:begin + n02 + 2] = [0,0,0]

    R0 = Vector{Int32}(undef, n0)
    SA0 = Vector{Int32}(undef, n0)

    @show n0, n1, n2
    @show n + n0 - n1

    j::Int32 = 0
    for i in 0:n+n0-n1-1
        if i % 3 != 0
            R[begin + j] = i
            j += 1
        end
    end

    @show T
    @show R

    radix_pass(R, SA12, T[3:end], n02, K)
    radix_pass(SA12, R, T[2:end], n02, K)
    radix_pass(R, SA12, T[1:end], n02, K)

    @show SA12

    name = 0
    c0 = -1
    c1 = -1
    c2 = -1

    @show "before", R
    for i in 0:n02 - 1
        if T[begin + SA12[begin + i]] != c0 ||
           T[begin + SA12[begin + i] + 1] != c1 ||
           T[begin + SA12[begin + i] + 2] != c2

           name += 1
           c0 = T[begin + SA12[begin + i]]
           c1 = T[begin + SA12[begin + i] + 1]
           c2 = T[begin + SA12[begin + i] + 2]
        end

        if SA12[begin + i] % 3 == 1
            R[begin + SA12[begin + i] ÷ 3] = name
        else
            R[begin + SA12[begin + i] ÷ 3 + n0] = name
        end
    end
    @show "afer", R

    @show "before_suffix_array", R, c0, c1, c2, name
    if name < n02
        suffix_array(R, SA12, n02, name)
        @show "afer_suffix_array", SA12

        for i in 1:n02
            R[begin + SA12[i]] = i
        end
    else
        for i in 1:n02
            @show "mapping", R[i] - 1, i - 1

            SA12[begin + R[i] - 1] = i - 1
        end
    end

    j = 0
    for i in 1:n02
        if SA12[i] < n0
            R0[begin + j] = 3 * SA12[i]
            j += 1
        end
    end

    @show "before_radix", R0
    radix_pass(R0, SA0, T, n0, K)
    @show "after_radix", SA0

    @show "merge 0 prefix", SA0, R0
    @show "merge 12 prefix", SA12, R

    k = 0
    p = 0
    t = n0 - n1
    while k < n
        i = SA12[begin + t] < n0 ? SA12[begin + t] * 3 + 1 : (SA12[begin + t] - n0) * 3 +2
        j = SA0[begin + p]

        result = SA12[begin + t] < n0 ?
           leq(T[begin + i], R[begin + SA12[begin + t] + n0], T[begin + j], R[begin + j ÷ 3]) :
           leq(T[begin + i], T[begin + i + 1], R[begin + SA12[begin + t] - n0 + 1], T[begin + j], T[begin + j + 1], R[begin + j ÷ 3 + n0])

        if result == true
            SA[begin + k] = i
            t += 1

            if t == n02
                @show "end n02", k, p
                k += 1
                while p < n0
                    SA[begin + k] = SA0[begin + p]
                    p += 1
                    k += 1
                end
            end
        else
            SA[begin + k] = j
            p += 1

            if p == n0
                @show "end n0", k
                k += 1

                while t < n02
                    SA[begin + k] = SA12[begin + t] < n0 ? SA12[begin + t] * 3 + 1 : (SA12[begin + t] - n0) * 3 + 2
                    t += 1
                    k += 1
                end
            end
        end

        k += 1
    end

    @show "end_suffix_array", SA
end

end # module
