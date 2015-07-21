# Timo Beller, et al. "Space-Efficient Construction
# of the Burrows-Wheeler Transform" (2013)

include("sais.jl")
include("arraybuffer.jl")

function accumulate!(count, startpoint=true)
    σ = length(count)
    sum = 0
    if startpoint
        for char in 0:σ-1
            tmp = count[char+1]
            count[char+1] = sum + 1
            sum += tmp
        end
    else
        for char in 0:σ-1
            tmp = count[char+1]
            count[char+1] = sum + tmp
            sum += tmp
        end
    end
    sum
end

function find_next_LMS(t, start)
    #@assert start == 1 || t[start] && !t[start-1] "first index or LMS position"
    #l_type = false
    l_type = !t[start]
    for i in start+1:endof(t)
        if l_type && t[i]
            return i
        end
        l_type = !t[i]
    end
    return 0
end

function sais_se(s, SA, σ)
    # Step 1
    n = length(s)
    t = falses(n)
    t[n] = true
    count = zeros(Int, σ)
    count_s = zeros(Int, σ + 1)
    count_l = zeros(Int, σ + 1)
    count_lms = zeros(Int, σ + 1)
    count[s[n]+1] += 1
    count_s[s[n]+1] += 1
    for i in n-1:-1:1
        char = s[i]
        count[char+1] += 1
        t[i] = s[i] == s[i+1] ? t[i+1] : s[i] < s[i+1]
        if t[i]
            # S-type
            count_s[char+1] += 1
        else
            # L-type
            count_l[char+1] += 1
            if t[i+1]
                # i+1 is LMS-type
                count_lms[s[i+1]+1] += 1
            end
        end
    end
    accumulate!(count)
    n_s = accumulate!(count_s, false)
    n_l = accumulate!(count_l)
    n_lms = accumulate!(count_lms)
    hoge = copy(count_lms)

    tmp = copy(count_lms)
    lms_positions = open("lms_positions", "w+")
    for i in 2:n
        if i == n || t[i] && !t[i-1]
            # LMS-type
            char = s[i]
            write(lms_positions, i)
            write(lms_positions, count_lms[char+1])
            count_lms[char+1] += 1
        end
    end
    count_lms = tmp
    seek(lms_positions, 0)
    A_lms_left = Array(Int, n_lms)
    for _ in 1:n_lms
        i = read(lms_positions, Int)
        suf = read(lms_positions, Int)
        A_lms_left[suf] = i
    end
    rm("lms_positions")

    # Step 2
    #A_L = ArrayBuffer{Int}(n_l)
    #A_lms_right = ArrayBuffer{Int}(n_lms)
    # debug
    A_L = zeros(Int, n_l)
    A_lms_right = zeros(Int, n_lms)
    count_lms′ = similar(count_lms)
    i = li = ri = 1
    count_l′ = copy(count_l)
    for char in 0:σ-1
        # advance suffixes in A_L and write suffixes to A_lms_right
        while i < count_l[char+1]
            k = A_L[i]
            @assert k > 0
            if k >= 2 && !t[k-1]
                # L-type
                char′ = s[k-1]
                A_L[count_l[char′+1]] = k - 1
                count_l[char′+1] += 1
            else
                A_lms_right[ri] = k
                ri += 1
            end
            i += 1
        end
        # read suffixes from A_lms_left to A_L
        while li ≤ endof(A_lms_left) && s[A_lms_left[li]] == char
            j = A_lms_left[li]
            li += 1
            @assert !t[j-1]
            char′ = s[j-1]
            A_L[count_l[char′+1]] = j - 1
            count_l[char′+1] += 1
        end
        # TODO: in-place update would be OK
        count_lms′[char+1] = ri - 1
    end
    count_lms′[σ+1] = n_lms
    count_lms = count_lms′
    count_l = count_l′
    #@show A_L A_lms_right

    # Step 3
    # debug
    A_S = zeros(Int, n_s)
    count_s′ = copy(count_s)
    i = n_s
    ri = n_lms
    li = n_lms
    for char in σ-1:-1:0
        while i > count_s[char+1]
            k = A_S[i]
            @assert k > 0
            if k == 1
                char′ = s[end]
                A_S[count_s[char′+1]] = n
                count_s[char′+1] -= 1
            elseif t[k-1]
                # S-type
                char′ = s[k-1]
                A_S[count_s[char′+1]] = k - 1
                count_s[char′+1] -= 1
            else
                A_lms_left[li] = k
                li -= 1
            end
            i -= 1
        end
        while ri ≥ 1 && s[A_lms_right[ri]] == char
            j = A_lms_right[ri]
            ri -= 1
            k = j == 1 ? n : j - 1
            @assert t[k]
            char′ = s[k]
            A_S[count_s[char′+1]] = k
            count_s[char′+1] -= 1
        end
    end
    count_s = count_s′

    # Step 4
    B = trues(n_lms)
    for i in 2:n_lms
        lo = A_lms_left[i]
        hi = find_next_LMS(t, lo)
        lo′ = A_lms_left[i-1]
        hi′ = find_next_LMS(t, lo′)
        if hi - lo == hi′ - lo′
            while lo <= hi && s[lo] == s[lo′]
                lo += 1
                lo′ += 1
            end
            if lo > hi
                B[i] = false
            end
        end
    end

    S′ = zeros(Int, div(n, 2) + 1)
    n′ = n_lms
    σ′ = 0
    for i in 1:n′
        if B[i]
            σ′ += 1
        end
        j = A_lms_left[i]
        S′[div(j-1,2)+1] = σ′
    end

    j = 1
    for i in 1:endof(S′)
        if S′[i] > 0
            S′[j] = S′[i]
            j += 1
        end
    end
    resize!(S′, n′)

    SA′ = Array(Int, n′)
    ISA′ = Array(Int, n′)
    if all(B)
        copy!(ISA′, S′)
    else
        sais(S′, SA′, 0, n′, nextpow2(σ′), false)
        SA′ += 1
        for i in 1:n′
            ISA′[SA′[i]] = i
        end
        i = j = 1
        while (i = find_next_LMS(t, i)) > 0
            A_lms_left[ISA′[j]] = i
            j += 1
        end
    end
    #@show ISA′

    # Step 7
    A_L = zeros(Int, n_l)
    A_lms_right = zeros(Int, n_lms)
    count_lms′ = similar(count_lms)
    i = 1
    li = 1
    ri = 1
    count_l′ = copy(count_l)
    for char in 0:σ-1
        # advance suffixes in A_L and write suffixes to A_lms_right
        while i < count_l[char+1]
            k = A_L[i]
            @assert k > 0
            if k >= 2 && !t[k-1]
                # L-type
                char′ = s[k-1]
                A_L[count_l[char′+1]] = k - 1
                count_l[char′+1] += 1
            else
                A_lms_right[ri] = k
                ri += 1
            end
            i += 1
        end
        # read suffixes from A_lms_left to A_L
        while li ≤ endof(A_lms_left) && s[A_lms_left[li]] == char
            j = A_lms_left[li]
            li += 1
            @assert !t[j-1]
            char′ = s[j-1]
            A_L[count_l[char′+1]] = j - 1
            count_l[char′+1] += 1
        end
        count_lms′[char+1] = ri - 1
    end
    count_lms′[σ+1] = n_lms
    count_lms = count_lms′
    count_l = count_l′

    # Step 8
    A_S = zeros(Int, n_s)
    count_s′ = copy(count_s)
    i = n_s
    ri = n_lms
    li = n_lms
    for char in σ-1:-1:0
        while i > count_s[char+1]
            k = A_S[i]
            @assert k > 0
            if k == 1
                char′ = s[end]
                A_S[count_s[char′+1]] = n
                count_s[char′+1] -= 1
            elseif t[k-1]
                # S-type
                char′ = s[k-1]
                A_S[count_s[char′+1]] = k - 1
                count_s[char′+1] -= 1
            else
                A_lms_left[li] = k
                li -= 1
            end
            i -= 1
        end
        while ri ≥ 1 && s[A_lms_right[ri]] == char
            j = A_lms_right[ri]
            ri -= 1
            k = j == 1 ? n : j - 1
            @assert t[k]
            char′ = s[k]
            A_S[count_s[char′+1]] = k
            count_s[char′+1] -= 1
        end
    end

    i = 1
    li = 1
    si = 1
    for char in 0:σ-1
        while li ≤ n_l && s[A_L[li]] == char
            SA[i] = A_L[li]
            li += 1
            i += 1
        end
        while si ≤ n_s && s[A_S[si]] == char
            SA[i] = A_S[si]
            si += 1
            i += 1
        end
    end
    #@show SA
end

let
    for s in ["amammmasasmasassaara", "abra", "baba", "dabraaadad", "abaadad",
        "aaaaa"]
        println(s)
        s = vcat(s.data, 0)
        n = length(s)
        SA = Array(Int, n)
        sais(s, SA, 0, n, 256, false)
        SA += 1
        SA1 = ArrayBuffer{Int}(Int, n, ".")
        sais_se(s, SA1, 128)
        @assert SA == SA1
    end
end

let
    for (s, σ) in [([1,2,3,2,1],4), ([3,2,1,2,3,2,1],4)]
        println(s)
        s = vcat(s, 0)
        n = length(s)
        SA = Array(Int, n)
        sais(s, SA, 0, n, σ, false)
        SA += 1
        SA1 = ArrayBuffer{Int}(Int, n, ".")
        sais_se(s, SA1, σ)
        @assert SA == SA1
    end
end

let
    srand(12345)
    n = 2
    for _ in 1:100
        s = randstring(n)
        s = vcat(s.data, 0)
        n = length(s)
        SA = Array(Int, n)
        sais(s, SA, 0, n, 256, false)
        SA += 1
        SA1 = Array(Int, n)
        sais_se(s, SA1, 128)
        @assert SA == SA1
    end
end
