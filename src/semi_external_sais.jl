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
    s_type = t[start]
    for i in start+1:endof(t)
        if !s_type && t[i]
            return i
        end
        s_type = t[i]
    end
    return 0
end

function left_to_right!(A_lms_left, A_lms_right, count_l, n_l, s, t, σ)
    A_L = ArrayBuffer{Int}(n_l)
    # debug
    #A_L = zeros(Int, n_l)
    tmp = copy(count_l)
    i = li = ri = 1
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
    end
    copy!(count_l, tmp)
    return A_L
end

function right_to_left!(A_lms_left, A_lms_right, count_s, n_s, n_lms, s, t, σ)
    # debug
    #A_S = zeros(Int, n_s)
    A_S = ArrayBuffer{Int}(n_s)
    tmp = copy(count_s)
    i = n_s
    ri = li = n_lms
    for char in σ-1:-1:0
        while i > count_s[char+1]
            k = A_S[i]
            @assert k > 0
            if k == 1
                char′ = s[end]
                A_S[count_s[char′+1]] = length(s)
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
            k = j == 1 ? length(s) : j - 1
            @assert t[k]
            char′ = s[k]
            A_S[count_s[char′+1]] = k
            count_s[char′+1] -= 1
        end
    end
    copy!(count_s, tmp)
    return A_S
end

function sais_se(s, SA, σ)
    # Step 1
    println("Step 1")
    tic()
    n = length(s)
    t = falses(n)
    t[n] = true
    count_s = zeros(Int, σ)
    count_l = zeros(Int, σ)
    count_lms = zeros(Int, σ + 1)
    count_s[s[n]+1] += 1
    for i in n-1:-1:1
        char = s[i]
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
    n_lms = accumulate!(count_lms)

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
    copy!(count_lms, tmp)
    seek(lms_positions, 0)
    A_lms_left = Array(Int, n_lms)
    for _ in 1:n_lms
        i = read(lms_positions, Int)
        suf = read(lms_positions, Int)
        A_lms_left[suf] = i
    end
    close(lms_positions)
    rm("lms_positions")
    n_l = accumulate!(count_l)
    n_s = accumulate!(count_s, false)

    # Step 2
    toc()
    println("Step 2")
    tic()
    A_lms_right = zeros(Int, n_lms)
    A_L = left_to_right!(A_lms_left, A_lms_right, count_l, n_l, s, t, σ)
    finalize(A_L)

    # Step 3
    toc()
    println("Step 3")
    tic()
    A_S = right_to_left!(A_lms_left, A_lms_right, count_s, n_s, n_lms, s, t, σ)
    finalize(A_S)

    # Step 4
    toc()
    println("Step 4")
    tic()
    B = trues(n_lms)
    for i in 2:n_lms
        lo  = A_lms_left[i]
        hi  = find_next_LMS(t, lo)
        lo′ = A_lms_left[i-1]
        hi′ = find_next_LMS(t, lo′)
        if hi - lo == hi′ - lo′
            while lo <= hi && s[lo] == s[lo′]
                lo  += 1
                lo′ += 1
            end
            if lo > hi
                # duplicated
                B[i] = false
            end
        end
    end

    # Step 5
    toc()
    println("Step 5")
    tic()
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

    # Step 6
    toc()
    println("Step 6")
    tic()
    ISA′ = Array(Int, n′)
    if all(B)
        copy!(ISA′, S′)
    else
        SA′ = Array(Int, n′)
        sais(S′, SA′, 0, n′, nextpow2(σ′), false)
        for i in 1:n′
            ISA′[SA′[i]+1] = i
        end
        i = j = 1
        while (i = find_next_LMS(t, i)) > 0
            A_lms_left[ISA′[j]] = i
            j += 1
        end
    end
    #@show ISA′

    # Step 7
    toc()
    println("Step 7")
    tic()
    A_L = left_to_right!(A_lms_left, A_lms_right, count_l, n_l, s, t, σ)

    # Step 8
    toc()
    println("Step 8")
    tic()
    A_S = right_to_left!(A_lms_left, A_lms_right, count_s, n_s, n_lms, s, t, σ)

    # Step 9
    toc()
    println("Step 9")
    tic()
    i = li = si = 1
    for char in 0:σ-1
        # L-type suffixes come first
        while li ≤ n_l && s[A_L[li]] == char
            SA[i] = A_L[li]
            i += 1
            li += 1
        end
        # then S-type ones
        while si ≤ n_s && s[A_S[si]] == char
            SA[i] = A_S[si]
            i += 1
            si += 1
        end
    end
    finalize(A_L)
    finalize(A_S)
    toc()

    return SA
end

#=
let
    for s in ["amammmasasmasassaara", "abra", "baba", "dabraaadad", "abaadad",
        "aaaaa", "acacac", "abababa"]
        println(s)
        s = vcat(s.data, 0)
        n = length(s)
        SA = Array(Int, n)
        sais(s, SA, 0, n, 256, false)
        SA += 1
        SA1 = ArrayBuffer{Int}(n)
        sais_se(s, SA1, 128)
        @assert SA == SA1
        finalize(SA1)
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
        SA1 = ArrayBuffer{Int}(n)
        sais_se(s, SA1, σ)
        @assert SA == SA1
        finalize(SA1)
    end
end

let
    srand(12345)
    n = 1
    for _ in 1:100
        s = randstring(n)
        #println(s)
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

let
    srand(12345)
    n = 10000
    s = rand(1:4, n)
    push!(s, 0)
    n += 1
    SA = Array(Int, n)
    sais(s, SA, 0, n, 256, false)
    SA += 1
    SA1 = ArrayBuffer{Int}(n)
    sais_se(s, SA1, 5)
    @assert SA == SA1
end
=#
