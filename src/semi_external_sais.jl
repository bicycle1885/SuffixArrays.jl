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

function invert(SA)
    n = length(SA)
    ISA = ArrayBuffer{Int}(n)
    for i in 1:n
        ISA[SA[i]+1] = i - 1
    end
    return ISA
end

function left_to_right!(A_lms_left, A_lms_right, count_l, n_l, s, t, σ)
    #A_L = ArrayBuffer{Int,true}(n_l)
    A_L = ArrayBuffer{Int}(n_l)
    # debug
    #A_L = zeros(Int, n_l)
    tmp = copy(count_l)
    i = li = ri = 1

    # sentinel: '$' (< 0)
    char′ = s[end]
    A_L[count_l[char′+1]] = endof(s)
    count_l[char′+1] += 1
    # Note: there is no need to move `li` index
    # because '$' is not a member of `A_lms_left`
    # <del>li += 1</del>

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

function right_to_left!(A_lms_left, A_lms_right, count_s, n_s, s, t, σ)
    A_S = ArrayBuffer{Int}(n_s)
    # debug
    #A_S = zeros(Int, n_s)
    tmp = copy(count_s)
    i = n_s
    ri = length(A_lms_right)
    li = length(A_lms_left)
    for char in σ-1:-1:0
        while i > count_s[char+1]
            k = A_S[i]
            @assert k > 0
            if k != 1
                if t[k-1]
                    # S-type
                    char′ = s[k-1]
                    A_S[count_s[char′+1]] = k - 1
                    count_s[char′+1] -= 1
                else
                    A_lms_left[li] = k
                    li -= 1
                end
            end
            i -= 1
        end
        while ri ≥ 1 && s[A_lms_right[ri]] == char
            k = A_lms_right[ri]
            ri -= 1
            if k != 1
                @assert t[k-1]
                char′ = s[k-1]
                A_S[count_s[char′+1]] = k - 1
                count_s[char′+1] -= 1
            end
        end
    end
    copy!(count_s, tmp)
    return A_S
end

# a thin wrapper of sequence object
immutable Sequence{T} <: AbstractVector{UInt}
    seq::T
end

Base.getindex(seq::Sequence, i) = convert(UInt, seq.seq[i])
Base.length(seq::Sequence) = length(seq.seq)

function sais_se(s, SA, σ)
    # Step 1: Scan sequence and determine suffix types
    println("Step 1")
    tic()
    n = length(s)
    t = falses(n)
    count_s = zeros(Int, σ)
    count_l = zeros(Int, σ)
    count_lms = zeros(Int, σ + 1)
    count_l[s[n]+1] += 1
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
    A_lms_left = ArrayBuffer{Int}(n_lms)
    for i in 2:n
        if t[i] && !t[i-1]
            # LMS-type
            char = s[i]
            A_lms_left[count_lms[char+1]] = i
            count_lms[char+1] += 1
        end
    end
    copy!(count_lms, tmp)
    n_l = accumulate!(count_l)
    n_s = accumulate!(count_s, false)

    # Step 2: Induced sort (stage 1)
    toc()
    #quit()
    println("Step 2")
    tic()
    A_lms_right = ArrayBuffer{Int}(n_lms + 1)
    A_L = left_to_right!(A_lms_left, A_lms_right, count_l, n_l, s, t, σ)
    #Profile.print()
    finalize(A_L)

    # Step 3: Induced sort (stage 2)
    toc()
    println("Step 3")
    tic()
    init!(A_lms_left)
    A_S = right_to_left!(A_lms_left, A_lms_right, count_s, n_s, s, t, σ)
    finalize(A_S)

    # Step 4: Check uniqueness of LMS-type suffixes
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

    # Step 5: Naming LMS-type suffixes
    toc()
    println("Step 5")
    tic()
    S″ = ArrayBuffer{Int}(div(n, 2) + 1)
    n′ = n_lms
    σ′ = 0
    for i in 1:n′
        if B[i]
            σ′ += 1
        end
        j = A_lms_left[i]
        S″[div(j-1,2)+1] = σ′
    end
    j = 1
    S′ = ArrayBuffer{Int}(n′)
    for i in 1:endof(S″)
        if S″[i] > 0
            S′[j] = S″[i] - 1
            j += 1
        end
    end
    if !isempty(S′)
        @assert extrema(S′) == (0, σ′ - 1)
    end

    # Step 6: Recursion
    toc()
    println("Step 6")
    tic()
    # note: values of ISA′ and S′ start from zero
    if all(B)
        ISA′ = S′
    else
        if n′ ≤ div(1024^3, sizeof(Int))
            # bounded by 1GiB
            SA′ = Array(Int, n′)
            sais(S′, SA′, 0, n′, nextpow2(σ′), false)
        else
            SA′ = ArrayBuffer{Int}(n′)
            sais_se(S′, SA′, σ′)
        end
        ISA′ = invert(SA′)
    end
    init!(A_lms_left)
    i = j = 1
    while (i = find_next_LMS(t, i)) > 0
        A_lms_left[ISA′[j]+1] = i
        j += 1
    end

    # Step 7: Induced sort (stage 1)
    toc()
    println("Step 7")
    tic()
    init!(A_lms_right)
    A_L = left_to_right!(A_lms_left, A_lms_right, count_l, n_l, s, t, σ)

    # Step 8: Induced sort (stage 2)
    toc()
    println("Step 8")
    tic()
    init!(A_lms_left)
    A_S = right_to_left!(A_lms_left, A_lms_right, count_s, n_s, s, t, σ)

    # Step 9: Write output
    toc()
    println("Step 9")
    tic()
    i = li = si = 1
    for char in 0:σ-1
        # L-type suffixes come first
        while li ≤ n_l && s[A_L[li]] == char
            SA[i] = A_L[li] - 1
            i += 1
            li += 1
        end
        # then S-type ones
        while si ≤ n_s && s[A_S[si]] == char
            SA[i] = A_S[si] - 1
            i += 1
            si += 1
        end
    end
    finalize(A_L)
    finalize(A_S)
    toc()

    return SA
end
