# Timo Beller, et al. "Space-Efficient Construction
# of the Burrows-Wheeler Transform" (2013)

include("sais.jl")
include("arraybuffer.jl")

function accumulate!(count, startpoint)
    σ = length(count)
    sum = 0
    for char in 0:σ-1
        tmp = count[char+1]
        count[char+1] = sum + ifelse(startpoint, 1, tmp)
        sum += tmp
    end
    sum
end

function scan!(t, S, σ)
    n = length(S)
    count_l = zeros(Int, σ)
    count_s = zeros(Int, σ)
    count_lms = zeros(Int, σ + 1)
    count_l[S[n]+1] += 1
    for i in n-1:-1:1
        char = S[i]
        t[i] = S[i] == S[i+1] ? t[i+1] : S[i] < S[i+1]
        if t[i]
            # S-type
            count_s[char+1] += 1
        else
            # L-type
            count_l[char+1] += 1
            if t[i+1]
                # i+1 is LMS-type
                count_lms[S[i+1]+1] += 1
            end
        end
    end
    n_l = accumulate!(count_l, true)
    n_s = accumulate!(count_s, false)
    n_lms = accumulate!(count_lms, true)
    (count_l, n_l), (count_s, n_s), (count_lms, n_lms)
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

function invert{T}(SA, ::Type{T})
    n = length(SA)
    ISA = ArrayBuffer{T}(n)
    for i in 1:n
        ISA[SA[i]+1] = i - 1
    end
    return ISA
end

function destroy!(a::Vector)
    empty!(a)
end

function left_to_right!{T}(A_lms_left::AbstractVector{T}, A_lms_right::AbstractVector{T}, count_l, n_l, s, t, σ)
    A_L = ArrayBuffer{T}(n_l)
    let count_l = copy(count_l)
        i = li = ri = 1

        # sentinel: '$' (< 0)
        char′ = s[end]
        A_L[count_l[char′+1]] = endof(s)
        count_l[char′+1] += 1
        # Note: there is no need to move `li` index
        # because '$' is not a member of `A_lms_left`
        # <del>li += 1</del>

        @inbounds for char in 0:σ-1
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
    end
    return A_L
end

function right_to_left!{T}(A_lms_left::AbstractVector{T}, A_lms_right::AbstractVector{T}, count_s, n_s, s, t, σ)
    A_S = ArrayBuffer{T}(n_s)
    let count_s = copy(count_s)
        i = n_s
        ri = length(A_lms_right)
        li = length(A_lms_left)
        @inbounds for char in σ-1:-1:0
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
    end
    return A_S
end

# a thin wrapper of sequence object, which encodes sequence values with UInt
immutable Sequence{S} <: AbstractVector{UInt}
    seq::S
end

@inbounds Base.getindex(seq::Sequence, i) = convert(UInt, seq.seq[i])
Base.length(seq::Sequence) = length(seq.seq)

function sais_se(S, SA, σ)
    n = length(S)
    if n ≤ typemax(UInt32)
        sais_se(S, SA, σ, UInt32)
    else
        sais_se(S, SA, σ, UInt64)
    end
    return SA
end

function sais_se{T<:Integer}(S, SA, σ, ::Type{T})
    # Step 1: Scan sequence and determine suffix types
    println("Step 1")
    tic()
    n = length(S)
    t = falses(n)
    (count_l, n_l), (count_s, n_s), (count_lms, n_lms) = scan!(t, S, σ)

    A_lms_left = ArrayBuffer{T}(n_lms)
    let count_lms = copy(count_lms)
        for i in 2:n
            if t[i] && !t[i-1]
                # LMS-type
                char = S[i]
                A_lms_left[count_lms[char+1]] = i
                count_lms[char+1] += 1
            end
        end
    end

    # Step 2: Induced sort (stage 1)
    toc()
    println("Step 2")
    tic()
    A_lms_right = ArrayBuffer{T}(n_lms + 1)
    A_L = left_to_right!(A_lms_left, A_lms_right, count_l, n_l, S, t, σ)
    destroy!(A_L)
    #finalize(A_L)
    gc()

    # Step 3: Induced sort (stage 2)
    toc()
    println("Step 3")
    tic()
    A_S = right_to_left!(A_lms_left, A_lms_right, count_s, n_s, S, t, σ)
    #finalize(A_S)
    destroy!(A_S)
    #finalize(A_lms_right)
    destroy!(A_lms_right)
    gc()

    # Step 4: Check uniqueness of LMS-type suffixes
    toc()
    println("Step 4")
    tic()
    B = trues(n_lms)
    for i in 2:n_lms
        l  = A_lms_left[i]
        u  = find_next_LMS(t, l)
        l′ = A_lms_left[i-1]
        u′ = find_next_LMS(t, l′)
        if u - l == u′ - l′
            # note: comparing backwards is much faster
            while l <= u && S[u] == S[u′]
                u  -= 1
                u′ -= 1
            end
            if l > u
                # duplicated
                B[i] = false
            end
        end
    end

    # Step 5: Naming LMS-type suffixes
    toc()
    println("Step 5")
    tic()
    S″ = ArrayBuffer{T}(div(n, 2) + 1)
    n′ = n_lms
    σ′ = 0
    for i in 1:n′
        if B[i]
            σ′ += 1
        end
        j = A_lms_left[i]
        S″[div(j-1,2)+1] = σ′
    end
    #finalize(A_lms_left)
    destroy!(A_lms_left)
    S′ = ArrayBuffer{T}(n′)
    let j = 1
        for i in 1:endof(S″)
            if S″[i] > 0
                S′[j] = S″[i] - 1
                j += 1
            end
        end
    end
    #finalize(S″)
    destroy!(S″)
    gc()
    #@assert isempty(S′) || extrema(S′) == (0, σ′ - 1)

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
            SA′ = Vector{Int}(n′)
            sais(S′, SA′, 0, n′, nextpow2(σ′), false)
        else
            SA′ = ArrayBuffer{T}(n′)
            sais_se(S′, SA′, σ′)
        end
        #finalize(S′)
        destroy!(S′)
        ISA′ = invert(SA′, T)
        destroy!(SA′)
    end
    A_lms_left = ArrayBuffer{T}(n_lms)
    let i = j = 1
        while (i = find_next_LMS(t, i)) > 0
            A_lms_left[ISA′[j]+1] = i
            j += 1
        end
    end
    destroy!(ISA′)
    gc()

    # Step 7: Induced sort (stage 1)
    toc()
    println("Step 7")
    tic()
    A_lms_right = ArrayBuffer{T}(n_lms + 1)
    A_L = left_to_right!(A_lms_left, A_lms_right, count_l, n_l, S, t, σ)
    gc()

    # Step 8: Induced sort (stage 2)
    toc()
    println("Step 8")
    tic()
    A_S = right_to_left!(A_lms_left, A_lms_right, count_s, n_s, S, t, σ)
    #finalize(A_lms_left)
    destroy!(A_lms_left)
    #finalize(A_lms_right)
    destroy!(A_lms_right)
    gc()

    # Step 9: Write output
    toc()
    println("Step 9")
    tic()
    let i = li = si = 1
        for char in 0:σ-1
            # L-type suffixes come first
            while li ≤ n_l && S[A_L[li]] == char
                SA[i] = A_L[li] - 1
                i += 1
                li += 1
            end
            # then S-type ones
            while si ≤ n_s && S[A_S[si]] == char
                SA[i] = A_S[si] - 1
                i += 1
                si += 1
            end
        end
    end
    #finalize(A_L)
    destroy!(A_L)
    #finalize(A_S)
    destroy!(A_S)
    gc()
    toc()

    return SA
end
