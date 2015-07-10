# Timo Beller, et al. "Space-Efficient Construction
# of the Burrows-Wheeler Transform" (2013)

include("sais.jl")

type ArrayBuffer
    path::String
    io::IO
    function ArrayBuffer()
        path, io = mktemp()
        buf = new(path, io)
        finalizer(buf, x -> rm(x.path))
        buf
    end
end


Base.write(ab::ArrayBuffer, x) = write(ab.io, x)

type AlphabetCounter
    count::Vector{Int}
    AlphabetCounter(σ) = new(zeros(Int, σ + 1))
end

add!(c::AlphabetCounter, ch) = c.count[ch+2] += 1

type Bucket
    count::Vector{Int}
end

Bucket(c::AlphabetCounter) = Bucket(cumsum(c.count))

Base.getindex(b::Bucket, ch) = b.count[ch+1] + 1
Base.copy(b::Bucket) = Bucket(copy(b.count))
move_right!(b::Bucket, ch) = b.count[ch+1] += 1

# s: sequence
# sa: suffix array
# σ: alphabet size
function sais_se(s, sa, σ; dir=tempdir())
    n = length(s)

    # Step 1
    C   = AlphabetCounter(σ)
    C_s = AlphabetCounter(σ)
    C_l = AlphabetCounter(σ)
    C_lms = AlphabetCounter(σ)
    add!(C, s[n])
    lms_poss = falses(n)
    t = falses(n)
    t[n] = true
    n_lms = 0
    n_l = 0
    n_s = 1
    s_type_prev = true
    for i in n-1:-1:1
        add!(C, s[i])
        s_type = s[i] == s[i+1] ? s_type_prev : s[i] < s[i+1]
        if s_type  # S-type
            add!(C_s, s[i])
            n_s += 1
            t[i] = true
        else  # L-type
            add!(C_l, s[i])
            n_l += 1
            t[i] = false
            if s_type_prev
                add!(C_lms, s[i+1])
                n_lms += 1
                lms_poss[i+1] = true
            end
        end
        s_type_prev = s_type
    end
    C   = Bucket(C)
    C_s = Bucket(C_s)
    C_l = Bucket(C_l)
    C_lms = Bucket(C_lms)

    # TODO: use disk
    A_lms_left = zeros(Int, n_lms)
    tmp = copy(C_lms)
    for i in 1:n
        if lms_poss[i]
            A_lms_left[C_lms[s[i]]] = i
            move_right!(C_lms, s[i])
        end
    end
    C_lms = tmp
    @show A_lms_left

    # Step 2
    A_l = zeros(Int, n_l)
    # point to the current character in A_l
    c_l = 0
    # point to the current character in A_lms
    c_lms = 0
    tmp = copy(C_l)
    for i in 1:n_l
        @assert countnz(A_l) ≤ n_lms
        # fill LMS-type suffixes to A_l when entering a new bucket
        if i == C_l[c_l]
            c_l, c_lms = fill_suffixes!(A_l, i, A_lms_left, C_l, c_l, C_lms, c_lms, s, σ)
        end
        j = A_l[i]
        @assert j > 0
        if j > 1 && isL(t, j - 1)
            k = C_l[s[j-1]]
            A_l[k] = j - 1
            A_l[i] = 0  # remove
            move_right!(C_l, s[j-1])
        end
    end
    @show A_l
    C_l = tmp

    A_lms_right = Array(Int, n_lms)
    j = 1
    for i in 1:n_l
        if A_l[i] > 0
            A_lms_right[j] = A_l[i]
            j += 1
        end
    end
    @show A_lms_right

    # Step 3
    A_s = zeros(Int, n_s)
    c_s = 0
    c_lms = 0
    tmp = copy(C_s)
    for i in 1:n_s
        @assert countnz(A_s) ≤ n_lms
        # fill LMS-type suffixes to A_l when entering a new bucket
        if i == C_s[c_s]
            c_s, c_lms = fill_suffixes!(A_s, i, A_lms_right, C_s, c_s, C_lms,
            c_lms, s, σ, false)
        end

        j = A_s[i]
        @assert j > 0
        if j == 1
            A_s[end] = n
            A_s[i] = 0
            move_right!(C_s, s[n])
        elseif isS(t, j - 1)
            k = C_s[s[j-1]]
            A_s[k] = j - 1
            A_s[i] = 0  # remove
            move_right!(C_s, s[j-1])
        end
    end
    reverse!(A_s)
    C_s = tmp
    @show A_s

    j = 1
    for i in 1:n_s
        if A_s[i] > 0
            A_lms_left[j] = A_s[i]
            j += 1
        end
    end
    @show A_lms_left

    # Step 4
    B = falses(n_lms)
    B[1] = true
    for i in 2:n_lms
        lo₁ = A_lms_left[i-1]
        lo₂ = A_lms_left[i  ]
        hi₁ = findnext(lms_poss, lo₁ + 1)
        hi₂ = findnext(lms_poss, lo₂ + 1)
        if (len = hi₁ - lo₁) == hi₂ - lo₂
            # two LMS substrings have the same length
            for d in 0:len
                if s[lo₁+d] != s[lo₂+d]
                    B[i] = true
                    break
                end
            end
        else
            B[i] = true
        end
    end
    @show B
    R = Int[]
    i = 1
    while i ≤ n_lms
        if B[i]
            while i < n_lms && !B[i+1]
                i += 1
            end
            lo = A_lms_left[i]
            if lo == n
                hi = findnext(lms_poss, 1)
            else
                hi = findnext(lms_poss, lo + 1)
            end
            push!(R, hi)
        end
        i += 1
    end
    @show R
    S′ = zeros(Int, div(n, 2))
    name = 0
    for i in 1:n_lms
        if B[i]
            name += 1
        end
        j = A_lms_left[i]
        S′[div(j, 2)] = name
    end
    filter!(x -> x > 0, S′)
    @show S′

    # Step 5
    n′ = length(S′)
    SA′ = Array(Int, n′)
    σ′ = name
    if all(B)
        copy!(SA′, S′)
    else
        if n′ < 10_000
            sais(S′, SA′, 0, length(S′), nextpow2(σ′), false)
            SA′ += 1
        else
            sais_se(S′, SA′, σ′)
        end
    end
    @show SA′

    # Step 6
    ISA′ = Array(Int, n′)
    for i in 1:n′
        ISA′[SA′[i]] = i
    end
    @show ISA′
    i = 0
    j = 0
    while (i = findnext(lms_poss, i + 1)) > 0
        j += 1
        A_lms_left[ISA′[j]] = i
    end
    @show A_lms_left

    # Step 7
    fill!(A_l, 0)
    # point to the current character in A_l
    c_l = 0
    # point to the current character in A_lms
    c_lms = 0
    for i in 1:n_l
        # fill LMS-type suffixes to A_l when entering a new bucket
        if i == C_l[c_l]
            c_l, c_lms = fill_suffixes!(A_l, i, A_lms_left, C_l, c_l, C_lms, c_lms, s, σ)
        end
        j = A_l[i]
        @assert j > 0
        if j > 1 && isL(t, j - 1)
            k = C_l[s[j-1]]
            A_l[k] = j - 1
            #A_l[i] = 0  # remove
            move_right!(C_l, s[j-1])
        end
    end
    @show A_l

    # Step 8
    fill!(A_s, 0)
    c_s = 0
    c_lms = 0
    for i in 1:n_s
        # fill LMS-type suffixes to A_l when entering a new bucket
        if i == C_s[c_s]
            c_s, c_lms = fill_suffixes!(A_s, i, A_lms_right, C_s, c_s, C_lms,
            c_lms, s, σ, false)
        end

        j = A_s[i]
        @assert j > 0
        if j == 1
            A_s[end] = n
            #A_s[i] = 0
            move_right!(C_s, s[n])
        elseif isS(t, j - 1)
            k = C_s[s[j-1]]
            A_s[k] = j - 1
            #A_s[i] = 0
            move_right!(C_s, s[j-1])
        end
    end
    reverse!(A_s)
    @show A_s

    c = 0
    i = 0
    i_l = 1
    i_s = 1
    while c < σ
        if i_l ≤ endof(A_l) && s[A_l[i_l]] == c
            sa[i+=1] = A_l[i_l]
            i_l += 1
        elseif i_s ≤ endof(A_s) && s[A_s[i_s]] == c
            sa[i+=1] = A_s[i_s]
            i_s += 1
        else
            c += 1
        end
    end
    sa
end

function fill_suffixes!(A, i, A_lms, C, c, C_lms, c_lms, s, σ, forward=true)
    @assert i == C[c]
    while c < σ && i == C[c]
        c += 1
    end
    c -= 1
    while c_lms ≤ c
        lo = C_lms[c_lms]
        hi = c_lms < σ ? C_lms[c_lms+1] - 1 : length(s) + 1
        for j in lo:hi
            if forward
                k = A_lms[j]
            else
                # backward
                k = A_lms[end-j+1]
            end
            A[C[s[k-1]]] = k - 1
            move_right!(C, s[k-1])
        end
        c_lms += 1
    end
    c, c_lms
end

function fill_A_l!(A_l, C_l, C_lms, n_lms)
    n_l = length(A_l)
    # point to the current character in A_l
    c_l = 0
    # point to the current character in A_lms
    c_lms = 0
    tmp = copy(C_l)
    for i in 1:n_l
        @assert countnz(A_l) ≤ n_lms
        # fill LMS-type suffixes to A_l when entering a new bucket
        if i == C_l[c_l]
            c_l, c_lms = fill_suffixes!(A_l, i, A_lms_left, C_l, c_l, C_lms, c_lms, s, σ)
        end
        j = A_l[i]
        @assert j > 0
        if j > 1 && isL(t, j - 1)
            k = C_l[s[j-1]]
            A_l[k] = j - 1
            A_l[i] = 0  # remove
            move_right!(C_l, s[j-1])
        end
    end
end

isS(t, i) =  t[i]
isL(t, i) = !t[i]

#sais_se([0x01, 0x00, 0x00, 0x01, 0x00, 0x00, 0x01], Int[], 2)
ama = [1, 2, 1, 2, 2, 2, 1, 4, 1, 4, 2, 1, 4, 1, 4, 4, 1, 1, 3, 1, 0]
n = length(ama)
sa = zeros(Int, n)
sais_se(ama, sa, 5)
@show sa
sa = zeros(Int, n)
sais(ama, sa, 0, n, nextpow2(5), false)
sa += 1
@show sa
