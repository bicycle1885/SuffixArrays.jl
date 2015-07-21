type Block{T}
    # immutable
    size::Int
    # mutable
    used::Int
    offset::Int
    io::IO
    buffer::Vector{T}
end

function Base.getindex{T}(block::Block{T}, i::Integer)
    buffersize = length(block.buffer)
    j = i - block.offset
    if 1 ≤ j ≤ buffersize
        return block.buffer[j]
    end
    # write buffer
    seekend(block.io)
    for j in 1:buffersize
        write(block.io, j + block.offset)
        write(block.io, block.buffer[j])
    end
    flush(block.io)
    # read buffer
    pairsize = sizeof(Int) + sizeof(T)
    block.offset = div(i - 1, buffersize) * buffersize
    n = min(buffersize, block.used - block.offset)
    s = IntSet(1:n)
    while !isempty(s) && position(block.io) != 0
        skip(block.io, -pairsize)
        idx = read(block.io, Int)
        elm = read(block.io, T)
        skip(block.io, -pairsize)
        j = idx - block.offset
        # NOTE: adding necessary condition may make it faster
        # i.e. if 1 ≤ j ≤ buffersize && j ∈ s
        if j ∈ s
            block.buffer[j] = elm
            delete!(s, j)
        end
    end
    j = i - block.offset
    # if not yet initialized, return zero(T)
    return j ∈ s ? zero(T) : block.buffer[j]
end

function Base.setindex!{T}(block::Block{T}, x::T, i::Integer)
    j = i - block.offset
    if 1 ≤ j ≤ length(block.buffer)
        block.buffer[j] = x
    else
        seekend(block.io)
        write(block.io, i)
        write(block.io, x)
        flush(block.io)
    end
    return x
end

select_n_blocks(len) = len == 0 ? 0 : len ≤ 32 ? 1 : 32
function select_blocksize(len, n_blocks)
    if n_blocks == 0
        return 64
    end
    return max(64, div(len - 1, n_blocks) + 1)
end
select_buffersize(blocksize) = div(blocksize - 1, 16) + 1

type ArrayBuffer{T} <: AbstractVector{T}
    len::Int
    dirname::String
    blocksize::Int
    blocks::Vector{Block}
    # space: buffersize * n_blocks * sizeof(T) bytes
    function ArrayBuffer{T}(::Type{T}, len::Integer, parent::String,
        n_blocks=select_n_blocks(len),
        blocksize=select_blocksize(len, n_blocks),
        buffersize=select_buffersize(blocksize))
        @assert 0 ≤ n_blocks
        @assert 1 ≤ buffersize ≤ blocksize
        @assert 0 ≤ len ≤ blocksize * n_blocks
        dirname = mktempdir(parent)
        arr = new(len, dirname, blocksize)
        finalizer(arr, x -> begin
            for block in x.blocks
                if isopen(block.io)
                    close(block.io)
                end
            end
            rm(dirname, recursive=true)
        end)
        arr.blocks = Block[]
        l = len
        for i in 1:n_blocks
            _, io = mktemp(dirname)
            buffer = zeros(T, buffersize)
            used = l ≥ blocksize ? blocksize : l
            l -= used
            block = Block(blocksize, used, 0, io, buffer)
            push!(arr.blocks, block)
        end
        return arr
    end
end

if VERSION >= v"0.4-"
    function Base.call{T}(::Type{ArrayBuffer{T}}, len::Int, n_blocks=select_n_blocks(len))
        ArrayBuffer{T}(T, len, pwd(), n_blocks)
    end
end

function Base.getindex(arr::ArrayBuffer, i::Integer)
    j, o = divrem(i - 1, arr.blocksize)
    return arr.blocks[j+1][o+1]
end

function Base.setindex!(arr::ArrayBuffer, x, i::Integer)
    j, o = divrem(i - 1, arr.blocksize)
    arr.blocks[j+1][o+1] = x
end

Base.size(arr::ArrayBuffer) = (arr.len,)
Base.length(arr::ArrayBuffer) = arr.len
Base.endof(arr::ArrayBuffer)  = arr.len

function Base.fill!(arr::ArrayBuffer, x)
    for i in 1:length(arr)
        arr[i] = x
    end
    arr
end

#=
let
    # array is zero-initialized
    for n in [0, 5, 10, 50, 100, 500, 1000]
        arr = ArrayBuffer{Int}(n)
        for i in 1:n
            @assert arr[i] == zero(Int)
        end
    end
end

let
    for n in [0, 5, 10, 50, 100, 500, 1000]
        arr = ArrayBuffer{Int}(n)
        for i in 1:n
            arr[i] = 1
        end
        for i in 1:n
            @assert arr[i] == 1
        end
        finalize(arr)
    end
end

let
    srand(1234)
    arr = Array{Int}(1000)
    barr = ArrayBuffer{Int}(1000)
    for i in 1:length(arr)
        x = rand(1:100)
        arr[i] = x
        barr[i] = x
    end
    for _ in 1:1000
        i = rand(1:1000)
        x = rand(1:100)
        arr[i] = x
        barr[i] = x
    end
    for i in 1:1000
        @assert arr[i] == barr[i]
    end
    for i in 1000:-1:1
        @assert arr[i] == barr[i]
    end
    fill!(arr, 0)
    fill!(barr, 0)
    for i in 1:1000
        @assert arr[i] == barr[i] == 0
    end
    for i in 1000:-1:1
        @assert arr[i] == barr[i] == 0
    end
end

let
    srand(1234)
    arr = ArrayBuffer{Int}(100_000, 53)
    for i in 1:length(arr)
        arr[i] = i
    end
    for _ in 1:1000
        arr[rand(1:100_000)] = 0
    end
    for i in 1:length(arr)
        arr[i] = i
    end
    for i in 1:length(arr)
        x = arr[i]
        @assert i == x
    end
end
=#
