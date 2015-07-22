type Block{T}
    # immutable
    size::Int
    # mutable
    used::Int
    offset::Int
    io::IO
    isdirtybuffer::Bool
    buffer::Vector{T}
    dirtybit::BitVector
end

function Block{T}(::Type{T}, size, used, buffersize, io)
    buffer = zeros(T, buffersize)
    dirtybit = falses(buffersize)
    Block(size, used, 0, io, false, buffer, dirtybit)
end

function Base.getindex{T}(block::Block{T}, i::Integer)
    # invariant:
    #   i (absolute index in block) = j (relative index in buffer)
    #                               + block.offset (buffer offset in block)
    buffersize = length(block.buffer)
    j = i - block.offset
    if 1 ≤ j ≤ buffersize
        return block.buffer[j]
    end
    if block.isdirtybuffer
        seekend(block.io)
        # write dirty elements
        j = 0
        while (j = findnext(block.dirtybit, j + 1)) > 0
            write(block.io, j + block.offset)
            write(block.io, block.buffer[j])
        end
        flush(block.io)
    end
    # read buffer
    block.offset = div(i - 1, buffersize) * buffersize
    s = IntSet(1:min(buffersize, block.used - block.offset))
    seekstart(block.io)
    while !isempty(s) && !eof(block.io)
        idx = read(block.io, Int)
        elm = read(block.io, T)
        j = idx - block.offset
        if j ∈ s
            block.buffer[j] = elm
            delete!(s, j)
        end
    end
    fill!(block.dirtybit, false)
    block.isdirtybuffer = false
    j = i - block.offset
    # if not yet initialized, return zero(T)
    return j ∈ s ? zero(T) : block.buffer[j]
end

function Base.setindex!{T}(block::Block{T}, x::T, i::Integer)
    j = i - block.offset
    if 1 ≤ j ≤ length(block.buffer)
        block.buffer[j] = x
        block.dirtybit[j] = true
        block.isdirtybuffer = true
    else
        seekend(block.io)
        write(block.io, i)
        write(block.io, x)
        flush(block.io)
    end
    return x
end

# parameter selectors
select_n_blocks(len) = len == 0 ? 0 : len ≤ 16 ? 1 : 16
select_blocksize(len, n_blocks) = n_blocks == 0 ? 64 : max(64, div(len - 1, n_blocks) + 1)
select_buffersize(blocksize) = div(blocksize - 1, 16) + 1

type ArrayBuffer{T} <: AbstractVector{T}
    len::Int
    dirname::String
    blocksize::Int
    blocks::Vector{Block}
    # space: buffersize * n_blocks * sizeof(T) bytes
    function ArrayBuffer{T}(::Type{T}, len::Integer;
        parent::String=pwd(),
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
        arr.blocks = Array(Block, n_blocks)
        l = len
        for i in 1:n_blocks
            _, io = mktemp(dirname)
            used = l ≥ blocksize ? blocksize : l
            l -= used
            arr.blocks[i] = Block(T, blocksize, used, buffersize, io)
        end
        return arr
    end
end

if VERSION >= v"0.4-"
    function Base.call{T}(::Type{ArrayBuffer{T}}, len::Int; kwargs...)
        ArrayBuffer{T}(T, len; kwargs...)
    end
end

function Base.getindex(arr::ArrayBuffer, i::Integer)
    j, k = divrem(i - 1, arr.blocksize)
    return arr.blocks[j+1][k+1]
end

function Base.setindex!(arr::ArrayBuffer, x, i::Integer)
    j, k = divrem(i - 1, arr.blocksize)
    arr.blocks[j+1][k+1] = x
end

Base.size(arr::ArrayBuffer) = (arr.len,)
Base.length(arr::ArrayBuffer) = arr.len

function Base.fill!(arr::ArrayBuffer, x)
    for i in 1:length(arr)
        arr[i] = x
    end
    arr
end
