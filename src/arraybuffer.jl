type Block{T,a}
    # immutable
    size::Int
    used::Int
    # mutable
    offset::Int
    io::IO
    isdirtybuffer::Bool
    ownbuffer::Bool
    buffer::Vector{T}
    dirtybit::BitVector
end

function Block{T}(::Type{T}, size, used, buffersize, io)
    @assert buffersize ≤ size
    @assert size % buffersize == 0
    buffer = zeros(T, buffersize)
    dirtybit = falses(0)  # not used
    truncate(io, size * sizeof(T))
    Block{T,false}(size, used, 0, io, false, true, buffer, dirtybit)
end

function Block{T}(::Type{T}, used, buffer, io)
    @assert used ≤ length(buffer)
    size = length(buffer)
    dirtybit = falses(size)
    Block{T,true}(size, used, 0, io, false, false, buffer, dirtybit)
end

function init!{T,a}(block::Block{T,a})
    block.offset = 0
    seekstart(block.io)
    truncate(block.io, 0)
    if !a
        truncate(block.io, block.size * sizeof(T))
    end
    block.isdirtybuffer = false
    fill!(block.buffer, zero(T))
    fill!(block.dirtybit, false)
end

function Base.getindex{T}(block::Block{T,false}, i::Integer)
    buffersize = length(block.buffer)
    j = i - block.offset
    if 1 ≤ j ≤ buffersize
        @inbounds return block.buffer[j]
    end
    if block.isdirtybuffer
        seek(block.io, block.offset * sizeof(T))
        write(block.io, block.buffer)
    end
    block.offset = div(i - 1, buffersize) * buffersize
    seek(block.io, block.offset * sizeof(T))
    read!(block.io, block.buffer)
    block.isdirtybuffer = false
    j = i - block.offset
    return block.buffer[j]
end

function Base.setindex!{T}(block::Block{T,false}, x::T, i::Integer)
    block[i]  # load buffer
    j = i - block.offset
    block.buffer[j] = x
    #block.dirtybit[j] = true
    block.isdirtybuffer = true
    return x
end

function Base.getindex{T}(block::Block{T,true}, i::Integer)
    @assert block.ownbuffer
    return block.buffer[i]
end

function Base.setindex!{T}(block::Block{T,true}, x::T, i::Integer)
    if block.ownbuffer
        block.buffer[i] = x
    else
        seekend(block.io)
        write(block.io, i)
        write(block.io, x)
    end
    return x
end

#global ndumps = 0
function dump!{T}(block::Block{T,true})
    #global ndumps
    #@show ndumps += 1
    @assert block.ownbuffer
    for i in 1:block.size
        write(block.io, i)
        write(block.io, block.buffer[i])
    end
end

function load!{T}(block::Block{T,true})
    @assert block.ownbuffer
    seekstart(block.io)
    while !eof(block.io)
        idx = read(block.io, Int)
        elm = read(block.io, T)
        block.buffer[idx] = elm
    end
end

#=
function Base.getindex{T}(block::Block{T}, i::Integer)
    # invariant:
    #   i (absolute index in block) = j (relative index in buffer)
    #                               + block.offset (buffer offset in block)
    buffersize = length(block.buffer)
    j = i - block.offset
    if 1 ≤ j ≤ buffersize
        @inbounds return block.buffer[j]
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
=#

# parameter selectors
select_n_blocks(len) = len == 0 ? 0 : len ≤ 1024 ? 1 : 16
#select_blocksize(len, n_blocks) = n_blocks == 0 ? 64 : max(64, div(len - 1, n_blocks) + 1)
function select_blocksize(len, n_blocks)
    if n_blocks == 0
        return 16
    end
    n = div(len - 1, n_blocks) + 1
    r = n % 16
    n + 16 - r
end
#select_blocksize(len, n_blocks) = n_blocks == 0 ? 64 : 1024 * 1024
select_buffersize(blocksize) = blocksize ≤ 1024 ? blocksize : (div(blocksize - 1, 16) + 1)

type ArrayBuffer{T,a} <: AbstractVector{T}
    len::Int
    dirname::String
    blocksize::Int
    bufferowner::Int
    blocks::Vector{Block{T,a}}
    # space: buffersize * n_blocks * sizeof(T) bytes
    function ArrayBuffer{T}(::Type{T}, len::Integer;
        parent::String=pwd(),
        n_blocks=select_n_blocks(len),
        blocksize=select_blocksize(len, n_blocks),
        buffersize=select_buffersize(blocksize))
        #@show len, n_blocks, blocksize, buffersize
        @assert 0 ≤ n_blocks
        @assert 1 ≤ buffersize ≤ blocksize
        @assert 0 ≤ len ≤ blocksize * n_blocks
        dirname = mktempdir(parent)
        arr = new(len, dirname, blocksize, a ? 1 : 0)
        finalizer(arr, x -> begin
            for block in x.blocks
                if isopen(block.io)
                    close(block.io)
                end
            end
            rm(dirname, recursive=true)
        end)
        arr.blocks = Array(Block{T,a}, n_blocks)
        if a
            buffer = Array(T, blocksize)
        end
        l = len
        for i in 1:n_blocks
            _, io = mktemp(dirname)
            used = l ≥ blocksize ? blocksize : l
            l -= used
            if a
                arr.blocks[i] = Block(T, used, buffer, io)
                if i == arr.bufferowner
                    arr.blocks[i].ownbuffer = true
                end
            else
                arr.blocks[i] = Block(T, blocksize, used, buffersize, io)
            end
        end
        return arr
    end
end

function Base.call{T}(::Type{ArrayBuffer{T}}, len::Int; kwargs...)
    ArrayBuffer{T,false}(T, len; kwargs...)
end
function Base.call{T,a}(::Type{ArrayBuffer{T,a}}, len::Int; kwargs...)
    ArrayBuffer{T,a}(T, len; kwargs...)
end
#function Base.call{T}(::Type{ArrayBuffer{T,true}}, len::Int; kwargs...)
#    ArrayBuffer{T,true}(T, len; kwargs...)
#end

function move_buffer_owner!{T}(arr::ArrayBuffer{T,true}, j)
    if j == arr.bufferowner
        return
    end
    old = arr.blocks[arr.bufferowner]
    dump!(old)
    old.ownbuffer = false
    arr.bufferowner = j
    new = arr.blocks[arr.bufferowner]
    new.ownbuffer = true
    load!(new)
end

function Base.getindex{T,a}(arr::ArrayBuffer{T,a}, i::Integer)
    j, k = divrem(i - 1, arr.blocksize)
    if a
        move_buffer_owner!(arr, j + 1)
    end
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

function init!(arr::ArrayBuffer)
    for block in arr.blocks
        init!(block)
    end
end
