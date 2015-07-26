type ArrayBuffer{T} <: AbstractVector{T}
    array::Vector{T}
    io::IOStream
    function ArrayBuffer(len)
        if len == 0
            array = Vector{T}(len)
            a = new(array)
        else
            dirname = "."
            path, f = mktemp(dirname)
            array = Mmap.mmap(f, Vector{T}, (len,), shared=false)
            a = new(array, f)
            finalizer(a, x -> begin
                close(x.io)
                rm(path)
            end)
        end
        return a
    end
end

function Base.getindex(a::ArrayBuffer, i::Integer)
    return a.array[i]
end

function Base.setindex!{T}(a::ArrayBuffer{T}, x::T, i::Integer)
    a.array[i] = x
end

function Base.length(a::ArrayBuffer)
    length(a.array)
end

function Base.size(a::ArrayBuffer)
    size(a.array)
end

function init!(a::ArrayBuffer)
    fill!(a.array, 0)
end
