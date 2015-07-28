const MiB = 1024^2

type ArrayBuffer{T} <: AbstractVector{T}
    array::Vector{T}
    path::String
    io::IOStream
    function ArrayBuffer(len)
        if len * sizeof(T) ≤ 256MiB
            a = new(zeros(T, len))
            finalizer(a, x -> begin
                x.array = T[]
            end)
        else
            dirname = pwd()
            path, io = mktemp(dirname)
            array = Mmap.mmap(io, Vector{T}, (len,), shared=false)
            a = new(array, path, io)
            finalizer(a, x -> begin
                isopen(x.io) && close(x.io)
                isfile(x.path) && rm(x.path)
                isempty(x.array) || (x.array = T[])
            end)
        end
        return a
    end
end

@inline function Base.getindex(a::ArrayBuffer, i::Integer)
    @inbounds return a.array[i]
end

@inline function Base.setindex!(a::ArrayBuffer, x, i::Integer)
    @inbounds return a.array[i] = x
end

function Base.length(a::ArrayBuffer)
    length(a.array)
end

function Base.size(a::ArrayBuffer)
    size(a.array)
end

function destroy!{T}(a::ArrayBuffer{T})
    a.array = T[]
    isdefined(a, :io) && close(a.io)
    isdefined(a, :path) && rm(a.path)
end
