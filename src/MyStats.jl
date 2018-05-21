__precompile__()
module MyStats

export hist, hist_indices, histND_indices, min_max, condmean, Bins, dbkhist_indices, dbkBins

function hist(field::AbstractArray, min::Real, max::Real, nbins::Integer=100, pdf::Bool=true)
    dx = max - min
    indices = Vector{Float64}(nbins)
    
    @inbounds for i in linearindices(field)
        n = trunc(Int, ((field[i] - min)/dx - eps())*nbins) + 1
        indices[n] += 1.
    end
    
    x = Bins(min,max,nbins)
    if pdf
        area = sum(indices)*step(x)
        indices ./= area
    end
    
    return [x indices]
end

hist(field::AbstractArray, nbins::Integer=100, pdf::Bool=true) = hist(field, min_max(field)..., nbins, pdf)

hist(field::AbstractArray, indices::AbstractVector{<:Integer}, nbins::Integer=100, pdf::Bool=true) = hist(view(field,indices), nbins, pdf)

# return a Vector of Vectors of indices corresponding to the bin.
function hist_indices(field::AbstractArray,min::Real,max::Real,nbins::Integer=30)
    dx = max - min
    indices = Vector{Vector{Int}}(nbins)
    
    @inbounds for i in linearindices(indices)
        indices[i] = Vector{Int}()
    end

    @inbounds for i in linearindices(field)
        n = trunc(Int, ((field[i] - min)/dx - eps())*nbins) + 1
        push!(indices[n],i)
    end

    return indices
end

hist_indices(field::AbstractArray,nbins::Integer=30) = hist_indices(field,min_max(field)...,nbins)

function histND_indices(field::NTuple{N,AbstractArray},min::NTuple{N,Real},max::NTuple{N,Real},nbins::NTuple{N,Integer}) where {N}
    dx = max .- min
    indices = Array{Vector{Int},N}(nbins)
    
    @inbounds for i in linearindices(indices)
        indices[i] = Vector{Int}()
    end

    @inbounds for i in linearindices(field[1])
        n = @. trunc(Int, ((getindex(field,i) - min)/dx - eps())*nbins) + 1
        push!(indices[n...],i)
    end


    return indices
end

histND_indices(field::NTuple{3,AbstractArray},nbins::NTuple{3,Integer}=(5,5,5)) = histND_indices(field,minimum.(field),maximum.(field),nbins)

function min_max(f::AbstractArray)
    nt::Int = Threads.nthreads()
    const pmin = zeros(eltype(f),nt)
    const pmax = zeros(eltype(f),nt)

    Threads.@threads for i in 1:nt
        get_min_max(f,pmin,pmax,i)
    end
    
    return minimum(pmin), maximum(pmax)
end

function get_min_max(f,pmi,pma,i)
    nt::Int = Threads.nthreads()

    r = divide_range(length(f), nt)[i]

    minv = f[r[1]]
    maxv = f[r[1]]

    @inbounds for j in r
        if f[j] > maxv
            maxv = f[j]
        elseif f[j] < minv
            minv = f[j]
        end
    end

    pmi[i] = minv
    pma[i] = maxv
    return nothing 
end

function divide_range(l::Integer,np::Integer)
    p = l÷np

    m = p
    a = Vector{UnitRange{Int}}(np)
    a[1] = 1:m
    @inbounds for i=2:(np-1)
        a[i] = (m+1):(m+p)
        m = m+p
    end
    a[np] = (m+1):l
    return a
end

function condmean(field::AbstractArray,condindices)
    const result = zeros(length(condindices))
    Threads.@threads for i in linearindices(condindices)
        @inbounds begin
            ind = condindices[i]
            l = length(ind)
            for j in ind
                result[i] += field[j]
            end
            result[i] = result[i]/l
        end
    end
    return result
end

struct Bins <: AbstractVector{Float64}
    minv::Float64
    maxv::Float64
    n::Int
    dx::Float64
    Bins(minv,maxv,n) = new(minv,maxv,n,(maxv-minv)/n)
end


Base.length(a::Bins) = a.n
Base.size(a::Bins) = (a.n,)
Base.step(a::Bins) = a.dx
Base.IndexStyle(::Type{Bins}) = IndexLinear()

@inline function Base.getindex(a::Bins,i::Integer) 
    @boundscheck((1 <= i <= a.n) || throw(BoundsError(a,i)))
    dx = a.dx
    return (a.minv + 0.5*dx) + (i-1)*dx 
end

function dbkhist_indices(field::AbstractArray, nbins::Integer=30)

    indicesp = sortperm(vec(field))

    indices = Vector{Vector{Int}}(nbins)
  
    ranges = divide_range(length(field),nbins)
    @inbounds for i in 1:nbins
        indices[i] = indicesp[ranges[i]]
    end
    
    return indices

end

struct dbkBins <: AbstractVector{Float64}
    edges::Vector{Float64}
    n::Int
    dbkBins(edges) = new(edges,length(endges)-1)
end

Base.length(a::dbkBins) = a.n
Base.size(a::dbkBins) = (a.n,)
Base.IndexStyle(::Type{dbkBins}) = IndexLinear()

@inline function Base.getindex(a::dbkBins,i::Integer) 
    @boundscheck((1 <= i <= a.n) || throw(BoundsError(a,i)))
    ed = a.edges
    return 0.5*(ed[i] + ed[i+1])
end

end # module
