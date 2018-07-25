__precompile__()
module MyStats

export hist, hist_indices, histND_indices, min_max, min_max_ind, condmean, Bins, dbkhist_indices, dbkBins, threaded_sum, tmean, tstd, histstd

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
    hint = length(field)÷nbins 
    @inbounds for i in linearindices(indices)
        indices[i] = Vector{Int}()
        sizehint!(indices[i],hint)
    end
    
    @inbounds for i in linearindices(field)
        n = trunc(Int, ((field[i] - min)/dx - eps())*nbins) + 1
        (1 <= n <= nbins) && push!(indices[n],i)
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
        all(1 .<= n .<= nbins) && push!(indices[n...],i)
    end


    return indices
end

histND_indices(field::NTuple{N,AbstractArray},nbins::NTuple{N,Integer}) where {N} = histND_indices(field,min_max(field)...,nbins)

histND_indices(field::NTuple{3,AbstractArray}) = histND_indices(field,(5,5,5))

histND_indices(field::NTuple{2,AbstractArray}) = histND_indices(field,(12,12))

function min_max(f::NTuple{N,AbstractArray}) where N
    aux = min_max.(f)
    return (getindex.(aux,1), getindex.(aux,2))
end

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

function min_max_ind(f::AbstractArray)
    nt::Int = Threads.nthreads()
    const imin = zeros(Int,nt)
    const imax = zeros(Int,nt)

    Threads.@threads for i in 1:nt
        get_min_max_ind(f,imin,imax,i)
    end
    
    return compare_ind(f,imin,Base.:<), compare_ind(f,imax,Base.:>)
end

function get_min_max_ind(f,pmi,pma,i)
    nt::Int = Threads.nthreads()

    r = divide_range(length(f), nt)[i]

    minv = r[1]
    maxv = r[1]

    @inbounds for j in r
        if f[j] > f[maxv]
            maxv = j
        elseif f[j] < f[minv]
            minv = j
        end
    end

    pmi[i] = minv
    pma[i] = maxv
    return nothing 
end

function compare_ind(v::AbstractArray,ind,f::Function)
    tv = ind[1]

    @inbounds for j in ind
        if f(v[j],v[tv])
            tv = j
        end
    end

    return tv
end

function condmean(field::AbstractArray,condindices)
    const result = zeros(length(condindices))
    const err = zeros(length(condindices))
    Threads.@threads for i in linearindices(condindices)
        @inbounds begin
            ind = condindices[i]
            l = length(ind)
            for j in ind
                result[i] += field[j]
            end
            m = result[i]/l
            result[i] = m
            for j in ind
                err[i] += (field[j] - m)^2
            end
            err[i] = sqrt(err[i]/(l*(l-1)))
        end
    end
    return result, err
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

function threaded_sum(v,func::Function=Base.identity)
    nt::Int = Threads.nthreads()
    const result = zeros(eltype(v),nt)
    ranges = divide_range(length(v),nt)
    Threads.@threads for i in 1:nt
        get_sum!(result,v,func,ranges[i],i)
    end

    return sum(result)
end

function get_sum!(result::AbstractVector,v::AbstractArray,f::Function,ran,j::Integer)
    @inbounds begin
        r = zero(eltype(v))
        @simd for i in ran
            r += f(v[i])
        end
        result[j] = r
    end
    return nothing
end

tmean(v::AbstractArray,func::Function=Base.identity) = threaded_sum(v,func)/length(v)

function tstd(v::AbstractArray,func::Function=Base.identity,m::Number=tmean(v,func))
    f = x->((func(x)-m)^2)
    return sqrt(threaded_sum(v,f)/(length(v)-1))
end

tstd(v::AbstractArray,m::Number) = tstd(v,Base.identity,m)

function histstd(field::AbstractArray, min::Real, max::Real, nbins::Integer=100, pdf::Bool=true)
    dx = max - min
    indices = Vector{Float64}(nbins)
    @inbounds for i in linearindices(field)
        n = trunc(Int, ((field[i] - min)/dx - eps())*nbins) + 1
        (1 <= n <= nbins) && (indices[n] += 1.)
    end
    
    x = Bins(min,max,nbins)
    if pdf
        area = sum(indices)*step(x)
        indices ./= area
    end
    
    return [x indices]
end

histstd(field::AbstractArray,m::Real,stand::Real,ntimes::Real,nbins::Integer=100,pdf::Bool=true) = histstd(field,m-ntimes*stand,m+ntimes*stand,nbins,pdf)

function histstd(field::AbstractArray,ntimes::Real,nbins::Integer=100,pdf::Bool=true)
    m = tmean(field)
    stand = tstd(field,m)
    histstd(field,m-ntimes*stand,m+ntimes*stand,nbins,pdf)
end

end # module
