__precompile__()
module MyStats
using StatsBase: fit, Histogram

export hist, hist_indices, histND_indices

function hist(field::S,nbins::Int=250,pdf::Bool=true) where S<:AbstractArray
    array = @view(field[:])
    h = fit(Histogram,array,nbins=nbins,closed=:left)
    x = [(h.edges[1][i] + h.edges[1][i+1])/2 for i=1:(length(h.edges[1])-1)]
    if pdf
        y = h.weights/(sum(h.weights)*step(h.edges[1]))
    else
        y = h.weights
    end
    return [x y]
end

function hist(field::S, indices::Array{T,1}, nbins::Int=250,pdf::Bool=true) where {S<:AbstractArray,T<:Integer}
    array = @view(field[indices])
    h = fit(Histogram,array,nbins=nbins,closed=:left)
    x = [(h.edges[1][i] + h.edges[1][i+1])/2 for i=1:(length(h.edges[1])-1)]
    if pdf
        y = h.weights/(sum(h.weights)*step(h.edges[1]))
    else
        y = h.weights
    end
    return [x y]
end

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

hist_indices(field::AbstractArray,nbins::Integer=30) = hist_indices(field,minimum(field),maximum(field),nbins)

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

end # module
