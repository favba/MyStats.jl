__precompile__()
module MyStats
using StatsBase: fit, Histogram

export hist

function hist(field::S,nbins::Int=250,pdf::Bool=true) where S<:AbstractArray
  array = @view(field[:])
  h = fit(Histogram,array,nbins=nbins,closed=:right)
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
  h = fit(Histogram,array,nbins=nbins,closed=:right)
  x = [(h.edges[1][i] + h.edges[1][i+1])/2 for i=1:(length(h.edges[1])-1)]
  if pdf
    y = h.weights/(sum(h.weights)*step(h.edges[1]))
  else
    y = h.weights
  end
  return [x y]
end

end # module
