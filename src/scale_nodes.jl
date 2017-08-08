function scale_nodes(nodes::Array{T,1},domain::Array{T,1}) where {T<:AbstractFloat}

  for i = 1:length(nodes)
    nodes[i] = domain[2] + (1.0+nodes[i])*(domain[1]-domain[2])/2
  end

  return nodes

end

function scale_nodes(nodes::Array{T,2},domain::Array{T,2}) where {T<:AbstractFloat}

  for i = 1:size(nodes,1)
    for j = 1:size(nodes,2)
      nodes[i,j] = domain[2,j] + (1.0+nodes[i,j])*(domain[1,j]-domain[2,j])/2
    end
  end

  return nodes

end
