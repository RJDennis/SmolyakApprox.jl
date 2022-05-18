function scale_nodes(nodes::Array{R,1},domain::Array{T,1}) where {T<:AbstractFloat,R<:Number}

  nodes = copy(nodes)
  @inbounds for i in eachindex(nodes)
    nodes[i] = domain[2] + (1.0+nodes[i])*(domain[1]-domain[2])/2
  end

  return nodes

end

function scale_nodes(nodes::Array{R,2},domain::Array{T,2}) where {T<:AbstractFloat,R<:Number}

  nodes = copy(nodes)
  @inbounds for i in axes(nodes,1)
    @inbounds for j in axes(nodes,2)
      nodes[i,j] = domain[2,j] + (1.0+nodes[i,j])*(domain[1,j]-domain[2,j])/2
    end
  end

  return nodes

end

function normalize_node(node::R,domain::Array{T,1}) where {T<:AbstractFloat,R<:Number}

  if domain[1] == domain[2]
    norm_node = zero(T)
    return norm_node
  else
    norm_node = 2.0*(node-domain[2])/(domain[1]-domain[2])-1.0
    return norm_node
  end

end

function normalize_node(node::Array{R,1},domain::Array{T,1}) where {T<:AbstractFloat,R<:Number}

  norm_nodes = similar(node)
  @inbounds for i in eachindex(node)
    norm_nodes[i] = normalize_node(node[i],domain)
  end

  return norm_nodes

end