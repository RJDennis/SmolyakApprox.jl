# These functions relate to both an ansiotropic and an isotropic grid

function smolyak_grid(node_type::Function,d::S,mu::Union{S,Array{S,1}}) where {S<:Integer}

  T = typeof(1.0)

  multi_index        = generate_multi_index(d,mu)
  unique_multi_index = sort(unique(multi_index))
  unique_node_number = m_i(unique_multi_index)

  # Create base nodes to be used in the sparse grid

  base_nodes   = Array{Array{T,1}}(undef,length(unique_node_number))
  base_weights = Array{Array{T,1}}(undef,length(unique_node_number))
  for i = 1:length(unique_node_number)
    base_nodes[i], base_weights[i] = node_type(unique_node_number[i])
  end

  # Determine the unique nodes introduced at each higher level

  unique_base_nodes = Array{Array{T,1}}(undef,length(unique_node_number))
  for i = length(unique_node_number):-1:2
    temp_array = Array{T}(undef,1)
    for j = 1:length(base_nodes[i])
      if !in(base_nodes[i][j],base_nodes[i-1])
        push!(temp_array,base_nodes[i][j])
      end
    end
    unique_base_nodes[i] = temp_array[2:end]
  end
  unique_base_nodes[1] = base_nodes[1]

  # Construct the sparse grid from the unique nodes

  new_nodes = unique_base_nodes[multi_index[1,1]]  # Here new_nodes is a 1d array
  for i = 2:d
    new_nodes = combine_nodes(new_nodes,unique_base_nodes[multi_index[1,i]])  # Here new_nodes becomes a 2d array
  end

  nodes = copy(new_nodes)

  for j = 2:size(multi_index,1)
    new_nodes = unique_base_nodes[multi_index[j,1]]
    for i = 2:d
      new_nodes = combine_nodes(new_nodes,unique_base_nodes[multi_index[j,i]])
    end
    nodes = [nodes; new_nodes]
  end

  # Eventually this function should also return the weights at each node on the grid
  # so that it can be used for numerical integration.

  return nodes, multi_index

end

function smolyak_grid(node_type::Function,d::S,mu::Union{S,Array{S,1}},domain::Union{Array{T,1},Array{T,2}}) where {S<:Integer, T<:AbstractFloat}

  if size(domain,2) != d
    error("domain is inconsistent with the number of dimensions")
  end

  (nodes, multi_index) = smolyak_grid(node_type,d,mu)

  # Now scale the nodes to the desired domain

  nodes = scale_nodes(nodes,domain)

  return nodes, multi_index

  # Eventually this function should also return the weights at each node on the grid
  # so that it can be used for numerical integration.

end
