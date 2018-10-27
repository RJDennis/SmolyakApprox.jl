function smolyak_pl_weights(y::AbstractArray{T,1},nodes::Array{T,2},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  interpolation_matrix = zeros(size(nodes,1),size(nodes,1))
  for l = 1:size(nodes,1)
    k = 1
    x = nodes[l,:]
    for i = 1:size(multi_index,1)
      m_node_number = m_i(multi_index[i,:])
      if prod(m_node_number) == 1
        interpolation_matrix[l,k] = 1.0
        k += 1
      else
        extra_nodes = 1
        for j = 1:length(m_node_number)
          if m_node_number[j] > 1
            extra_nodes *= m_node_number[j] - m_i(multi_index[i,j]-1)
          end
        end
        for h = 1:extra_nodes
          a = 1.0
          for j = 1:length(m_node_number)
            if m_node_number[j] > 1
              if abs(x[j] - nodes[k,j]) > 2/(m_node_number[j]-1)
                a *= 0.0
              else
                a *= 1.0 - ((m_node_number[j]-1)/2)*abs(x[j]-nodes[k,j])
              end
            end
          end
          interpolation_matrix[l,k] = a
          k += 1
        end
      end
    end
  end

  weights = interpolation_matrix\y

  return weights

end

function smolyak_pl_weights(y::AbstractArray{T,1},nodes::Array{T,2},multi_index::Array{S,2},domain::Array{T,2}) where {T<:AbstractFloat,S<:Integer}

  # Normalize nodes to the [-1.0 1.0] interval

  nodes = copy(nodes)
  for i = 1:size(nodes,1)
    for j = 1:size(domain,2)
      if domain[1,j] == domain[2,j]
        nodes[i,j] = (domain[1,j]+domain[2,j])/2
      else
        nodes[i,j] = 2*(nodes[i,j]-domain[2,j])/(domain[1,j]-domain[2,j])-one(T)
      end
    end
  end

  weights = smolyak_pl_weights(y,nodes,multi_index)

  return weights

end

function smolyak_pl_evaluate(weights::Array{T,1},point::Array{T,1},nodes::Array{T,2},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  basis = zeros(size(nodes,1))
  k = 1
  for i = 1:size(multi_index,1)
    m_node_number = m_i(multi_index[i,:])
    if prod(m_node_number) == 1
      basis[k] = 1.0
      k += 1
    else
      extra_nodes = 1
      for j = 1:length(m_node_number)
        if m_node_number[j] > 1
          extra_nodes *= m_node_number[j] - m_i(multi_index[i,j]-1)
        end
      end
      for h = 1:extra_nodes
        a = 1.0
        for j = 1:length(m_node_number)
          if m_node_number[j] > 1
            if abs(point[j] - nodes[k,j]) > 2/(m_node_number[j]-1)
              a *= 0.0
            else
              a *= 1.0 - ((m_node_number[j]-1)/2)*abs(point[j]-nodes[k,j])
            end
          end
        end
        basis[k] = a
        k += 1
      end
    end
  end

  estimate = zero(T)
  for i = 1:length(basis)
    estimate += basis[i]*weights[i]
  end

  return estimate

end

function smolyak_pl_evaluate(weights::Array{T,1},point::Array{T,1},nodes::Array{T,2},multi_index::Array{S,2},domain::Array{T,2}) where {T<:AbstractFloat,S<:Integer}

  nodes = copy(nodes)
  for i = 1:size(nodes,1)
    for j = 1:size(domain,2)
      if domain[1,j] == domain[2,j]
        nodes[i,j] = (domain[1,j]+domain[2,j])/2
      else
        nodes[i,j] = 2*(nodes[i,j]-domain[2,j])/(domain[1,j]-domain[2,j])-one(T)
      end
    end
  end

  point = copy(point)
  for j = 1:size(domain,2)
    if domain[1,j] == domain[2,j]
      point[j] = (domain[1,j]+domain[2,j])/2
    else
      point[j] = 2*(point[j]-domain[2,j])/(domain[1,j]-domain[2,j])-one(T)
    end
  end

  estimate = smolyak_pl_evaluate(weights,point,nodes,multi_index)

  return estimate

end
