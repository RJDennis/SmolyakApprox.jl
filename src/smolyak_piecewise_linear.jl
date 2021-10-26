function smolyak_pl_weights(y::AbstractArray{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  interpolation_matrix = zeros(size(nodes,1),size(nodes,1))

  @inbounds for l = 1:size(nodes,1)
    k = 1
    x = nodes[l,:]
    @inbounds for i = 1:size(multi_index,1)
      m_node_number = m_i(multi_index[i,:])
      if prod(m_node_number) == 1
        interpolation_matrix[l,k] = 1.0
        k += 1
      else
        extra_nodes = 1
        @inbounds for j = 1:length(m_node_number)
          if m_node_number[j] > 1
            extra_nodes *= m_node_number[j] - m_i(multi_index[i,j]-1)
          end
        end
        for h = 1:extra_nodes
          a = 1.0
          @inbounds for j = 1:length(m_node_number)
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

function smolyak_pl_weights(y::AbstractArray{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  # Normalize nodes to the [-1.0 1.0] interval

  d = size(multi_index,2)
  nodes = copy(nodes)
  @inbounds for i = 1:d
    nodes[:,i] = normalize_node(nodes[:,i],domain[:,i])
  end

  weights = smolyak_pl_weights(y,nodes,multi_index)

  return weights

end

function smolyak_pl_weights_threaded(y::AbstractArray{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  interpolation_matrix = zeros(size(nodes,1),size(nodes,1))

  @inbounds @sync @qthreads for l = 1:size(nodes,1)
    k = 1
    x = nodes[l,:]
    @inbounds for i = 1:size(multi_index,1)
      m_node_number = m_i(multi_index[i,:])
      if prod(m_node_number) == 1
        interpolation_matrix[l,k] = 1.0
        k += 1
      else
        extra_nodes = 1
        @inbounds for j = 1:length(m_node_number)
          if m_node_number[j] > 1
            extra_nodes *= m_node_number[j] - m_i(multi_index[i,j]-1)
          end
        end
        for h = 1:extra_nodes
          a = 1.0
          @inbounds for j = 1:length(m_node_number)
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

function smolyak_pl_weights_threaded(y::AbstractArray{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  # Normalize nodes to the [-1.0 1.0] interval

  d = size(multi_index,2)
  nodes = copy(nodes)
  @inbounds for i = 1:d
    nodes[:,i] = normalize_node(nodes[:,i],domain[:,i])
  end

  weights = smolyak_pl_weights_threaded(y,nodes,multi_index)

  return weights

end

function smolyak_pl_evaluate(weights::Array{T,1},point::Array{T,1},nodes::Array{T,2},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  basis = zeros(size(nodes,1))
  k = 1
  @inbounds for i = 1:size(multi_index,1)
    m_node_number = m_i(multi_index[i,:])
    if prod(m_node_number) == 1
      basis[k] = 1.0
      k += 1
    else
      extra_nodes = 1
      @inbounds for j = 1:length(m_node_number)
        if m_node_number[j] > 1
          extra_nodes *= m_node_number[j] - m_i(multi_index[i,j]-1)
        end
      end
      @inbounds for h = 1:extra_nodes
        a = 1.0
        @inbounds for j = 1:length(m_node_number)
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
  @inbounds for i = 1:length(basis)
    estimate += basis[i]*weights[i]
  end

  return estimate

end

function smolyak_pl_evaluate(weights::Array{T,1},point::Array{T,1},nodes::Array{T,2},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  d = size(multi_index,2)
  nodes = copy(nodes)
  pont = copy(point)
  @inbounds for i = 1:d
    nodes[:,i] = normalize_node(nodes[:,i],domain[:,i])
    point[i] = normalize_node(point[i],domain[:,i])
  end

  estimate = smolyak_pl_evaluate(weights,point,nodes,multi_index)

  return estimate

end

###########################################################################

function smolyak_pl_evaluate(weights::Array{T,1},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  function goo(x::Array{T,1}) where {T <: AbstractFloat}

    return smolyak_pl_evaluate(weights,x,multi_index)

  end

  return goo

end

function smolyak_pl_evaluate(weights::Array{T,1},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  function goo(x::Array{T,1}) where {T <: AbstractFloat}

    return smolyak_pl_evaluate(weights,x,multi_index,domain)

  end

  return goo

end
