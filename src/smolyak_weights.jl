function smolyak_weights(y::Array{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  interpolation_matrix = zeros(size(nodes,1),size(nodes,1))

  unique_multi_index = sort(unique(multi_index))
  unique_orders      = m_i(unique_multi_index).-1

  # Below we do the following things:

  #   Generate the polynomial terms for each order
  #   Generate the unique polynomial terms introduced at each higher order
  #   Combine the polynomial terms to construct a row of the interpolation matrix
  #   Iterate over the nodes, doing the above for steps at each iteration, to compute all rows of the interpolation matrix

  base_polynomials        = Array{Array{T,2},1}(undef,length(unique_orders))
  unique_base_polynomials = Array{Array{T,2},1}(undef,length(unique_orders))

  @inbounds for k in axes(nodes,1)

    # Construct the base polynomials

    for i in eachindex(unique_orders)
      base_polynomials[i] = chebyshev_polynomial(unique_orders[i],nodes[k,:])
    end

    # Compute the unique polynomial terms from the base polynomials

    for i = length(unique_orders):-1:2
      unique_base_polynomials[i] = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
    end
    unique_base_polynomials[1] = base_polynomials[1]

    # Construct a row of the interplation matrix

    l = 1
    @inbounds @views for j in axes(multi_index,1)
      new_polynomials = unique_base_polynomials[multi_index[j,1]][1,:]
      for i = 2:size(multi_index,2)
        new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][i,:])
      end
      m = length(new_polynomials)
      interpolation_matrix[k,l:l+m-1] = new_polynomials
      l += m
    end

  end

  weights = interpolation_matrix\y

  return weights

end

function smolyak_weights(y::Array{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  # Normalize nodes to the [-1.0 1.0] interval

  d = size(multi_index,2)
  nodes = copy(nodes)
  for i = 1:d
    nodes[:,i] = normalize_node(nodes[:,i],domain[:,i])
  end

  weights = smolyak_weights(y,nodes,multi_index)

  return weights

end

function smolyak_weights_threaded(y::Array{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  interpolation_matrix = zeros(size(nodes,1),size(nodes,1))

  unique_multi_index = sort(unique(multi_index))
  unique_orders      = m_i(unique_multi_index).-1

  # Below we do the following things:

  #   Generate the polynomial terms for each order
  #   Generate the unique polynomial terms introduced at each higher order
  #   Combine the polynomial terms to construct a row of the interpolation matrix
  #   Iterate over the nodes, doing the above for steps at each iteration, to compute all rows of the interpolation matrix

  @inbounds @sync @qthreads for k in axes(nodes,1)

    base_polynomials        = Array{Array{T,2},1}(undef,length(unique_orders))
    unique_base_polynomials = Array{Array{T,2},1}(undef,length(unique_orders))

    # Construct the base polynomials

    for i in eachindex(unique_orders)
      base_polynomials[i] = chebyshev_polynomial(unique_orders[i],nodes[k,:])
    end

    # Compute the unique polynomial terms from the base polynomials

    for i = length(unique_orders):-1:2
      unique_base_polynomials[i] = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
    end
    unique_base_polynomials[1] = base_polynomials[1]

    # Construct a row of the interplation matrix

    l = 1
    @inbounds @views for j in axes(multi_index,1)
      new_polynomials = unique_base_polynomials[multi_index[j,1]][1,:]
      for i = 2:size(multi_index,2)
        new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][i,:])
      end
      m = length(new_polynomials)
      interpolation_matrix[k,l:l+m-1] = new_polynomials
      l += m
    end

  end

  weights = interpolation_matrix\y

  return weights

end

function smolyak_weights_threaded(y::Array{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  # Normalize nodes to the [-1.0 1.0] interval

  d = size(multi_index,2)
  nodes = copy(nodes)
  for i = 1:d
    nodes[:,i] = normalize_node(nodes[:,i],domain[:,i])
  end

  weights = smolyak_weights_threaded(y,nodes,multi_index)

  return weights

end

function smolyak_inverse_interpolation_matrix(nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  interpolation_matrix = zeros(size(nodes,1),size(nodes,1))

  unique_multi_index = sort(unique(multi_index))
  unique_orders      = m_i(unique_multi_index).-1

  # Below we do the following things:

  #   Generate the polynomial terms for each order
  #   Generate the unique polynomial terms introduced at each higher order
  #   Combine the polynomial terms to construct a row of the interpolation matrix
  #   Iterate over the nodes, doing the above for steps at each iteration, to compute all rows of the interpolation matrix

  base_polynomials        = Array{Array{T,2},1}(undef,length(unique_orders))
  unique_base_polynomials = Array{Array{T,2},1}(undef,length(unique_orders))

  @inbounds for k in axes(nodes,1)

    # Construct the base polynomials

    for i in eachindex(unique_orders)
      base_polynomials[i] = chebyshev_polynomial(unique_orders[i],nodes[k,:])
    end

    # Compute the unique polynomial terms from the base polynomials

    for i = length(unique_orders):-1:2
      unique_base_polynomials[i] = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
    end
    unique_base_polynomials[1] = base_polynomials[1]

    # Construct the first row of the interplation matrix

    l = 1
    @inbounds @views for j in axes(multi_index,1)
      new_polynomials = unique_base_polynomials[multi_index[j,1]][1,:]
      for i = 2:size(multi_index,2)
        new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][i,:])
      end
      m = length(new_polynomials)
      interpolation_matrix[k,l:l+m-1] = new_polynomials
      l += m
    end

  end

  inverse_interpolation_matrix = inv(interpolation_matrix)

  return inverse_interpolation_matrix

end

function smolyak_inverse_interpolation_matrix(nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  # Normalize nodes to the [-1.0 1.0] interval

  d = size(multi_index,2)
  nodes = copy(nodes)
  for i = 1:d
    nodes[:,i] = normalize_node(nodes[:,i],domain[:,i])
  end

  inverse_interpolation_matrix = smolyak_inverse_interpolation_matrix(nodes,multi_index)

  return inverse_interpolation_matrix

end

function smolyak_inverse_interpolation_matrix_threaded(nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  interpolation_matrix = zeros(size(nodes,1),size(nodes,1))

  unique_multi_index = sort(unique(multi_index))
  unique_orders      = m_i(unique_multi_index).-1

  # Below we do the following things:

  #   Generate the polynomial terms for each order
  #   Generate the unique polynomial terms introduced at each higher order
  #   Combine the polynomial terms to construct a row of the interpolation matrix
  #   Iterate over the nodes, doing the above for steps at each iteration, to compute all rows of the interpolation matrix

  @inbounds @sync @qthreads for k in axes(nodes,1)

    base_polynomials        = Array{Array{T,2},1}(undef,length(unique_orders))
    unique_base_polynomials = Array{Array{T,2},1}(undef,length(unique_orders))

    # Construct the base polynomials

    for i in eachindex(unique_orders)
      base_polynomials[i] = chebyshev_polynomial(unique_orders[i],nodes[k,:])
    end

    # Compute the unique polynomial terms from the base polynomials

    for i = length(unique_orders):-1:2
      unique_base_polynomials[i] = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
    end
    unique_base_polynomials[1] = base_polynomials[1]

    # Construct the first row of the interplation matrix

    l = 1
    @inbounds @views for j in axes(multi_index,1)
      new_polynomials = unique_base_polynomials[multi_index[j,1]][1,:]
      for i = 2:size(multi_index,2)
        new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][i,:])
      end
      m = length(new_polynomials)
      interpolation_matrix[k,l:l+m-1] = new_polynomials
      l += m
    end

  end

  inverse_interpolation_matrix = inv(interpolation_matrix)

  return inverse_interpolation_matrix

end

function smolyak_inverse_interpolation_matrix_threaded(nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  # Normalize nodes to the [-1.0 1.0] interval

  d = size(multi_index,2)
  nodes = copy(nodes)
  for i = 1:d
    nodes[:,i] = normalize_node(nodes[:,i],domain[:,i])
  end

  inverse_interpolation_matrix = smolyak_inverse_interpolation_matrix_threaded(nodes,multi_index)

  return inverse_interpolation_matrix

end

function smolyak_weights(y::Array{T,1},inverse_interpolation_matrix::Array{T,2}) where {T<:AbstractFloat}

  weights = inverse_interpolation_matrix*y

  return weights

end