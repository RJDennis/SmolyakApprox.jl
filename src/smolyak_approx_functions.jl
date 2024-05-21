abstract type SApproximationPlan end

struct SApproxPlan{S<:Integer,T<:AbstractFloat} <: SApproximationPlan

  node_type::Symbol
  grid::Union{Array{T,1},Array{T,2}}
  multi_index::Union{Array{S,1},Array{S,2}}
  domain::Union{Array{T,1},Array{T,2}}

end

const chebyshev_gauss_lobatto = chebyshev_extrema

function clenshaw_curtis_equidistant(n::S,domain = [1.0,-1.0]) where {S<:Integer}

  # Construct the nodes on the [-1.0,1.0] interval

  if n <= 0
    error("The number of nodes must be positive.")
  end

  if n == 1
    nodes   = [0.0]
  else
    nodes    = zeros(n)
    nodes[1] = -1.0
    nodes[n] = 1.0

    for i = 2:div(n,2)
      nodes[i]       = 2*(i-1)/(n-1)-1.0
      nodes[end-i+1] = -2*(i-1)/(n-1)+1.0
    end

    if isodd(n)
      nodes[div(n+1,2)] = 0.0
    end
  end

  # Scale the nodes to the desired domain

  nodes = scale_nodes(nodes,domain)

  return nodes

end

# This function relates to the ansiotropic case

function generate_multi_index(d::S,mu::Array{S,1}) where {S<:Integer}

  nt = num_terms(mu,d)
  multi_index = Array{S,2}(undef,nt,d)
  multi_index[1,:] = ones(S,1,d)

  max_mu = maximum(mu)

  w = Tuple(repeat([max_mu+1],inner = d))
  pos = 0
  @inbounds for i = 2:(max_mu+1)^d
    candidate_index = Tuple(CartesianIndices(w)[i])
    if sum(candidate_index) <= d+max_mu && sum(candidate_index .<= mu.+1) == d
        pos += 1
      if pos > nt # handles the case where nt is under-estimated
        multi_index = [multi_index; collect(candidate_index)']
      else
        multi_index[pos,:] .= candidate_index
      end
    end
  end

  if pos < nt # handles case where nt is over-estimated
    multi_index = multi_index[1:pos,:]
  end

  return multi_index

end

# The function below relates to the isotropic case

function generate_multi_index(d::S,mu::S) where {S<:Integer}

  if d < 1
    error("d must be positive")
  end

  if mu < 0
    error("mu must be non-negative")
  end

  if d == 1
    multi_index = [i for i in 1:mu+1]
    return multi_index
  else
    multi_index_base = generate_multi_index(d-1,mu)
    N = size(multi_index_base,1)
    multi_index = zeros(S,N*(mu+1),d)
    pos = 0
    @inbounds @views for j = 1:N
      for i = 1:mu+1
        if sum(multi_index_base[j,:]) + i <= d+mu
          pos += 1
          multi_index[pos,2:d] .= multi_index_base[j,:]
          multi_index[pos,1] = i
        end
      end
    end
    return multi_index[1:pos,:]
  end
end    
  
# Following function computes the number of terms in the multi-index for the
# isotropic case (it also computes the number of terms in a complete
# polynominal based on the order and the number of dimensions.

function num_terms(order::S,d::S) where {S <: Integer}

  if d == 1
    return order+1
  else
    return div(num_terms(order,d-1)*(order+d),d)
  end

end

# The following function is a poor approximation to the number of terms in
# the multi-index for the ansiotropic case.

function num_terms(order::Array{S,1},d::S) where {S<:Integer}

  max_mu = maximum(order)
  nt = num_terms(max_mu,d) # Deliberate over-estimate of the number of terms
    
  return nt
  
end

m_i(x::Integer) = (x == 1 ? 1 : 2^(x-1) + 1)
m_i(x::Array{Integer,1}) = m_i.(x)

function combine_nodes(nodes1::Union{Array{R,1},Array{R,2}},nodes2::Array{R,1}) where {R<:Number}  # nodes1 can be a 1d or 2d array; nodes2 is a 1d array

  n1 = size(nodes1,1)
  n2 = size(nodes1,2)
  n3 = length(nodes2)

  combined_nodes = Array{R,2}(undef,n1*n3,n2+1)

  @inbounds for i = 1:n3
    combined_nodes[(i-1)*n1+1:i*n1,1:n2] = nodes1
  end
  @inbounds for i = 1:n1
    @inbounds for j = 1:n3
      combined_nodes[(j-1)*n1+i,n2+1] = nodes2[j]
    end
  end

  return combined_nodes

end

function scale_nodes(nodes::Array{R,1},domain::Array{T,1}) where {T<:AbstractFloat,R<:Number}

  nodes = copy(nodes)
  @inbounds for i in eachindex(nodes)
    nodes[i] = domain[2] + (1.0+nodes[i])*(domain[1]-domain[2])*0.5
  end

  return nodes

end

function scale_nodes(nodes::Array{R,2},domain::Array{T,2}) where {T<:AbstractFloat,R<:Number}

  nodes = copy(nodes)
  @inbounds for i in CartesianIndices(nodes)
    nodes[i] = domain[2,i[2]] + (1.0+nodes[i])*(domain[1,i[2]]-domain[2,i[2]])*0.5
  end

  return nodes

end

# These functions relate to both an ansiotropic and an isotropic grid

function smolyak_grid(node_type::Function,d::S,mu::Union{S,Array{S,1}}) where {S<:Integer}

  T = typeof(1.0)

  multi_index        = generate_multi_index(d,mu)
  unique_multi_index = sort(unique(multi_index))
  unique_node_number = m_i.(unique_multi_index)

  # Create base nodes to be used in the sparse grid

  base_nodes   = Array{Array{T,1},1}(undef,length(unique_node_number))
  for i in eachindex(unique_node_number)
    base_nodes[i] = node_type(unique_node_number[i])
  end

  # Determine the unique nodes introduced at each higher level

  unique_base_nodes = Array{Array{T,1},1}(undef,length(unique_node_number))
  unique_base_nodes[1] = base_nodes[1]
  for i = 2:length(unique_base_nodes)
    unique_base_nodes[i] = setdiff(base_nodes[i],base_nodes[i-1])
  end

  # Construct the sparse grid from the unique nodes

  nodes = Array{T,2}(undef,determine_grid_size(multi_index))
  l = 1
  @inbounds for j in axes(multi_index,1)
    new_nodes = unique_base_nodes[multi_index[j,1]] # Here new_nodes is a 1d array
    for i = 2:d
      new_nodes = combine_nodes(new_nodes,unique_base_nodes[multi_index[j,i]])  # Here new_nodes becomes a 2d array
    end
    m = size(new_nodes,1)
    nodes[l:l+m-1,:] = new_nodes
    l += m
  end

  # Eventually this function should also return the weights at each node on the grid
  # so that it can be used for numerical integration.

  if d == 1
    nodes = nodes[:]
  end

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

function determine_grid_size(mi)

  temp = similar(mi)

  for i in axes(mi,1)
    for j in axes(mi,2)
      if mi[i,j] == 1
        temp[i,j] = 1
      elseif mi[i,j] == 2
        temp[i,j] = 2^(mi[i,j]-1)
      else
        temp[i,j] = 2^(mi[i,j]-1)+1 - (2^(mi[i,j]-2)+1)
      end
    end
  end

  s = 0
  for i in axes(mi,1)
    t = 1
    for j in axes(mi,2)
      t *= temp[i,j]
    end
    s += t
  end

  return (s, size(mi,2))

end

function smolyak_plan(node_type::Function,d::S,mu::Union{S,Array{S,1}},domain::Union{Array{T,1},Array{T,2}}) where {S<:Integer, T<:AbstractFloat}

  g, mi = smolyak_grid(node_type,d,mu,domain)

  plan = SApproxPlan(Symbol(node_type),g,mi,domain)

  return plan

end

function smolyak_weights(y::Array{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}}) where {T<:AbstractFloat,S<:Integer}

  interpolation_matrix = zeros(size(nodes,1),size(nodes,1))

  unique_multi_index = sort(unique(multi_index))
  unique_orders      = m_i.(unique_multi_index) .- 1

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

function smolyak_weights(y::Array{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  # Normalize nodes to the [-1.0 1.0] interval

  d = size(multi_index,2)
  nodes = copy(nodes)
  for i = 1:d
    nodes[:,i] = normalize_node(nodes[:,i],domain[:,i])
  end

  weights = smolyak_weights(y,nodes,multi_index)

  return weights

end

function smolyak_weights_threaded(y::Array{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}}) where {T<:AbstractFloat,S<:Integer}

  interpolation_matrix = zeros(size(nodes,1),size(nodes,1))

  unique_multi_index = sort(unique(multi_index))
  unique_orders      = m_i.(unique_multi_index) .- 1

  # Below we do the following things:

  #   Generate the polynomial terms for each order
  #   Generate the unique polynomial terms introduced at each higher order
  #   Combine the polynomial terms to construct a row of the interpolation matrix
  #   Iterate over the nodes, doing the above for steps at each iteration, to compute all rows of the interpolation matrix

  @inbounds @sync Threads.@threads for k in axes(nodes,1)

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

function smolyak_weights_threaded(y::Array{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  # Normalize nodes to the [-1.0 1.0] interval

  d = size(multi_index,2)
  nodes = copy(nodes)
  for i = 1:d
    nodes[:,i] = normalize_node(nodes[:,i],domain[:,i])
  end

  weights = smolyak_weights_threaded(y,nodes,multi_index)

  return weights

end

function smolyak_weights(y::Array{T,1},inverse_interpolation_matrix::Array{T,2}) where {T<:AbstractFloat}

  weights = inverse_interpolation_matrix*y

  return weights

end

function smolyak_inverse_interpolation_matrix(nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}}) where {T<:AbstractFloat,S<:Integer}

  interpolation_matrix = zeros(size(nodes,1),size(nodes,1))

  unique_multi_index = sort(unique(multi_index))
  unique_orders      = m_i.(unique_multi_index) .- 1

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

function smolyak_inverse_interpolation_matrix(nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  # Normalize nodes to the [-1.0 1.0] interval

  d = size(multi_index,2)
  nodes = copy(nodes)
  for i = 1:d
    nodes[:,i] = normalize_node(nodes[:,i],domain[:,i])
  end

  inverse_interpolation_matrix = smolyak_inverse_interpolation_matrix(nodes,multi_index)

  return inverse_interpolation_matrix

end

function smolyak_inverse_interpolation_matrix_threaded(nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}}) where {T<:AbstractFloat,S<:Integer}

  interpolation_matrix = zeros(size(nodes,1),size(nodes,1))

  unique_multi_index = sort(unique(multi_index))
  unique_orders      = m_i.(unique_multi_index) .- 1

  # Below we do the following things:

  #   Generate the polynomial terms for each order
  #   Generate the unique polynomial terms introduced at each higher order
  #   Combine the polynomial terms to construct a row of the interpolation matrix
  #   Iterate over the nodes, doing the above for steps at each iteration, to compute all rows of the interpolation matrix

  @inbounds @sync Threads.@threads for k in axes(nodes,1)

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

function smolyak_inverse_interpolation_matrix_threaded(nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  # Normalize nodes to the [-1.0 1.0] interval

  d = size(multi_index,2)
  nodes = copy(nodes)
  for i = 1:d
    nodes[:,i] = normalize_node(nodes[:,i],domain[:,i])
  end

  inverse_interpolation_matrix = smolyak_inverse_interpolation_matrix_threaded(nodes,multi_index)

  return inverse_interpolation_matrix

end

function smolyak_pl_weights(y::AbstractArray{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}}) where {T<:AbstractFloat,S<:Integer}

  interpolation_matrix = zeros(size(nodes,1),size(nodes,1))

  @inbounds for l in axes(nodes,1)
    k = 1
    x = nodes[l,:]
    @inbounds for i in axes(multi_index,1)
      m_node_number = m_i.(multi_index[i,:])
      if prod(m_node_number) == 1
        interpolation_matrix[l,k] = 1.0
        k += 1
      else
        extra_nodes = 1
        @inbounds for j in eachindex(m_node_number)
          if m_node_number[j] > 1
            extra_nodes *= m_node_number[j] - m_i(multi_index[i,j] - 1)
          end
        end
        for h = 1:extra_nodes
          a = 1.0
          @inbounds for j in eachindex(m_node_number)
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

function smolyak_pl_weights(y::AbstractArray{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  # Normalize nodes to the [-1.0 1.0] interval

  d = size(multi_index,2)
  nodes = copy(nodes)
  @inbounds for i = 1:d
    nodes[:,i] = normalize_node(nodes[:,i],domain[:,i])
  end

  weights = smolyak_pl_weights(y,nodes,multi_index)

  return weights

end

function smolyak_pl_weights_threaded(y::AbstractArray{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}}) where {T<:AbstractFloat,S<:Integer}

  interpolation_matrix = zeros(size(nodes,1),size(nodes,1))

  @inbounds @sync Threads.@threads for l in axes(nodes,1)
    k = 1
    x = nodes[l,:]
    @inbounds for i in axes(multi_index,1)
      m_node_number = m_i.(multi_index[i,:])
      if prod(m_node_number) == 1
        interpolation_matrix[l,k] = 1.0
        k += 1
      else
        extra_nodes = 1
        @inbounds for j in eachindex(m_node_number)
          if m_node_number[j] > 1
            extra_nodes *= m_node_number[j] - m_i(multi_index[i,j] - 1)
          end
        end
        for h = 1:extra_nodes
          a = 1.0
          @inbounds for j in eachindex(m_node_number)
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

function smolyak_pl_weights_threaded(y::AbstractArray{T,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  # Normalize nodes to the [-1.0 1.0] interval

  d = size(multi_index,2)
  nodes = copy(nodes)
  @inbounds for i = 1:d
    nodes[:,i] = normalize_node(nodes[:,i],domain[:,i])
  end

  weights = smolyak_pl_weights_threaded(y,nodes,multi_index)

  return weights

end

function smolyak_polynomial(node::AbstractArray{R,1},multi_index::Union{Array{S,1},Array{S,2}}) where {R<:Number,S<:Integer}

  unique_multi_index = sort(unique(multi_index))
  unique_orders = m_i.(unique_multi_index) .- 1

  # Below we do the following things:

  #   Generate the polynomial terms for each order
  #   Generate the unique polynomial terms introduced at each higher order
  #   Combine the polynomial terms to construct the Smolyak polynomial

  # Here we construct the base polynomials

  base_polynomials = Array{Array{R,2}}(undef,length(unique_orders))
  for i in eachindex(unique_orders)
    base_polynomials[i] = chebyshev_polynomial(unique_orders[i],node)
  end

  # Compute the unique polynomial terms from the base polynomials

  unique_base_polynomials = Array{Array{R,2}}(undef,length(unique_orders))
  for i = length(unique_orders):-1:2
    unique_base_polynomials[i] = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
  end
  unique_base_polynomials[1] = base_polynomials[1]

  # Construct the first row of the interplation matrix

  n = determine_grid_size(multi_index)
  polynomial = Array{R,1}(undef,n[1])

  # Iterate over nodes, doing the above three steps at each iteration

  l = 1
  @inbounds for j in axes(multi_index,1)
    new_polynomials = unique_base_polynomials[multi_index[j,1]][1,:]
    for i = 2:size(multi_index,2)
      new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][i,:])
    end
    m = length(new_polynomials)
    polynomial[l:l+m-1] = new_polynomials
    l += m
  end

  return polynomial

end

function smolyak_polynomial(node::AbstractArray{R,1},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,R<:Number,S<:Integer}

  node = copy(node)

  if size(domain,2) != length(node)
    error("domain is inconsistent with the number of dimensions")
  end

  d = length(node)
  for i = 1:d
    node[i] = normalize_node(node[i],domain[:,i])
  end

  poly = smolyak_polynomial(node,multi_index)

  return poly

end

function smolyak_evaluate(weights::Array{T,1},node::AbstractArray{R,1},multi_index::Union{Array{S,1},Array{S,2}}) where {T<:AbstractFloat,R<:Number,S<:Integer}

  unique_multi_index = sort(unique(multi_index))
  unique_orders = m_i.(unique_multi_index) .- 1
  
  # Below we do the following things:
  
  #   Generate the polynomial terms for each order
  #   Generate the unique polynomial terms introduced at each higher order
  #   Combine the polynomial terms to construct the Smolyak polynomial
  
  # Here we construct the base polynomials
  
  base_polynomials = Array{Array{R,2}}(undef,length(unique_orders))
  for i in eachindex(unique_orders)
    base_polynomials[i] = chebyshev_polynomial(unique_orders[i],node)
  end
  
  # Compute the unique polynomial terms from the base polynomials
  
  unique_base_polynomials = Array{Array{R,2}}(undef,length(unique_orders))
  for i = length(unique_orders):-1:2
    unique_base_polynomials[i] = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
  end
  unique_base_polynomials[1] = base_polynomials[1]
  
  # Construct the first row of the interplation matrix
  
  polynomials = Array{R,1}(undef,length(weights))
  
  # Iterate over nodes, doing the above three steps at each iteration
  
  l = 1
    @inbounds for j in axes(multi_index,1)
    new_polynomials = unique_base_polynomials[multi_index[j,1]][1,:]
    for i = 2:size(multi_index,2)
      new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][i,:])
    end
    m = length(new_polynomials)
    polynomials[l:l+m-1] = new_polynomials
    l += m
  end
  
  estimate = zero(T)
  for i in eachindex(polynomials)
    estimate += polynomials[i]*weights[i]
  end
  
  return estimate
  
end
  
function smolyak_evaluate(weights::Array{T,1},node::AbstractArray{R,1},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  node = copy(node)
  
  if size(domain,2) != length(node)
    error("domain is inconsistent with the number of dimensions")
  end
  
  d = length(node)
  for i = 1:d
    node[i] = normalize_node(node[i],domain[:,i])
  end
  
  estimate = smolyak_evaluate(weights,node,multi_index)
  
  return estimate
  
end

function smolyak_evaluate(weights::Array{T,1},polynomial::Array{R,1}) where {T<:AbstractFloat,R<:Number}

  estimate = weights'polynomial
  
  return estimate

end

function smolyak_pl_evaluate(weights::Array{T,1},point::Array{R,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}}) where {T<:AbstractFloat,R<:Number,S<:Integer}

  basis = Array{R,1}(undef,size(nodes,1))
  k = 1
  @inbounds for i in axes(multi_index,1)
    m_node_number = m_i.(multi_index[i,:])
    if prod(m_node_number) == 1
      basis[k] = one(R)
      k += 1
    else
      extra_nodes = 1
      @inbounds for j in eachindex(m_node_number)
        if m_node_number[j] > 1
          extra_nodes *= m_node_number[j] - m_i(multi_index[i,j] - 1)
        end
      end
      @inbounds for h = 1:extra_nodes
        a = 1.0
        @inbounds for j in eachindex(m_node_number)
          if m_node_number[j] > 1
            if abs(point[j] - nodes[k,j]) > 2/(m_node_number[j]-1)
              a *= zero(R)
            else
              a *= one(R) - ((m_node_number[j]-1)/2)*abs(point[j]-nodes[k,j])
            end
          end
        end
        basis[k] = a
        k += 1
      end
    end
  end

  estimate = zero(R)
  @inbounds for i in eachindex(basis)
    estimate += basis[i]*weights[i]
  end

  return estimate

end

function smolyak_pl_evaluate(weights::Array{T,1},point::Array{R,1},nodes::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,R<:Number,S<:Integer}

  d = size(multi_index,2)
  nodes = copy(nodes)
  point = copy(point)
  @inbounds for i = 1:d
    nodes[:,i] = normalize_node(nodes[:,i],domain[:,i])
    point[i] = normalize_node(point[i],domain[:,i])
  end

  estimate = smolyak_pl_evaluate(weights,point,nodes,multi_index)

  return estimate

end

function smolyak_interp(y::Array{T,1},plan::P) where {T<:AbstractFloat,P<:SApproximationPlan}

  if plan.node_type == :chebyshev_extrema || plan.node_type == :chebyshev_gauss_lobatto
    weights = smolyak_weights(y,plan.grid,plan.multi_index,plan.domain)
  elseif plan.node_type == :clenshaw_curtis_equidistant
    weights = smolyak_pl_weights(y,plan.grid,plan.multi_index,plan.domain)
  end

  function interp(x::Array{R,1}) where {R<:Number}

    if plan.node_type == :chebyshev_extrema || plan.node_type == :chebyshev_gauss_lobatto
      return smolyak_evaluate(weights,x,plan.multi_index,plan.domain)
    elseif plan.node_type == :clenshaw_curtis_equidistant
      return smolyak_pl_evaluate(weights,x,plan.grid,plan.multi_index,plan.domain)
    end

  end

  return interp

end

function smolyak_interp_threaded(y::Array{T,1},plan::P) where {T<:AbstractFloat,P<:SApproximationPlan}

  if plan.node_type == :chebyshev_extrema || plan.node_type == :chebyshev_gauss_lobatto
    weights = smolyak_weights_threaded(y,plan.grid,plan.multi_index,plan.domain)
  elseif plan.node_type ==:clenshaw_curtis_equidistant
    weights = smolyak_pl_weights_threaded(y,plan.grid,plan.multi_index,plan.domain)
  end

  function interp(x::Array{R,1}) where {R<:Number}

    if plan.node_type == :chebyshev_extrema || plan.node_type == :chebyshev_gauss_lobatto
      return smolyak_evaluate(weights,x,plan.multi_index,plan.domain)
    elseif plan.node_type == :clenshaw_curtis_equidistant
      return smolyak_pl_evaluate(weights,x,plan.grid,plan.multi_index,plan.domain)
    end

  end

  return interp

end

function smolyak_derivative(weights::Array{T,1},node::Array{R,1},multi_index::Union{Array{S,1},Array{S,2}},pos::S) where {T<:AbstractFloat,R<:Number,S<:Integer}

  unique_multi_index = sort(unique(multi_index))
  unique_orders = m_i.(unique_multi_index) .- 1

  # Here we construct the base polynomials

  base_polynomials            = Array{Array{R,2},1}(undef,length(unique_orders))
  base_polynomial_derivatives = Array{Array{R,2},1}(undef,length(unique_orders))
  for i in eachindex(unique_orders)
    base_polynomials[i]            = chebyshev_polynomial(unique_orders[i],node)
    base_polynomial_derivatives[i] = chebyshev_polynomial_deriv(unique_orders[i],node)
  end

  # Compute the unique polynomial terms from the base polynomials

  unique_base_polynomials = Array{Array{R,2},1}(undef,length(unique_orders))
  unique_base_polynomial_derivatives = Array{Array{R,2},1}(undef,length(unique_orders))
  for i = length(unique_orders):-1:2
    unique_base_polynomials[i] = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
    unique_base_polynomial_derivatives[i] = base_polynomial_derivatives[i][:,size(base_polynomial_derivatives[i-1],2)+1:end]
  end
  unique_base_polynomials[1] = base_polynomials[1]
  unique_base_polynomial_derivatives[1] = base_polynomial_derivatives[1]

  # Construct the first row of the interplation matrix

  polynomials = Array{R,1}(undef,length(weights))

  # Iterate over nodes, doing the above three steps at each iteration

  l = 1
  @inbounds for j in axes(multi_index,1)
    if pos == 1
      new_polynomials = unique_base_polynomial_derivatives[multi_index[j,1]][1,:]
    else
      new_polynomials = unique_base_polynomials[multi_index[j,1]][1,:]
    end
    for i = 2:size(multi_index,2)
      if pos == i
        new_polynomials = kron(new_polynomials,unique_base_polynomial_derivatives[multi_index[j,i]][i,:])
      else
        new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][i,:])
      end
    end
    m = length(new_polynomials)
    polynomials[l:l+m-1] = new_polynomials
    l += m
  end

  evaluated_derivative = zero(T)

  for i in eachindex(polynomials)
    evaluated_derivative += polynomials[i]*weights[i]
  end

  return evaluated_derivative

end

function smolyak_derivative(weights::Array{T,1},node::Array{R,1},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}},pos::S) where {T<:AbstractFloat,R<:Number,S<:Integer}

  node = copy(node)

  if size(domain,2) != length(node)
    error("domain is inconsistent with the number of dimensions")
  end

  d = length(node)
  for i = 1:d
    node[i] = normalize_node(node[i],domain[:,i])
  end

  evaluated_derivative = smolyak_derivative(weights,node,multi_index,pos)

  return evaluated_derivative*(2.0/(domain[1,pos]-domain[2,pos]))

end

function smolyak_gradient(weights::Array{T,1},node::Array{R,1},multi_index::Union{Array{S,1},Array{S,2}}) where {T<:AbstractFloat,R<:Number,S<:Integer}

  d = length(node)
  gradient = Array{R,2}(undef,1,d)

  for i = 1:d
    gradient[i] = smolyak_derivative(weights,node,multi_index,i)
  end

  return gradient

end

function smolyak_gradient(weights::Array{T,1},node::Array{R,1},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,R<:Number,S<:Integer}

  d = length(node)
  gradient = Array{R,2}(undef,1,d)

  for i = 1:d
    gradient[i] = smolyak_derivative(weights,node,multi_index,domain,i)
  end

  return gradient

end

function smolyak_gradient(y::AbstractArray{T,1},plan::P) where {T<:AbstractFloat,P<:SApproximationPlan}

  if plan.node_type == :clenshaw_curtis_equidistant
    error("Not implemented for clenshaw_curtis_equidistant nodes")
  end

  weights = smolyak_weights(y,plan.grid,plan.multi_index,plan.domain)
  
  function smolyak_grad(x::Array{R,1}) where {R<:Number}
  
    return smolyak_gradient(weights,x,plan.multi_index,plan.domain)
  
  end
  
  return smolyak_grad
  
end
  
function smolyak_gradient_threaded(y::AbstractArray{T,1},plan::P) where {T<:AbstractFloat,P<:SApproximationPlan}
  
  if plan.node_type == :clenshaw_curtis_equidistant
    error("Not implemented for clenshaw_curtis_equidistant nodes")
  end

  weights = smolyak_weights_threaded(y,plan.grid,plan.multi_index,plan.domain)
  
  function smolyak_grad(x::Array{R,1}) where {R<:Number}
  
    return smolyak_gradient(weights,x,plan.multi_index,plan.domain)
  
  end
  
  return smolyak_grad
  
end

function smolyak_hessian(weights::Array{T,1},point::Array{R,1},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  point = copy(point)

  if size(domain,2) != length(point)
    error("domain is inconsistent with the number of dimensions")
  end
  
  d = length(point)
  for i = 1:d
    point[i] = normalize_node(point[i],domain[:,i])
  end

  unique_multi_index = sort(unique(multi_index))
  unique_orders = m_i.(unique_multi_index) .- 1

  hess = Array{T,2}(undef,d,d)

  unique_multi_index = sort(unique(multi_index))
  unique_orders = m_i.(unique_multi_index) .- 1

  # Here we construct the base polynomials

  base_polynomials                = Array{Array{R,2},1}(undef,length(unique_orders))
  base_polynomial_derivatives     = Array{Array{R,2},1}(undef,length(unique_orders))
  base_polynomial_sec_derivatives = Array{Array{R,2},1}(undef,length(unique_orders))
  for i in eachindex(unique_orders)
    base_polynomials[i]                = chebyshev_polynomial(unique_orders[i],point)
    base_polynomial_derivatives[i]     = chebyshev_polynomial_deriv(unique_orders[i],point)
    base_polynomial_sec_derivatives[i] = chebyshev_polynomial_sec_deriv(unique_orders[i],point)
  end

  # Compute the unique polynomial terms from the base polynomials

  unique_base_polynomials                = Array{Array{R,2},1}(undef,length(unique_orders))
  unique_base_polynomial_derivatives     = Array{Array{R,2},1}(undef,length(unique_orders))
  unique_base_polynomial_sec_derivatives = Array{Array{R,2},1}(undef,length(unique_orders))
  for i = length(unique_orders):-1:2
    unique_base_polynomials[i]                = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
    unique_base_polynomial_derivatives[i]     = base_polynomial_derivatives[i][:,size(base_polynomial_derivatives[i-1],2)+1:end]
    unique_base_polynomial_sec_derivatives[i] = base_polynomial_sec_derivatives[i][:,size(base_polynomial_sec_derivatives[i-1],2)+1:end]
  end
  unique_base_polynomials[1]                = base_polynomials[1]
  unique_base_polynomial_derivatives[1]     = base_polynomial_derivatives[1]
  unique_base_polynomial_sec_derivatives[1] = base_polynomial_sec_derivatives[1]

  # Construct the first row of the interplation matrix

  polynomials = Array{R,1}(undef,length(weights))

  # Iterate over nodes, doing the above three steps at each iteration

  @inbounds for c in CartesianIndices(hess)
    l = 1
    @inbounds for j in axes(multi_index,1)
      if 1 == c[1] == c[2]
        new_polynomials = unique_base_polynomial_sec_derivatives[multi_index[j,1]][1,:]
      elseif 1 == c[1] || 1 == c[2]
        new_polynomials = unique_base_polynomial_derivatives[multi_index[j,1]][1,:]
      else
        new_polynomials = unique_base_polynomials[multi_index[j,1]][1,:]
      end
      for i = 2:size(multi_index,2)
        if i == c[1] == c[2]
          new_polynomials = kron(new_polynomials,unique_base_polynomial_sec_derivatives[multi_index[j,i]][i,:])
        elseif i == c[1] || i == c[2]
          new_polynomials = kron(new_polynomials,unique_base_polynomial_derivatives[multi_index[j,i]][i,:])
        else
          new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][i,:])
        end
      end
      m = length(new_polynomials)
      polynomials[l:l+m-1] = new_polynomials
      l += m
    end

    evaluated_derivative = zero(T)

    for i in eachindex(polynomials)
      evaluated_derivative += polynomials[i]*weights[i]
    end

    hess[c] = evaluated_derivative*(2.0/(domain[1,c[1]]-domain[2,c[1]]))*(2.0/(domain[1,c[2]]-domain[2,c[2]]))

  end

  return hess

end

function smolyak_hessian(y::AbstractArray{T,1},plan::P) where {T<:AbstractFloat,P<:SApproximationPlan}

  weights = smolyak_weights(y,plan.grid,plan.multi_index,plan.domain)
  
  function smolyak_hess(x::Array{R,1}) where {R<:Number}
  
    return smolyak_hessian(weights,x,plan.multi_index,plan.domain)
  
  end
  
  return smolyak_hess
  
end
  
function smolyak_hessian_threaded(y::AbstractArray{T,1},plan::P) where {T<:AbstractFloat,P<:SApproximationPlan}
  
  weights = smolyak_weights_threaded(y,plan.grid,plan.multi_index,plan.domain)
  
  function smolyak_hess(x::Array{R,1}) where {R<:Number}
  
    return smolyak_hessian(weights,x,plan.multi_index,plan.domain)
  
  end
  
  return smolyak_hess
  
end

function integrate_cheb_polys(order::S) where {S <: Integer}

  # Integrates Chebyshev polynomials over the domain [-1,1]

  p = zeros(order+1)
    
  for i in 1:order+1
    if i == 2
      p[i] = 0.0
    else
      p[i] = ((-1)^(i-1)+1)/(1-(i-1)^2)
    end
  end

  return p

end

function smolyak_integrate(f::Function,plan::SApproxPlan,method::Symbol)

    if method == :clenshaw_curtis
        integral = smolyak_clenshaw_curtis(f,plan)
    elseif method == :gauss_chebyshev_quad
        integral = smolyak_gauss_chebyshev_quad(f,plan)
    else
        error("Integration not implemented for that method")
    end

    return integral

end

function smolyak_clenshaw_curtis(f::Function,plan::SApproxPlan)

  grid        = plan.grid
  multi_index = plan.multi_index
  domain      = plan.domain

  y = zeros(size(grid,1))
  for i in eachindex(y)
    y[i] = f(grid[i,:])
  end

  weights = smolyak_weights(y,grid,multi_index,domain)

  # Uses Clenshaw-Curtis to integrate over all dimensions

  unique_multi_index = sort(unique(multi_index))
  unique_orders = m_i.(unique_multi_index) .- 1

  # Here we construct the base polynomials

  T = eltype(grid)

  base_polynomial_integrals = Array{Array{T,1},1}(undef,length(unique_orders))
  for i in eachindex(unique_orders)
    base_polynomial_integrals[i] = integrate_cheb_polys(unique_orders[i])
  end

  # Compute the unique polynomial terms from the base polynomials

  unique_base_polynomial_integrals = Array{Array{T,1},1}(undef,length(unique_orders))
  for i = length(unique_orders):-1:2
    unique_base_polynomial_integrals[i] = base_polynomial_integrals[i][length(base_polynomial_integrals[i-1])+1:end]
  end
  unique_base_polynomial_integrals[1] = base_polynomial_integrals[1]

  # Construct the first row of the interplation matrix

  polynomials = Array{T,1}(undef,length(weights))

  # Iterate over nodes, doing the above three steps at each iteration

  l = 1
  @inbounds for j in axes(multi_index,1)
    new_polynomials = unique_base_polynomial_integrals[multi_index[j,1]][:]
    for i = 2:size(multi_index,2)
      new_polynomials = kron(new_polynomials,unique_base_polynomial_integrals[multi_index[j,i]][:])
    end
    m = length(new_polynomials)
    polynomials[l:l+m-1] = new_polynomials
    l += m
  end

  evaluated_integral = zero(T)

  for i in eachindex(polynomials)
    evaluated_integral += polynomials[i]*weights[i]
  end

  scale_factor = (domain[1,1]-domain[2,1])/2
  for i in 2:size(multi_index,2)
    scale_factor = scale_factor*(domain[1,i]-domain[2,i])/2
  end

  return evaluated_integral*scale_factor

end

function smolyak_clenshaw_curtis(f::Function,plan::SApproxPlan,pos::S) where {S<:Integer}

  # Uses Clenshaw-Curtis to integrate over all dimensions except for pos

  grid        = plan.grid
  multi_index = plan.multi_index
  domain      = plan.domain

  y = zeros(size(grid,1))
  for i in eachindex(y)
    y[i] = f(grid[i,:])
  end

  weights = smolyak_weights(y,grid,multi_index,domain)

  function smolyak_int(point::R) where {R <: Number}
    
    point = normalize_node(point,domain[:,pos])

    unique_multi_index = sort(unique(multi_index))
    unique_orders = m_i.(unique_multi_index) .- 1

    # Here we construct the base polynomials

    T = eltype(grid)

    base_polynomials          = Array{Array{T,1},1}(undef,length(unique_orders))
    base_polynomial_integrals = Array{Array{T,1},1}(undef,length(unique_orders))
    for i in eachindex(unique_orders)
      base_polynomials[i]          = chebyshev_polynomial(unique_orders[i],point)[:]
      base_polynomial_integrals[i] = integrate_cheb_polys(unique_orders[i])
    end

    # Compute the unique polynomial terms from the base polynomials

    unique_base_polynomials          = Array{Array{T,1},1}(undef,length(unique_orders))
    unique_base_polynomial_integrals = Array{Array{T,1},1}(undef,length(unique_orders))
    for i = length(unique_orders):-1:2
      unique_base_polynomials[i]          = base_polynomials[i][length(base_polynomials[i-1])+1:end]
      unique_base_polynomial_integrals[i] = base_polynomial_integrals[i][length(base_polynomial_integrals[i-1])+1:end]
    end
    unique_base_polynomials[1]          = base_polynomials[1]
    unique_base_polynomial_integrals[1] = base_polynomial_integrals[1]

    # Construct the first row of the interplation matrix

    polynomials = Array{T,1}(undef,length(weights))

    # Iterate over nodes, doing the above three steps at each iteration

    l = 1
    @inbounds for j in axes(multi_index,1)
      if pos == 1
        new_polynomials = unique_base_polynomials[multi_index[j,1]][:]
      else
        new_polynomials = unique_base_polynomial_integrals[multi_index[j,1]][:]
      end
      for i = 2:size(multi_index,2)
        if pos == i
          new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][:])
        else
          new_polynomials = kron(new_polynomials,unique_base_polynomial_integrals[multi_index[j,i]][:])
        end
      end
      m = length(new_polynomials)
      polynomials[l:l+m-1] = new_polynomials
      l += m
    end

    evaluated_integral = zero(T)

    for i in eachindex(polynomials)
      evaluated_integral += polynomials[i]*weights[i]
    end

    scale_factor = 1.0
    for i in 1:size(multi_index,2)
      if pos != i
        scale_factor = scale_factor*(domain[1,i]-domain[2,i])/2
      end
    end

    return evaluated_integral*scale_factor

  end

  return smolyak_int

end

function smolyak_gauss_chebyshev_quad(f::Function,plan::SApproxPlan)

  # Uses Gauss-Chebyshev quadrature to integrate over all dimensions
  
  grid        = plan.grid
  multi_index = plan.multi_index
  domain      = plan.domain

  iim = smolyak_inverse_interpolation_matrix(grid,multi_index,domain) 

  d = size(grid,2)

  e = zeros(1,size(grid,1))
  e[1] = π^d
  w = e*iim

  y = zeros(size(grid,1))
  for i in eachindex(y)
    integrating_weights = sqrt(1.0-normalize_node(grid[i,1],domain[:,1])^2)
    for j = 2:d
      integrating_weights *= sqrt(1.0-normalize_node(grid[i,j],domain[:,j])^2)
    end
    y[i] = f(grid[i,:])*integrating_weights
  end

  scale_factor = (domain[1,1]-domain[2,1])/2
  for i in 2:d
    scale_factor = scale_factor*(domain[1,i]-domain[2,i])/2
  end

  return (w*y)[1]*scale_factor

end

########################################################################
########################################################################

function smolyak_grid_full(node_type::Function,d::S,mu::S) where {S<:Integer}

  T = typeof(1.0)

  multi_index        = generate_multi_index(d,mu)
  unique_multi_index = sort(unique(multi_index))
  unique_node_number = m_i.(unique_multi_index)

  # Create base nodes to be used in the sparse grid

  base_nodes   = Array{Array{T,1},1}(undef,length(unique_node_number))
  base_weights = Array{Array{T,1},1}(undef,length(unique_node_number))
  @inbounds for i in eachindex(unique_node_number)
    base_nodes[i] = node_type(unique_node_number[i])
  end

  # Select the relevant polynomials from the multi index

  ii = (sum(multi_index,dims=2) .>= max(d,mu+1)).*(sum(multi_index,dims=2) .<= d+mu)
  multi_index_full = zeros(S,sum(ii),size(multi_index,2))
  j = 1
  @inbounds for i in axes(multi_index,1)
    if ii[i] == true
      multi_index_full[j,:] = multi_index[i,:]
      j += 1
    end
  end

  # Construct the sparse grid from the nodes

  nodes = Array{T,2}(undef,determine_grid_size_full(multi_index_full))
  l = 1
  for j in axes(multi_index_full,1)
    new_nodes = base_nodes[multi_index_full[j,1]]  # Here new_nodes is a 1d array
    for i = 2:d
      new_nodes = combine_nodes(new_nodes,base_nodes[multi_index_full[j,i]])  # Here new_nodes becomes a 2d array
    end
    m = size(new_nodes,1)
    nodes[l:l+m-1,:] = new_nodes
    l += m
  end

  return nodes, multi_index_full

end

function smolyak_grid_full(node_type::Function,d::S,mu::S,domain::Union{Array{T,1},Array{T,2}}) where {S<:Integer, T<:AbstractFloat}

  if size(domain,2) != d
    error("domain is inconsistent with the number of dimensions")
  end

  nodes, multi_index_full = smolyak_grid_full(node_type,d,mu)

  nodes = scale_nodes(nodes,domain)

  return nodes, multi_index_full

end

function determine_grid_size_full(mi)

  temp = similar(mi)

  @inbounds for i in axes(mi,1)
    @inbounds for j in axes(mi,2)
      if mi[i,j] == 1
        temp[i,j] = 1
      else
        temp[i,j] = 2^(mi[i,j]-1)+1
      end
    end
  end

  s = 0
  @inbounds for i in axes(mi,1)
    t = 1
    @inbounds for j in axes(mi,2)
      t *= temp[i,j]
    end
    s += t
  end

  return (s, size(mi,2))

end

function master_index(multi_index::Array{S,2}) where {S<:Integer}

  temp_ind   = similar(multi_index)
  master_ind = zeros(S,size(multi_index,1),2)

  @inbounds for i in eachindex(multi_index)
    if multi_index[i] == 1
      temp_ind[i] = 1
    else
      temp_ind[i] = m_i(multi_index[i])
    end
  end

  master_ind[1,1] = 1
  master_ind[1,2] = prod(temp_ind[1,:])

  @inbounds for i = 2:size(master_ind,1)
    master_ind[i,1] = master_ind[i-1,1] + master_ind[i-1,2]
    master_ind[i,2] = prod(temp_ind[i,:])
  end

  return master_ind

end

function cheb_poly(order::S,x::R) where {S<:Integer,R<:Number}

  p  = one(R)
  p1 = zero(R)
  p2 = zero(R)

  for i = 2:order+1
    if i == 2
      p1, p = p, x
    else
      p2, p1 = p1, p
        p  = 2*x*p1-p2
    end
  end

  return p

end

function prod_cjs(max_grid::Union{Array{T,1},Array{T,2}},min_grid::Union{Array{T,1},Array{T,2}},poly_grid::Array{T,2}) where {T<:AbstractFloat}

  cjs = ones(size(poly_grid))

  @inbounds for i in axes(poly_grid,1)
    @inbounds for j in axes(poly_grid,2)
      if poly_grid[i,j] == max_grid[j] || poly_grid[i,j] == min_grid[j]
        cjs[i,j] *= 2
      end
    end
  end

  return prod(cjs,dims=2)

end

function compute_scale_factor(multi_index::Array{S,1}) where {S<:Integer}

  scale_factor = 1.0

  @inbounds for i in eachindex(multi_index)
    if multi_index[i] > 1
      scale_factor *= 2.0/(m_i(multi_index[i]) - 1)
    end
  end

  return scale_factor

end

function smolyak_weights_full(y_f::Array{T,1},grid::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2}) where {S<:Integer, T<:AbstractFloat}

  mi = sum(multi_index,dims=2)
  d  = size(multi_index,2)
  mu = maximum(mi)-d

  max_grid = maximum(grid,dims=1)
  min_grid = minimum(grid,dims=1)

  weights = Array{Array{T,1},1}(undef,size(multi_index,1))
  g_ind = master_index(multi_index)

  @inbounds for i in axes(g_ind,1) # This loops over the number of polynomials
    ws = zeros(g_ind[i,2])
    poly_grid = grid[g_ind[i,1]:g_ind[i,1]+g_ind[i,2]-1,:]
    poly_y    = y_f[g_ind[i,1]:g_ind[i,1]+g_ind[i,2]-1]
      if size(grid,1) == 1 # This is to accommodate the mu = 0 case
        cjs_prod = ones(size(poly_grid,1))
      else
      cjs_prod = prod_cjs(max_grid,min_grid,poly_grid)
    end
    @inbounds for l = 1:g_ind[i,2] # This loops over the weights
      ll = CartesianIndices(Tuple(m_i.(multi_index[i,:])))[l]
      @inbounds for j = 1:g_ind[i,2] # This loops over the nodes
        rhs_term = cheb_poly(ll[1]-1,poly_grid[j,1])*poly_y[j]
        @inbounds for k = 2:size(poly_grid,2) # This loops over the polynomial terms in the product
          rhs_term *= cheb_poly(ll[k]-1,poly_grid[j,k])
        end
        ws[l] += rhs_term/cjs_prod[j]
      end
      scale_factor = compute_scale_factor(multi_index[i,:])
      ws[l] = scale_factor*(1/cjs_prod[l])*ws[l]
    end
    weights[i] = (-1)^(d+mu-mi[i])*factorial(d-1)/(factorial(d+mu-mi[i])*factorial(-1-mu+mi[i]))*ws
  end

  return weights

end

function smolyak_weights_full(y_f::Array{T,1},grid::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2},domain::Array{T,2}) where {S<:Integer, T<:AbstractFloat}

  d = size(multi_index,2)
  grid = copy(grid)
  for i = 1:d
    grid[:,i] = normalize_node(grid[:,i],domain[:,i])
  end

  weights = smolyak_weights_full(y_f,grid,multi_index)

  return weights

end

function smolyak_evaluate_full(weights::Array{Array{T,1},1},point::Array{R,1},multi_index::Array{S,2}) where {S<:Integer,R<:Number,T<:AbstractFloat}

  d = size(multi_index,2)

  evaluated_polynomials = zeros(size(multi_index,1))

  @inbounds for i in axes(multi_index,1) # This loops over the number of polynomials
    @inbounds for l in eachindex(weights[i])
      ll = CartesianIndices(Tuple(m_i.(multi_index[i:i,:])))[l]
      temp = weights[i][l]*cheb_poly(ll[1]-1,point[1])
      @inbounds for k = 2:d
        temp *= cheb_poly(ll[k]-1,point[k])
      end
      evaluated_polynomials[i] += temp
    end
    #evaluated_polynomials[i] *= (-1)^(d+mu-mi[i])*factorial(d-1)/(factorial(d+mu-mi[i])*factorial(-1-mu+mi[i]))
  end

  return sum(evaluated_polynomials)

end

function smolyak_evaluate_full(weights::Array{Array{T,1},1},point::Array{R,1},multi_index::Array{S,2},domain::Array{T,2}) where {S<:Integer,R<:Number,T<:AbstractFloat}

  d = size(multi_index,2)
  point = copy(point)
  for i = 1:d
    point[i] = normalize_node(point[i],domain[:,i])
  end

  estimate = smolyak_evaluate_full(weights,point,multi_index)

  return estimate

end

function deriv_cheb_poly(order::S,x::R) where {S<:Integer,R<:Number}

  p0 = one(R)
  p1 = zero(R)
  p2 = zero(R)
  pd0 = zero(R)
  pd1 = zero(R)
  pd2 = zero(R)

  for i = 2:order+1
    if i == 2
      p1, p0 = p0, x
      pd1, pd0 = pd0, one(R)
    else
      p2, p1 = p1, p0
      pd2, pd1 = pd1, pd0
        p0  = 2*x*p1-p2
      pd0 = 2*p1+2*x*pd1-pd2
    end
  end

  return pd0

end

function smolyak_derivative_full(weights::Array{Array{T,1},1},point::Array{R,1},multi_index::Array{S,2},pos::S) where {S<:Integer,R<:Number,T<:AbstractFloat}

  mi = sum(multi_index,dims=2)
  d  = size(multi_index,2)
  mu = maximum(mi)-d

  evaluated_polynomials = zeros(size(multi_index,1))

  for i in axes(multi_index,1) # This loops over the number of polynomials
    for l in eachindex(weights[i])
      ll = CartesianIndices(Tuple(m_i.(multi_index[i:i,:])))[l]
      if pos == 1
        temp = weights[i][l]*deriv_cheb_poly(ll[1]-1,point[1])
      else
        temp = weights[i][l]*cheb_poly(ll[1]-1,point[1])
      end
      for k = 2:d
        if k == pos
          temp *= deriv_cheb_poly(ll[k]-1,point[k])
        else
          temp *= cheb_poly(ll[k]-1,point[k])
        end
      end
      evaluated_polynomials[i] += temp
    end
    #evaluated_polynomials[i] *= (-1)^(d+mu-mi[i])*factorial(d-1)/(factorial(d+mu-mi[i])*factorial(-1-mu+mi[i]))
  end

  evaluated_derivative = sum(evaluated_polynomials)

  return evaluated_derivative

end

function smolyak_derivative_full(weights::Array{Array{T,1},1},point::Array{R,1},multi_index::Array{S,2},domain::Array{T,2},pos::S) where {S<:Integer,R<:Number,T<:AbstractFloat}

  d = size(multi_index,2)
  point = copy(point)
  for i = 1:d
    point[i] = normalize_node(point[i],domain[:,i])
  end

  evaluated_derivative = smolyak_derivative_full(weights,point,multi_index,pos)

  return evaluated_derivatives

end

function smolyak_gradient_full(weights::Array{Array{T,1},1},point::Array{R,1},multi_index::Array{S,2}) where {S<:Integer,R<:Number,T<:AbstractFloat}

  d = size(multi_index,2)

  gradient = Array{R,2}(undef,1,d)
  for i = 1:d
    gradient[i] = smolyak_derivative_full(weights,point,multi_index,i)
  end

  return gradient

end

function smolyak_gradient_full(weights::Array{Array{T,1},1},point::Array{R,1},multi_index::Array{S,2},domain::Array{T,2}) where {S<:Integer,R<:Number,T<:AbstractFloat}

  d = size(multi_index,2)

  gradient = Array{R,2}(undef,1,d)
  for i = 1:d
    gradient[i] = smolyak_derivative_full(weights,point,multi_index,domain,i)
  end

  return gradient

end

###########################################################################

function smolyak_evaluate_full(weights::Array{T,1},multi_index_full::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  function goo(x::Array{R,1}) where {R<:Number}

    return smolyak_evaluate_full(weights,x,multi_index_full)

  end

  return goo

end

function smolyak_evaluate_full(weights::Array{T,1},multi_index_full::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  function goo(x::Array{R,1}) where {R<:Number}

    return smolyak_evaluate_full(weights,x,multi_index_full,domain)

  end

  return goo

end

function smolyak_derivative_full(weights::Array{T,1},multi_index_full::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  function goo(x::Array{R,1}) where {R<:Number}

    return smolyak_derivative_full(weights,x,multi_index_full)

  end

  return goo

end

function smolyak_derivative_full(weights::Array{T,1},multi_index_full::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  function goo(x::Array{R,1}) where {R<:Number}

    return smolyak_derivative_full(weights,x,multi_index_full,domain)

  end

  return goo

end

function smolyak_derivative_full(weights::Array{T,1},multi_index_full::Array{S,2},pos::Array{S,1}) where {T<:AbstractFloat,S<:Integer}

  function goo(x::Array{R,1}) where {R<:Number}

    return smolyak_derivative_full(weights,x,multi_index_full,pos)

  end

  return goo

end

function smolyak_derivative_full(weights::Array{T,1},multi_index_full::Array{S,2},domain::Union{Array{T,1},Array{T,2}},pos::Array{S,1}) where {T<:AbstractFloat,S<:Integer}

  function goo(x::Array{R,1}) where {R<:Number}

    return smolyak_derivative_full(weights,x,multi_index_full,domain,pos)

  end

  return goo

end