function smolyak_grid_full(node_type::Function,d::S,mu::S) where {S<:Integer}

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

  # Select the relevant polynomials from the multi index

  ii = (sum(multi_index,dims=2) .>= max(d,mu+1)).*(sum(multi_index,dims=2) .<= d+mu)
  multi_index_full = zeros(Int64,sum(ii),size(multi_index,2))
  j = 1
  for i = 1:size(multi_index,1)
    if ii[i] == true
      multi_index_full[j:j,:] .= multi_index[i:i,:]
      j += 1
    end
  end

  # Construct the sparse grid from the nodes

  new_nodes = base_nodes[multi_index_full[1,1]]  # Here new_nodes is a 1d array
  for i = 2:d
    new_nodes = combine_nodes(new_nodes,base_nodes[multi_index_full[1,i]])  # Here new_nodes becomes a 2d array
  end

  nodes = copy(new_nodes)

  for j = 2:size(multi_index_full,1)
    new_nodes = base_nodes[multi_index_full[j,1]]
    for i = 2:d
      new_nodes = combine_nodes(new_nodes,base_nodes[multi_index_full[j,i]])
    end
    nodes = [nodes; new_nodes]
  end

  return nodes, multi_index_full

end

function smolyak_grid_full(node_type::Function,d::S,mu::S,domain::Array{T,2}) where {S<:Integer, T<:AbstractFloat}

  (nodes, multi_index_full) = smolyak_grid_full(node_type,d,mu)

  nodes = scale_nodes(nodes,domain)

  return nodes, multi_index_full

end

function master_index(multi_index::Array{S,2}) where {S<:Integer}

  temp_ind   = similar(multi_index)
  master_ind = zeros(S,size(multi_index,1),2)

  for i in eachindex(multi_index)
    if multi_index[i] == 1
      temp_ind[i] = 1
    else
      temp_ind[i] = m_i(multi_index[i])
    end
  end

  master_ind[1,1] = 1
  master_ind[1,2] = prod(temp_ind[1,:])

  for i = 2:size(master_ind,1)
    master_ind[i,1] = master_ind[i-1,1] + master_ind[i-1,2]
    master_ind[i,2] = prod(temp_ind[i,:])
  end

  return master_ind

end

function P(order::S,x::T) where {S<:Integer, T<:AbstractFloat}

  p  = 1.0
  p1 = 0.0
  p2 = 0.0

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

function prod_cjs(max_grid::Array{T,2},min_grid::Array{T,2},poly_grid::Array{T,2}) where {T<:AbstractFloat}

  cjs = ones(size(poly_grid))

  for i = 1:size(poly_grid,1)
    for j = 1:size(poly_grid,2)
      if poly_grid[i,j] == max_grid[j] || poly_grid[i,j] == min_grid[j]
        cjs[i,j] *= 2
      end
    end
  end

  return prod(cjs,dims=2)

end

function compute_scale_factor(poly_multi_index::Array{S,2}) where {S<:Integer}

  scale_factor = 1.0

  for i = 1:length(poly_multi_index)
    if poly_multi_index[i] > 1
      scale_factor *= 2.0/(m_i(poly_multi_index[i])-1)
    end
  end

  return scale_factor

end

function smolyak_weights_full(y_f::Array{T,1},full_grid::Array{T,2},multi_index_full::Array{S,2}) where {S<:Integer, T<:AbstractFloat}

  max_grid = maximum(full_grid,dims=1)
  min_grid = minimum(full_grid,dims=1)

  weights = Array{Array{T,1}}(undef,size(multi_index_full,1))
  g_ind = master_index(multi_index_full)

  for i = 1:size(g_ind,1) # This loops over the number of polynomials
    ws = zeros(g_ind[i,2])
    poly_grid = full_grid[g_ind[i,1]:g_ind[i,1]+g_ind[i,2]-1,:]
    poly_y    = y_f[g_ind[i,1]:g_ind[i,1]+g_ind[i,2]-1]
	  if size(full_grid,1) == 1 # This is to accommodate the mu = 0 case
	    cjs_prod = ones(size(poly_grid,1))
	  else
      cjs_prod = prod_cjs(max_grid,min_grid,poly_grid)
    end
    cls_prod = copy(cjs_prod)
    for l = 1:g_ind[i,2] # This loops over the weights
      ll = CartesianIndices(Tuple(m_i(multi_index_full[i:i,:])))[l]
      for j = 1:g_ind[i,2] # This loops over the nodes
        rhs_term = P(ll[1]-1,poly_grid[j,1])*poly_y[j]
        for k = 2:size(poly_grid,2) # This loops over the polynomial terms in the product
          rhs_term *= P(ll[k]-1,poly_grid[j,k])
        end
        ws[l] += rhs_term/cjs_prod[j]
      end
      scale_factor = compute_scale_factor(multi_index_full[i:i,:])
      ws[l] = scale_factor*(1/cjs_prod[l])*ws[l]
    end
    weights[i] = ws
  end

  return weights

end

function smolyak_weights_full(y_f::Array{T,1},full_grid::Array{T,2},multi_index_full::Array{S,2},domain::Array{T,2}) where {S<:Integer, T<:AbstractFloat}

  full_grid = copy(full_grid)
  for i = 1:size(full_grid,1)
    for j = 1:size(domain,2)
      if domain[1,j] == domain[2,j]
        full_grid[i,j] = (domain[1,j]+domain[2,j])/2
      else
        full_grid[i,j] = 2*(full_grid[i,j]-domain[2,j])/(domain[1,j]-domain[2,j])-one(T)
      end
    end
  end

  weights = smolyak_weights_full(y_f,full_grid,multi_index_full)

  return weights

end

function smolyak_evaluate_full(weights::Array{Array{T,1},1},node::Array{T,1},multi_index_full::Array{S,2}) where {S<:Integer, T<:AbstractFloat}

  mi = sum(multi_index_full,dims=2)
  d  = size(multi_index_full,2)
  mu = maximum(mi)-d

  evaluated_polynomials = zeros(size(multi_index_full,1))

  for i = 1:size(multi_index_full,1) # This loops over the number of polynomials
  	for l = 1:length(weights[i])
      ll = CartesianIndices(Tuple(m_i(multi_index_full[i:i,:])))[l]
      temp = weights[i][l]*P(ll[1]-1,node[1])
      for k = 2:d
  	    temp *= P(ll[k]-1,node[k])
      end
      evaluated_polynomials[i] += temp
    end
    evaluated_polynomials[i] *= (-1)^(d+mu-mi[i])*factorial(d-1)/(factorial(d+mu-mi[i])*factorial(-1-mu+mi[i]))
  end

  return sum(evaluated_polynomials)

end

function smolyak_evaluate_full(weights::Array{Array{T,1},1},node::Array{T,1},multi_index_full::Array{S,2},domain::Array{T,2}) where {S<:Integer, T<:AbstractFloat}

  node = copy(node)
  for j = 1:size(domain,2)
    if domain[1,j] == domain[2,j]
      node[j] = (domain[1,j]+domain[2,j])/2
    else
      node[j] = 2*(node[j]-domain[2,j])/(domain[1,j]-domain[2,j])-one(T)
    end
  end

  estimate = smolyak_evaluate_full(weights,node,multi_index_full)

  return estimate

end

function P_derivative(order::S,x::T) where {S<:Integer, T<:AbstractFloat}

  p  = 1.0
  p1 = 0.0
  p2 = 0.0
  p_deriv = 0.0

  for i = 2:order+1
    if i == 2
	  p1, p = p, x
	  p_deriv = 1.0
    else
      p2, p1 = p1, p
  	  p  = 2*x*p1-p2
	  p_deriv = ((i-1)*p1-(i-1)*x*p)/(1-x^2)
    end
  end

  return p_deriv

end

function smolyak_derivative_full(weights::Array{Array{T,1},1},node::Array{T,1},multi_index_full::Array{S,2},pos::Array{S,1}) where {S<:Integer, T<:AbstractFloat}

  mi = sum(multi_index_full,dims=2)
  d  = size(multi_index_full,2)
  mu = maximum(mi)-d

  evaluated_derivatives = zeros(1,length(pos))

  for m in pos

    evaluated_polynomials = zeros(size(multi_index_full,1))

    for i = 1:size(multi_index_full,1) # This loops over the number of polynomials
  	  for l = 1:length(weights[i])
        ll = CartesianIndices(Tuple(m_i(multi_index_full[i:i,:])))[l]
        temp = weights[i][l]*((m!==1)*P(ll[1]-1,node[1])+(m==1)*P_derivative(ll[1]-1,node[1]))
        for k = 2:d
  	      temp *= (m!==k)*P(ll[k]-1,node[k])+(m==k)*P_derivative(ll[k]-1,node[k])
	    end
	    evaluated_polynomials[i] += temp
      end
	  evaluated_polynomials[i] *= (-1)^(d+mu-mi[i])*factorial(d-1)/(factorial(d+mu-mi[i])*factorial(-1-mu+mi[i]))
    end

    evaluated_derivatives[m] = sum(evaluated_polynomials)

  end

  return evaluated_derivatives

end

function smolyak_derivative_full(weights::Array{Array{T,1},1},node::Array{T,1},multi_index_full::Array{S,2},domain::Array{T,2},pos::Array{S,1}) where {S<:Integer, T<:AbstractFloat}

  node = copy(node)
  for j = 1:size(domain,2)
    if domain[1,j] == domain[2,j]
      node[j] = (domain[1,j]+domain[2,j])/2
    else
      node[j] = 2*(node[j]-domain[2,j])/(domain[1,j]-domain[2,j])-one(T)
    end
  end

  evaluated_derivatives = smolyak_derivative_full(weights,node,multi_index,pos)

  return evaluated_derivatives

end
