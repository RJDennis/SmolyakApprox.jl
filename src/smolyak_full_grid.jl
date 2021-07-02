function smolyak_grid_full(node_type::Function,d::S,mu::S) where {S<:Integer}

  T = typeof(1.0)

  multi_index        = generate_multi_index(d,mu)
  unique_multi_index = sort(unique(multi_index))
  unique_node_number = m_i(unique_multi_index)

  # Create base nodes to be used in the sparse grid

  base_nodes   = Array{Array{T,1},1}(undef,length(unique_node_number))
  base_weights = Array{Array{T,1},1}(undef,length(unique_node_number))
  @inbounds for i = 1:length(unique_node_number)
    base_nodes[i], base_weights[i] = node_type(unique_node_number[i])
  end

  # Select the relevant polynomials from the multi index

  ii = (sum(multi_index,dims=2) .>= max(d,mu+1)).*(sum(multi_index,dims=2) .<= d+mu)
  multi_index_full = zeros(S,sum(ii),size(multi_index,2))
  j = 1
  @inbounds for i = 1:size(multi_index,1)
    if ii[i] == true
      multi_index_full[j,:] = multi_index[i,:]
      j += 1
    end
  end

  # Construct the sparse grid from the nodes

  nodes = Array{T,2}(undef,determine_grid_size_full(multi_index_full))
  l = 1
  for j = 1:size(multi_index_full,1)
    new_nodes = base_nodes[multi_index_full[j,1]]  # Here new_nodes is a 1d array
    for i = 2:d
      new_nodes = combine_nodes(new_nodes,base_nodes[multi_index_full[j,i]])  # Here new_nodes becomes a 2d array
    end
    m = size(new_nodes,1)
    nodes[l:l+m-1,:] = new_nodes
    l += m
  end

  if d == 1
    nodes = nodes[:]
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

  @inbounds for i = 1:size(mi,1)
    @inbounds for j = 1:size(mi,2)
      if mi[i,j] == 1
        temp[i,j] = 1
      else
        temp[i,j] = 2^(mi[i,j]-1)+1
      end
    end
  end

  s = 0
  @inbounds for i = 1:size(mi,1)
    t = 1
    @inbounds for j = 1:size(mi,2)
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

function cheb_poly(order::S,x::T) where {S<:Integer, T<:AbstractFloat}

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

function prod_cjs(max_grid::Union{Array{T,1},Array{T,2}},min_grid::Union{Array{T,1},Array{T,2}},poly_grid::Array{T,2}) where {T<:AbstractFloat}

  cjs = ones(size(poly_grid))

  @inbounds for i = 1:size(poly_grid,1)
    @inbounds for j = 1:size(poly_grid,2)
      if poly_grid[i,j] == max_grid[j] || poly_grid[i,j] == min_grid[j]
        cjs[i,j] *= 2
      end
    end
  end

  return prod(cjs,dims=2)

end

function compute_scale_factor(multi_index::Array{S,1}) where {S<:Integer}

  scale_factor = 1.0

  @inbounds for i = 1:length(multi_index)
    if multi_index[i] > 1
      scale_factor *= 2.0/(m_i(multi_index[i])-1)
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

  @inbounds for i = 1:size(g_ind,1) # This loops over the number of polynomials
    ws = zeros(g_ind[i,2])
    poly_grid = grid[g_ind[i,1]:g_ind[i,1]+g_ind[i,2]-1,:]
    poly_y    = y_f[g_ind[i,1]:g_ind[i,1]+g_ind[i,2]-1]
	  if size(grid,1) == 1 # This is to accommodate the mu = 0 case
	    cjs_prod = ones(size(poly_grid,1))
	  else
      cjs_prod = prod_cjs(max_grid,min_grid,poly_grid)
    end
    @inbounds for l = 1:g_ind[i,2] # This loops over the weights
      ll = CartesianIndices(Tuple(m_i(multi_index[i,:])))[l]
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
  for i = 1:d
    grid[:,i] = normalize_node(grid[:,i],domain[:,i])
  end

  weights = smolyak_weights_full(y_f,grid,multi_index)

  return weights

end

function smolyak_evaluate_full(weights::Array{Array{T,1},1},point::Array{T,1},multi_index::Array{S,2}) where {S<:Integer, T<:AbstractFloat}

  d = size(multi_index,2)

  evaluated_polynomials = zeros(size(multi_index,1))

  @inbounds for i = 1:size(multi_index,1) # This loops over the number of polynomials
  	@inbounds for l = 1:length(weights[i])
      ll = CartesianIndices(Tuple(m_i(multi_index[i:i,:])))[l]
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

function smolyak_evaluate_full(weights::Array{Array{T,1},1},point::Array{T,1},multi_index::Array{S,2},domain::Array{T,2}) where {S<:Integer, T<:AbstractFloat}

  d = size(multi_index,2)
  for i = 1:d
    point[i] = normalize_node(point[i],domain[:,i])
  end

  estimate = smolyak_evaluate_full(weights,point,multi_index)

  return estimate

end

function deriv_cheb_poly(order::S,x::T) where {S<:Integer, T<:AbstractFloat}

  p0 = one(T)
  p1 = zero(T)
  p2 = zero(T)
  pd0 = zero(T)
  pd1 = zero(T)
  pd2 = zero(T)

  for i = 2:order+1
    if i == 2
	  p1, p0 = p0, x
	  pd1, pd0 = pd0, one(T)
    else
      p2, p1 = p1, p0
      pd2, pd1 = pd1, pd0
  	  p0  = 2*x*p1-p2
      pd0 = 2*p1+2*x*pd1-pd2 
    end
  end

  return pd0

end

function smolyak_derivative_full(weights::Array{Array{T,1},1},point::Array{T,1},multi_index::Array{S,2},pos::S) where {S<:Integer, T<:AbstractFloat}

  mi = sum(multi_index,dims=2)
  d  = size(multi_index,2)
  mu = maximum(mi)-d

  evaluated_polynomials = zeros(size(multi_index,1))

  for i = 1:size(multi_index,1) # This loops over the number of polynomials
    for l = 1:length(weights[i])
      ll = CartesianIndices(Tuple(m_i(multi_index[i:i,:])))[l]
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

function smolyak_derivative_full(weights::Array{Array{T,1},1},point::Array{T,1},multi_index::Array{S,2},domain::Array{T,2},pos::S) where {S<:Integer, T<:AbstractFloat}

  d = size(multi_index,2)
  for i = 1:d
    point[i] = normalize_node(point[i],domain[:,i])
  end

  evaluated_derivative = smolyak_derivative_full(weights,point,multi_index,pos)

  return evaluated_derivatives

end

function smolyak_gradient_full(weights::Array{Array{T,1},1},point::Array{T,1},multi_index::Array{S,2}) where {S<:Integer, T<:AbstractFloat}

  d = size(multi_index,2)
  
  gradient = Array{T,2}(undef,1,d)
  for i = 1:d
    gradient[i] = smolyak_derivative_full(weights,point,multi_index,i)
  end

  return gradient

end

function smolyak_gradient_full(weights::Array{Array{T,1},1},point::Array{T,1},multi_index::Array{S,2},domain::Array{T,2}) where {S<:Integer, T<:AbstractFloat}

  d = size(multi_index,2)
  
  gradient = Array{T,2}(undef,1,d)
  for i = 1:d
    gradient[i] = smolyak_derivative_full(weights,point,multi_index,domain,i)
  end

  return gradient

end

###########################################################################

function smolyak_evaluate_full(weights::Array{T,1},multi_index_full::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  function goo(x::Array{T,1}) where {T <: AbstractFloat}

    return smolyak_evaluate_full(weights,x,multi_index_full)

  end

  return goo

end

function smolyak_evaluate_full(weights::Array{T,1},multi_index_full::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  function goo(x::Array{T,1}) where {T <: AbstractFloat}

    return smolyak_evaluate_full(weights,x,multi_index_full,domain)

  end

  return goo

end

function smolyak_derivative_full(weights::Array{T,1},multi_index_full::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  function goo(x::Array{T,1}) where {T <: AbstractFloat}

    return smolyak_derivative_full(weights,x,multi_index_full)

  end

  return goo

end

function smolyak_derivative_full(weights::Array{T,1},multi_index_full::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  function goo(x::Array{T,1}) where {T <: AbstractFloat}

    return smolyak_derivative_full(weights,x,multi_index_full,domain)

  end

  return goo

end

function smolyak_derivative_full(weights::Array{T,1},multi_index_full::Array{S,2},pos::Array{S,1}) where {T<:AbstractFloat,S<:Integer}

  function goo(x::Array{T,1}) where {T <: AbstractFloat}

    return smolyak_derivative_full(weights,x,multi_index_full,pos)

  end

  return goo

end

function smolyak_derivative_full(weights::Array{T,1},multi_index_full::Array{S,2},domain::Union{Array{T,1},Array{T,2}},pos::Array{S,1}) where {T<:AbstractFloat,S<:Integer}

  function goo(x::Array{T,1}) where {T <: AbstractFloat}

    return smolyak_derivative_full(weights,x,multi_index_full,domain,pos)

  end

  return goo

end
