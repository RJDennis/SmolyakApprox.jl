function smolyak_derivative(weights::Array{T,1},node::Array{T,1},multi_index::Array{S,2},pos::S) where {T<:AbstractFloat,S<:Integer}

  unique_multi_index = sort(unique(multi_index))
  unique_orders = m_i(unique_multi_index).-1

  # Here we construct the base polynomials

  base_polynomials = Array{Array{T,2},1}(undef,length(unique_orders))
  base_polynomial_derivatives = Array{Array{T,2},1}(undef,length(unique_orders))
  for i = 1:length(unique_orders)
    base_polynomials[i] = chebyshev_polynomial(unique_orders[i],node)
    base_polynomial_derivatives[i] = chebyshev_polynomial_derivative(unique_orders[i],node)
  end

  # Compute the unique polynomial terms from the base polynomials

  unique_base_polynomials = Array{Array{T,2},1}(undef,length(unique_orders))
  unique_base_polynomial_derivatives = Array{Array{T,2},1}(undef,length(unique_orders))
  for i = length(unique_orders):-1:2
    unique_base_polynomials[i] = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
    unique_base_polynomial_derivatives[i] = base_polynomial_derivatives[i][:,size(base_polynomial_derivatives[i-1],2)+1:end]
  end
  unique_base_polynomials[1] = base_polynomials[1]
  unique_base_polynomial_derivatives[1] = base_polynomial_derivatives[1]

  # Construct the first row of the interplation matrix

  polynomials = Array{T,1}(undef,length(weights))

  # Iterate over nodes, doing the above three steps at each iteration

  l = 1
  @inbounds for j = 1:size(multi_index,1)
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

  for i = 1:length(polynomials)
    evaluated_derivative += polynomials[i]*weights[i]
  end

  return evaluated_derivative

end

function smolyak_derivative(weights::Array{T,1},node::Array{T,1},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}},pos::S) where {T<:AbstractFloat,S<:Integer}

  nodes = copy(node)

  if size(domain,2) != length(node)
    error("domain is inconsistent with the number of dimensions")
  end

  d = length(node)
  for i = 1:d
    nodes[i] = normalize_node(node[i],domain[:,i])
  end

  evaluated_derivative = smolyak_derivative(weights,nodes,multi_index,pos)

  return evaluated_derivative

end

function smolyak_gradient(weights::Array{T,1},node::Array{T,1},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  d = length(node)
  gradient = Array{T,2}(undef,1,d)

  for i = 1:d
    gradient[i] = smolyak_derivative(weights,node,multi_index,i)
  end

  return gradient

end

function smolyak_gradient(weights::Array{T,1},node::Array{T,1},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  d = length(node)
  gradient = Array{T,2}(undef,1,d)

  for i = 1:d
    gradient[i] = smolyak_derivative(weights,node,multi_index,domain,i)
  end

  return gradient

end


#########################################################################

function smolyak_gradient(weights::Array{T,1},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  function goo(x::Array{T,1}) where {T <: AbstractFloat}

    return smolyak_gradient(weights,x,multi_index)

  end

  return goo

end

function smolyak_gradient(weights::Array{T,1},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  function goo(x::Array{T,1}) where {T <: AbstractFloat}

    return smolyak_gradient(weights,x,multi_index,domain)

  end

  return goo

end

function smolyak_derivative(weights::Array{T,1},multi_index::Array{S,2},pos::S) where {T<:AbstractFloat,S<:Integer}

  function goo(x::Array{T,1}) where {T <: AbstractFloat}

    return smolyak_derivative(weights,x,multi_index,pos)

  end

  return goo

end

function smolyak_derivative(weights::Array{T,1},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}},pos::S) where {T<:AbstractFloat,S<:Integer}

  function goo(x::Array{T,1}) where {T <: AbstractFloat}

    return smolyak_derivative(weights,x,multi_index,domain,pos)

  end

  return goo

end
