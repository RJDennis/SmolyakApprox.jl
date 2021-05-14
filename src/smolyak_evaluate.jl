function smolyak_evaluate(weights::Array{T,1},node::Array{T,1},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  unique_multi_index = sort(unique(multi_index))
  unique_orders = m_i(unique_multi_index).-1

  # Below we do the following things:

  #   Generate the polynomial terms for each order
  #   Generate the unique polynomial terms introduced at each higher order
  #   Combine the polynomial terms to construct the Smolyak polynomial

  # Here we construct the base polynomials

  base_polynomials = Array{Array{T,2}}(undef,length(unique_orders))
  for i = 1:length(unique_orders)
    base_polynomials[i] = chebyshev_polynomial(unique_orders[i],node)
  end

  # Compute the unique polynomial terms from the base polynomials

  unique_base_polynomials = Array{Array{T,2}}(undef,length(unique_orders))
  for i = length(unique_orders):-1:2
    unique_base_polynomials[i] = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
  end
  unique_base_polynomials[1] = base_polynomials[1]

  # Construct the first row of the interplation matrix

  polynomials = Array{T,1}(undef,length(weights))

  # Iterate over nodes, doing the above three steps at each iteration

  l = 1
  @inbounds for j = 1:size(multi_index,1)
    new_polynomials = unique_base_polynomials[multi_index[j,1]][1,:]
    for i = 2:size(multi_index,2)
      new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][i,:])
    end
    m = length(new_polynomials)
    polynomials[l:l+m-1] = new_polynomials
    l += m
  end

  estimate = zero(T)
  for i = 1:length(polynomials)
    estimate += polynomials[i]*weights[i]
  end

  return estimate

end

function smolyak_evaluate(weights::Array{T,1},node::Array{T,1},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

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

###########################################################################

function smolyak_evaluate(weights::Array{T,1},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  function goo(x::Array{T,1}) where {T <: AbstractFloat}

    return smolyak_evaluate(weights,x,multi_index)

  end

  return goo

end

function smolyak_evaluate(weights::Array{T,1},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  function goo(x::Array{T,1}) where {T <: AbstractFloat}

    return smolyak_evaluate(weights,x,multi_index,domain)

  end

  return goo

end
