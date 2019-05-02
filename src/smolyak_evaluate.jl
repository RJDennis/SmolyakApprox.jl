function smolyak_evaluate(weights::Array{T,1},node::Array{T,1},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  unique_multi_index = sort(unique(multi_index))
  unique_orders = m_i(unique_multi_index).-1

  # Below we do the following things:

  #   Generate the polynomial terms for each order
  #   Generate the unique polynomial terms introduced at each higher order
  #   Combine the polynomial terms to construct the first row of the interpolation matrix
  #   Iterate over the nodes, doing the above three steps at each iteration, to compute all rows of the interpolation matrix

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

  new_polynomials = unique_base_polynomials[multi_index[1,1]][1,:]
  for i = 2:size(multi_index,2)
    new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[1,i]][i,:])
  end

  polynomials = copy(new_polynomials)

  # Iterate over nodes, doing the above three steps at each iteration

  for j = 2:size(multi_index,1)
    new_polynomials = unique_base_polynomials[multi_index[j,1]][1,:]
    for i = 2:size(multi_index,2)
      new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][i,:])
    end
    polynomials = [polynomials; new_polynomials]
  end

  estimate = zero(T)
  for i = 1:length(polynomials)
    estimate += polynomials[i]*weights[i]
  end

  return estimate

end

function smolyak_evaluate(weights::Array{T,1},node::Array{T,1},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}

  if size(domain,2) != size(node,2)
    error("domain is inconsistent with the number of dimensions")
  end

  node = copy(node)
  for j = 1:size(domain,2)
    if domain[1,j] == domain[2,j]
      node[j] = (domain[1,j]+domain[2,j])/2
    else
      node[j] = 2*(node[j]-domain[2,j])/(domain[1,j]-domain[2,j])-one(T)
    end
  end

  estimate = smolyak_evaluate(weights,node,multi_index)

  return estimate

end
