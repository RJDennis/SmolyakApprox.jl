function smolyak_derivative(weights::Array{T,1},node::Array{T,1},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  unique_multi_index = sort(unique(multi_index))
  unique_orders = m_i(unique_multi_index)-1

#  m_node_number = m_i(multi_index)
#  multi_orders  = m_node_number-1
#  unique_orders = sort(unique(multi_orders))

  # Below we do the following things:

  #   Generate the polynomial terms for each order
  #   Generate the unique polynomial terms introduced at each higher order
  #   Combine the polynomial terms to construct the first row of the interpolation matrix
  #   Iterate over the nodes, doing the above three steps at each iteration, to compute all rows of the interpolation matrix

  # Here we construct the base polynomials

  base_polynomials = Array{Array{T,2}}(length(unique_orders))
  base_polynomial_derivatives = Array{Array{T,2}}(length(unique_orders))
  for i = 1:length(unique_orders)
    base_polynomials[i] = chebyshev_polynomial(unique_orders[i],node)
    base_polynomial_derivatives[i] = chebyshev_polynomial_derivative(unique_orders[i],node)
  end

  # Compute the unique polynomial terms from the base polynomials

  unique_base_polynomials = Array{Array{T,2}}(length(unique_orders))
  unique_base_polynomial_derivatives = Array{Array{T,2}}(length(unique_orders))
  for i = length(unique_orders):-1:2
    unique_base_polynomials[i] = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
    unique_base_polynomial_derivatives[i] = base_polynomial_derivatives[i][:,size(base_polynomial_derivatives[i-1],2)+1:end]
  end
  unique_base_polynomials[1] = base_polynomials[1]
  unique_base_polynomial_derivatives[1] = base_polynomial_derivatives[1]

  # Construct the first row of the interplation matrix

  evaluated_derivative = zeros(1,length(node))

  new_polynomials = Array{Array{T,1}}(length(node))
  polynomials     = Array{Array{T,1}}(length(node))

  for k = 1:length(node)
    new_polynomials[k] = (k!==1)*unique_base_polynomials[multi_index[1,1]][1,:]+(k==1)*unique_base_polynomial_derivatives[multi_index[1,1]][1,:]
    for i = 2:size(multi_index,2)
      new_polynomials[k] = kron(new_polynomials[k],(k!==i)*unique_base_polynomials[multi_index[1,i]][i,:]+(k==i)*unique_base_polynomial_derivatives[multi_index[1,i]][i,:])
    end

    polynomials[k] = copy(new_polynomials[k])

    # Iterate over nodes, doing the above three steps at each iteration

    for j = 2:size(multi_index,1)

      new_polynomials[k] = (k!==1)*unique_base_polynomials[multi_index[j,1]][1,:]+(k==1)*unique_base_polynomial_derivatives[multi_index[j,1]][1,:]
      for i = 2:size(multi_index,2)
        new_polynomials[k] = kron(new_polynomials[k],(k!==i)*unique_base_polynomials[multi_index[j,i]][i,:]+(k==i)*unique_base_polynomial_derivatives[multi_index[j,i]][i,:])
      end
      polynomials[k] = [polynomials[k]; new_polynomials[k]]

    end

    println(polynomials)

    for i = 1:length(polynomials[k])
      evaluated_derivative[k] += polynomials[k][i]*weights[i]
    end

  end

  return evaluated_derivative

end

function smolyak_derivative(weights::Array{T,1},node::Array{T,1},multi_index::Array{S,2},pos::S) where {T<:AbstractFloat,S<:Integer}

  unique_multi_index = sort(unique(multi_index))
  unique_orders = m_i(unique_multi_index)-1

  # Here we construct the base polynomials

  base_polynomials = Array{Array{T,2}}(length(unique_orders))
  base_polynomial_derivatives = Array{Array{T,2}}(length(unique_orders))
  for i = 1:length(unique_orders)
    base_polynomials[i] = chebyshev_polynomial(unique_orders[i],node)
    base_polynomial_derivatives[i] = chebyshev_polynomial_derivative(unique_orders[i],node)
  end

  # Compute the unique polynomial terms from the base polynomials

  unique_base_polynomials = Array{Array{T,2}}(length(unique_orders))
  unique_base_polynomial_derivatives = Array{Array{T,2}}(length(unique_orders))
  for i = length(unique_orders):-1:2
    unique_base_polynomials[i] = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
    unique_base_polynomial_derivatives[i] = base_polynomial_derivatives[i][:,size(base_polynomial_derivatives[i-1],2)+1:end]
  end
  unique_base_polynomials[1] = base_polynomials[1]
  unique_base_polynomial_derivatives[1] = base_polynomial_derivatives[1]

  # Construct the first row of the interplation matrix

  evaluated_derivative = 0.0

  new_polynomials = Array{Array{T,1}}(1)
  polynomials     = Array{Array{T,1}}(1)

  new_polynomials = (pos!==1)*unique_base_polynomials[multi_index[1,1]][1,:]+(pos==1)*unique_base_polynomial_derivatives[multi_index[1,1]][1,:]
  for i = 2:size(multi_index,2)
    new_polynomials = kron(new_polynomials,(pos!==i)*unique_base_polynomials[multi_index[1,i]][i,:]+(pos==i)*unique_base_polynomial_derivatives[multi_index[1,i]][i,:])
  end

  polynomials = copy(new_polynomials)

  # Iterate over nodes, doing the above three steps at each iteration

  for j = 2:size(multi_index,1)

    new_polynomials = (pos!==1)*unique_base_polynomials[multi_index[j,1]][1,:]+(pos==1)*unique_base_polynomial_derivatives[multi_index[j,1]][1,:]
    for i = 2:size(multi_index,2)
      new_polynomials = kron(new_polynomials,(pos!==i)*unique_base_polynomials[multi_index[j,i]][i,:]+(pos==i)*unique_base_polynomial_derivatives[multi_index[j,i]][i,:])
    end

    polynomials = [polynomials; new_polynomials]
    println(polynomials)

    for i = 1:length(polynomials)
      evaluated_derivative += polynomials[i]*weights[i]
    end

  end

  return evaluated_derivative

end

function smolyak_derivative(weights::Array{T,1},node::Array{T,1},multi_index::Array{S,2},domain::Array{T,2}) where {T<:AbstractFloat,S<:Integer}

  node = copy(node)
  for j = 1:size(domain,2)
    if domain[1,j] == domain[2,j]
      node[j] = (domain[1,j]+domain[2,j])/2
    else
      node[j] = 2*(node[j]-domain[2,j])/(domain[1,j]-domain[2,j])-one(T)
    end
  end

  evaluated_derivative = smolyak_derivative(weights,node,multi_index)

  return evaluated_derivative

end

function smolyak_derivative(weights::Array{T,1},node::Array{T,1},multi_index::Array{S,2},domain::Array{T,2},pos::S) where {T<:AbstractFloat,S<:Integer}

  node = copy(node)
  for j = 1:size(domain,2)
    if domain[1,j] == domain[2,j]
      node[j] = (domain[1,j]+domain[2,j])/2
    else
      node[j] = 2*(node[j]-domain[2,j])/(domain[1,j]-domain[2,j])-one(T)
    end
  end

  evaluated_derivative = smolyak_derivative(weights,node,multi_index,pos)

  return evaluated_derivative

end
