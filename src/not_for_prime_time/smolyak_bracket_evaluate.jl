function smolyak_find_bracket(x::Array{T,1},nodes::Array{T,2}) where T<:AbstractFloat

  n = size(nodes,1)
  num_points = 2^size(nodes,2)

  dist_vector    = zeros(n)
  index_vector_k = Array{Int64}(num_points)
  dist_vector_k  = Array{Float64}(num_points)

  for i = 1:n
    dist_vector[i] = norm(x-nodes[i,:])
  end

  (dist,index)       = findmin(dist_vector)
  dist_vector_k[1]   = dist
  index_vector_k[1]  = index
  dist_vector[index] = Inf

  bool_vector = collect(x .< nodes[index,:])'
  counter = 1

  while counter < num_points

    (dist,index) = findmin(dist_vector)
    trial_vector = collect(x .< nodes[index,:])'

    inclusion_check = 0
    for i = 1:size(bool_vector,1)

      if sum(trial_vector .== bool_vector[i,:]) == size(bool_vector,2)
        inclusion_check += 1
      end

    end

    if inclusion_check == 0
      counter += 1
      dist_vector_k[counter]  = dist
      index_vector_k[counter] = index
      dist_vector[index]    = Inf
      bool_vector = [bool_vector; trial_vector]
    else
      dist_vector[index]    = Inf
    end

  end

  return index_vector_k,dist_vector_k

end

function smolyak_bracket_evaluate(y::AbstractArray{T,1},nodes::Array{T,2},x::Array{T,1}) where T<:AbstractFloat

  (index_vector_k,dist_vector_k) = smolyak_find_bracket(x,nodes)
  k = size(index_vector_k,1)

  weights = Array{T}(k)
  for i = 1:k
    weights[i] = 1.0/(dist_vector_k[i]+eps(T))
  end

  weights = weights./sum(weights)

  y_hat = zero(T)

  for i = 1:k
    y_hat += y[index_vector_k[i]]*weights[i]
  end

  return y_hat

end
