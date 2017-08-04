function smolyak_knn_points{T<:AbstractFloat,S<:Integer}(nodes::Array{T,2},x::Array{T,1},num_points::S)

  n = size(nodes,1)

  if num_points > n
    println("smolyak_knn_points: Too many points")
    throw(ErrorException())
  end

  dist_vector    = zeros(T,n)
  index_vector_k = Array{S}(num_points)
  dist_vector_k  = Array{T}(num_points)

  for i = 1:n
    dist_vector[i] = norm(x-nodes[i,:])^2
  end

  for i = 1:num_points
    (dist,index)       = findmin(dist_vector)
    dist_vector_k[i]   = dist
    index_vector_k[i]  = index
    dist_vector[index] = Inf
  end

  return index_vector_k, dist_vector_k

end

function smolyak_knn_evaluate{T<:AbstractFloat}(y::AbstractArray{T,1},nodes::Array{T,2},x::Array{T,1},k = 2*size(nodes,2)+1)

  (index_vector_k,dist_vector_k) = smolyak_knn_points(nodes,x,k)

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
