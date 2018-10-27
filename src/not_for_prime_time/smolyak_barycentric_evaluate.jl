function smolyak_barycentric_points(nodes::Array{T,2},x::Array{T,1},num_points::S) where {T<:AbstractFloat,S<:Integer}

  n = size(nodes,1)

  if num_points > n
    println("smolyak_barycentric_points: Too many points")
    throw(ErrorException())
  end

  dist_vector    = zeros(T,n)
  index_vector_k = Array{S}(num_points)

  for i = 1:n
    dist_vector[i] = norm(x-nodes[i,:])
  end

  for i = 1:num_points
    (dist,index)       = findmin(dist_vector)
    index_vector_k[i]  = index
    dist_vector[index] = Inf
  end

  return index_vector_k

end

function smolyak_barycentric_evaluate(y::AbstractArray{T,1},nodes::Array{T,2},x::Array{T,1}) where T<:AbstractFloat

  k = size(nodes,2)+1
  index_vector_k = smolyak_barycentric_points(nodes,x,k)

  l = vcat(nodes[index_vector_k,:]',ones(1,k))
  println(l)
  weights = l\[x;1.0]

  y_hat = zero(T)

  for i = 1:k
    y_hat += y[index_vector_k[i]]*weights[i]
  end

  return y_hat

end
