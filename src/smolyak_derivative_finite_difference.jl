function smolyak_derivative_finite_difference(weights::Array{T,1},node::Array{R,1},multi_index::Array{S,2}) where {T<:AbstractFloat,R<:Number,S<:Integer}

  m = length(node)
  e  = eps(T)^(1/3)*maximum(abs,[node;one(R)])
  dh = Matrix{R}(I,m,m)*e

  evaluated_derivative = zeros(1,m)

  for i = 1:m
    f1 = smolyak_evaluate(weights,node+2*dh[:,i],multi_index)
    f2 = smolyak_evaluate(weights,node+dh[:,i],multi_index)
    f3 = smolyak_evaluate(weights,node-dh[:,i],multi_index)
    f4 = smolyak_evaluate(weights,node-2*dh[:,i],multi_index)
    evaluated_derivative[i] = (-f1+8*f2-8*f3+f4)/(12*e)
  end

  return evaluated_derivative

end

function smolyak_derivative_finite_difference(weights::Array{T,1},node::Array{R,1},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,R<:Number,S<:Integer}

  if size(domain,2) != length(node)
    error("domain is inconsistent with the number of dimensions")
  end

  m = length(node)
  e  = eps(T)^(1/3)*maximum(abs,[node;one(R)])
  dh = eye(m)*e

  evaluated_derivative = zeros(1,m)

  for i = 1:m
    f1 = smolyak_evaluate(weights,node+2*dh[:,i],multi_index,domain)
    f2 = smolyak_evaluate(weights,node+dh[:,i],multi_index,domain)
    f3 = smolyak_evaluate(weights,node-dh[:,i],multi_index,domain)
    f4 = smolyak_evaluate(weights,node-2*dh[:,i],multi_index,domain)
    evaluated_derivative[i] = (2.0/(domain[1,i]-domain[2,i]))*(-f1+8*f2-8*f3+f4)/(12*e)
  end

  return evaluated_derivative

end