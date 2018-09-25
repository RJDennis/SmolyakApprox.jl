function smolyak_derivative_finite_difference(weights::Array{T,1},node::Array{T,1},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

  m = length(node)
  e  = eps(T)^(1/3)*maximum(abs,[node;one(T)])
  dh = Matrix{T}(I,m,m)*e

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

function smolyak_derivative_finite_difference(weights::Array{T,1},node::Array{T,1},multi_index::Array{S,2},domain::Array{T,2}) where {T<:AbstractFloat,S<:Integer}

  m = length(node)
  e  = eps(T)^(1/3)*maximum(abs,[node;one(T)])
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
