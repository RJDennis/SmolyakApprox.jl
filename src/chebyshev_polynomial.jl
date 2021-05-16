function chebyshev_polynomial(order::S,x::Array{T,1}) where {T<:AbstractFloat,S<:Integer}

  polynomial      = Array{T}(undef,length(x),order+1)
  polynomial[:,1] = ones(T,length(x))
  
  for i = 2:order+1
      for j = 1:length(x)
      if i == 2
        polynomial[j,i] = x[j]
      else
        polynomial[j,i] = 2*x[j]*polynomial[j,i-1]-polynomial[j,i-2]
      end
    end
  end
  
  return polynomial
  
end

function chebyshev_polynomial(order::S,x::T) where {T<:AbstractFloat,S<:Integer}

  polynomial = ones(T,1,order+1)
  
  for i = 2:order+1
    if i == 2
      polynomial[1,i] = x
    else
      polynomial[1,i] = 2*x*polynomial[1,i-1]-polynomial[1,i-2]
    end
  end
  
  return polynomial
  
end
