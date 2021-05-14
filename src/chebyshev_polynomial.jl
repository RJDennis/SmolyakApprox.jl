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