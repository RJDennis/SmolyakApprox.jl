function chebyshev_polynomial(order::S,x::Array{R,1}) where {R<:Number,S<:Integer}

  polynomial      = Array{R}(undef,length(x),order+1)
  polynomial[:,1] = ones(R,length(x))

  for i = 2:order+1
      for j in eachindex(x)
      if i == 2
        polynomial[j,i] = x[j]
      else
        polynomial[j,i] = 2*x[j]*polynomial[j,i-1]-polynomial[j,i-2]
      end
    end
  end

  return polynomial

end

function chebyshev_polynomial(order::S,x::R) where {R<:Number,S<:Integer}

  polynomial = ones(R,1,order+1)

  for i = 2:order+1
    if i == 2
      polynomial[1,i] = x
    else
      polynomial[1,i] = 2*x*polynomial[1,i-1]-polynomial[1,i-2]
    end
  end

  return polynomial

end