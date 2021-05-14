function chebyshev_polynomial_derivative(order::S,x::T) where {T<:AbstractFloat,S<:Integer}

  polynomial    = Array{T}(undef,1,order+1)
  poly_deriv    = Array{T}(undef,1,order+1)
  polynomial[1] = one(T)
  poly_deriv[1] = zero(T)

  for i = 2:order+1
    if i == 2
      polynomial[i] = x
      poly_deriv[i] = one(T)
    else
      polynomial[i] = 2*x*polynomial[i-1]-polynomial[i-2]
      poly_deriv[i] = 2*polynomial[i-1]+2*x*poly_deriv[i-1]-poly_deriv[i-2] 
    end
  end

  return poly_deriv

end

function chebyshev_polynomial_derivative(order::S,x::Array{T,1}) where {T<:AbstractFloat,S<:Integer}

  poly_deriv = Array{T,2}(undef,length(x),order+1)

  for i = 1:length(x)
    poly_deriv[i,:] = chebyshev_polynomial_derivative(order,x[i])
  end

  return poly_deriv

end
