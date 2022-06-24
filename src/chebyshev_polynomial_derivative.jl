function chebyshev_polynomial_derivative(order::S,x::R) where {R<:Number,S<:Integer}

  polynomial    = Array{R,2}(undef,1,order+1)
  poly_deriv    = Array{R,2}(undef,1,order+1)
  polynomial[1] = one(R)
  poly_deriv[1] = zero(R)

  for i = 2:order+1
    if i == 2
      polynomial[i] = x
      poly_deriv[i] = one(R)
    else
      polynomial[i] = 2*x*polynomial[i-1]-polynomial[i-2]
      poly_deriv[i] = 2*polynomial[i-1]+2*x*poly_deriv[i-1]-poly_deriv[i-2]
    end
  end

  return poly_deriv

end

function chebyshev_polynomial_derivative(order::S,x::Array{R,1}) where {R<:Number,S<:Integer}

  poly_deriv = Array{R,2}(undef,length(x),order+1)

  for i in eachindex(x)
    poly_deriv[i,:] = chebyshev_polynomial_derivative(order,x[i])
  end

  return poly_deriv

end