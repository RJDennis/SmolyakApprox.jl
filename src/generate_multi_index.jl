# This function relates to the ansiotropic case

function generate_multi_index(d::S,mu::Array{S,1}) where {S<:Integer}

  nt = num_terms(mu,d)
  multi_index = Array{S,2}(undef,nt,d)
  multi_index[1,:] = ones(S,1,d)

  max_mu = maximum(mu)

  w = Tuple(repeat([max_mu+1],inner = d))
  j = 1
  @inbounds for i = 2:(max_mu+1)^d
    candidate_index = Tuple(CartesianIndices(w)[i])
    if sum(candidate_index) <= d+max_mu && sum(candidate_index .<= mu.+1) == d
      if j < nt
        j += 1
        multi_index[j,:] = [candidate_index...]
      else # handles case where nt is under-estimated
        multi_index = [multi_index; collect(candidate_index)']
      end
    end
  end
  if j < nt # handles case where nt is over-estimated
    multi_index = multi_index[1:j,:]
  end

  return multi_index

end

# The function below relates to the isotropic case

function generate_multi_index(d::S,mu::S) where {S<:Integer}

  nt = num_terms(mu,d)
  multi_index      = Array{S,2}(undef,nt,d)
  multi_index[1,:] = ones(S,1,d)

  w = Tuple(repeat([mu+1],inner = d))
  j = 1
  @inbounds for i = 2:(mu+1)^d
    candidate_index = Tuple(CartesianIndices(w)[i])
    if sum(candidate_index) <= d+mu
      j += 1
      multi_index[j,:] = [candidate_index...]
    end
  end

  return multi_index

end

# Following function computes the number of terms in the multi-index for the 
# isotropic case (it also computes the number of terms in a complete 
# polynominal based on the order and the number of dimensions. 

function num_terms(order::S,d::S) where {S <: Integer}

  if d == 1
    return order+1
  else
    return div(num_terms(order,d-1)*(order+d),d)
  end

end

# The following function is a poor approximation to the number of terms in 
# the multi-index for the ansiotropic case.

function num_terms(order::Array{S,1},d::S) where {S <: Integer}

  nt = div(prod(order.+2),2^(d-1)) # Deliberate over-estimate of the number of terms
  
  return nt

end
