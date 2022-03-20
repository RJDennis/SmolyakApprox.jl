# This function relates to the ansiotropic case

function generate_multi_index(d::S,mu::Array{S,1}) where {S<:Integer}

  nt = num_terms(mu,d)
  multi_index = Array{S,2}(undef,nt,d)
  multi_index[1,:] = ones(S,1,d)

  max_mu = maximum(mu)

  w = Tuple(repeat([max_mu+1],inner = d))
  pos = 0
  @inbounds for i = 2:(max_mu+1)^d
    candidate_index = Tuple(CartesianIndices(w)[i])
    if sum(candidate_index) <= d+max_mu && sum(candidate_index .<= mu.+1) == d
        pos += 1
      if pos > nt # handles the case where nt is under-estimated
        multi_index = [multi_index; collect(candidate_index)']
      else
        multi_index[pos,:] .= candidate_index
      end
    end
  end

  if pos < nt # handles case where nt is over-estimated
    multi_index = multi_index[1:pos,:]
  end

  return multi_index

end

# The function below relates to the isotropic case

function generate_multi_index(d::S,mu::S) where {S<:Integer}

  if d < 1
    error("d must be positive")
  end

  if mu < 0
    error("mu must be non-negative")
  end

  if d == 1
    multi_index = zeros(mu+1)
    for i = 1:mu+1
      multi_index[i] = i
    end
    return multi_index
  else
    multi_index_base = generate_multi_index(d-1,mu)
    N = size(multi_index_base,1)
    multi_index = zeros(S,N*(mu+1),d)
    pos = 0
    @inbounds @views for j = 1:N
      for i = 1:mu+1
        if sum(multi_index_base[j,:]) + i <= d+mu
          pos += 1
          multi_index[pos,2:d] .= multi_index_base[j,:]
          multi_index[pos,1] = i
        end
      end
    end
    return multi_index[1:pos,:]
  end
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

  max_mu = maximum(order)
  nt = num_terms(max_mu,d) # Deliberate over-estimate of the number of terms
    
  return nt
  
end
